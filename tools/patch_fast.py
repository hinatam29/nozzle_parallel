# -*- coding: utf-8 -*-
"""
計算の高速化パッチ（2026-10-22）
  python3 patch_fast.py <計算フォルダ>

変えるのは2つのファイルだけで、物理（式・条件）は変えない。
  intraction.f : 液滴どうしのクーロン力
      - 直接・左の鏡像・右の鏡像（並列のみ）の3つの二重ループを1つにまとめる
      - 飛行中（pout=0）の液滴だけで計算する（付着した液滴は電荷0で力も使われないため）
      - 内側のループをベクトル化（omp simd）
  diffusion1.f : 気体の拡散
      - 格子のループを OpenMP で並列化
元のファイルは intraction.f.orig / diffusion1.f.orig として残す。
2回目以降の実行では何もしない。
Makefile に -march=native を足す（そのパソコンの CPU に合わせた命令を使う）。
"""
import os, re, sys, shutil

D = sys.argv[1]

TEMPLATE = r'''      subroutine intraction
      include 'subcom.inc'
c     =================================================================
c     液滴どうしのクーロン力（高速化版、2026-10-22、patch_fast.py）
c     物理は旧版と同じ:
c       (1) 直接の力
c       (2) 左の鏡面 x=0 の鏡像の力（自分自身の鏡像は含めない＝旧版と同じ）
RIGHTCOMMENT
c     変えたところ（結果は旧版と一致することを確認済み）:
c       - 二重ループを1つにまとめた
c       - 飛行中の液滴（pout=0）だけを集めて計算する
c         （付着などした液滴は電荷 0 として扱われ、力も使われないため）
c       - 内側のループをベクトル化（omp simd）した
c     元のファイルは intraction.f.orig
c     =================================================================
      real xa(100000), ya(100000), da(100000), qa(100000)
      integer ia(100000)
      real fxs, fys, fx1, fy1, fx2, fy2
      real rx, ry, r, r3, dd, w, xk, yk, dk

      eps0 = 8.8542e-12
      fs = 42.9e-3
      pi = 4.*atan(1.)
XLENLINE

c     --- 飛行中の液滴を集める ---
      na = 0
      do 2 k=1, ip
         ftxt(k) = 0.
         ftyt(k) = 0.
         if( pout(k) .eq. 0. )then
            na = na + 1
            ia(na) = k
            xa(na) = xp(k)
            ya(na) = yp(k)
            da(na) = dp(k)
            qa(na) = pi*sqrt( 8.*eps0*fs*dp(k)*dp(k)*dp(k) )
         end if
 2    continue

c$omp parallel do schedule(static)
c$omp& private(kk,ll,xk,yk,dk,fxs,fys,fx1,fy1,fx2,fy2,
c$omp&         rx,ry,r,r3,dd,w)
      do 10 kk=1, na
         xk = xa(kk)
         yk = ya(kk)
         dk = da(kk)
         fxs = 0.
         fys = 0.
         fx1 = 0.
         fy1 = 0.
         fx2 = 0.
         fy2 = 0.
c$omp simd reduction(+:fxs,fys,fx1,fy1,fx2,fy2)
c$omp& private(rx,ry,r,r3,dd,w)
         do 20 ll=1, na
            dd = dk + da(ll)
            ry = yk - ya(ll)
c           自分自身（ll=kk）は重み 0
            w = merge( 0., 1., ll .eq. kk )
c           (1) 直接
            rx = xk - xa(ll)
            r  = sqrt( max( rx*rx + ry*ry, dd*dd ) )
            r3 = w*qa(ll)/(r*r*r)
            fxs = fxs + r3*rx
            fys = fys + r3*ry
c           (2) 左の鏡面 x=0 の鏡像
            rx = xk + xa(ll)
            r  = sqrt( max( rx*rx + ry*ry, dd*dd ) )
            r3 = w*qa(ll)/(r*r*r)
            fx1 = fx1 + r3*rx
            fy1 = fy1 + r3*ry
RIGHTBLOCK
 20      continue
         ftxt(ia(kk)) = fxs + fx1 + fx2
         ftyt(ia(kk)) = fys + fy1 + fy2
 10   continue
c$omp end parallel do

      return
      end
'''

RIGHT = r'''c           (3) 右の鏡面 x=xlen の鏡像
            rx = xk + xa(ll) - 2.*xlen
            r  = sqrt( max( rx*rx + ry*ry, dd*dd ) )
            r3 = w*qa(ll)/(r*r*r)
            fx2 = fx2 + r3*rx
            fy2 = fy2 + r3*ry'''


def patch_intraction():
    p = os.path.join(D, 'intraction.f')
    s = open(p, encoding='utf-8', errors='replace').read()
    if 'patch_fast.py' in s:
        print('  intraction.f: すでに高速化済み')
        return
    # 右の鏡像のループが有効（コメントでない）か
    right = re.search(r'^\s+rx2 = xp\(k\) \+ xp\(l\) - 2\.\*xlen', s, re.M) is not None
    m = re.search(r'^      xlen = .*$', s, re.M)
    xlenline = m.group(0) if (right and m) else 'c     （右の鏡面なし）'
    if right and not m:
        sys.exit('intraction.f に xlen の行が見つかりません')
    t = TEMPLATE.replace('XLENLINE', xlenline)
    t = t.replace('RIGHTBLOCK', RIGHT if right else 'c           （右の鏡面なし: 1本ノズル・円筒）')
    t = t.replace('RIGHTCOMMENT', 'c       (3) 右の鏡面 x=xlen の鏡像の力（自分自身の鏡像は含めない）'
                  if right else 'c       （右の鏡面はなし）')
    shutil.copy(p, p + '.orig')
    open(p, 'w', encoding='utf-8').write(t)
    print('  intraction.f: 高速化（右の鏡面 %s）' % ('あり, ' + xlenline.strip() if right else 'なし'))


def patch_diffusion1():
    p = os.path.join(D, 'diffusion1.f')
    s = open(p, encoding='utf-8', errors='replace').read()
    if 'patch_fast.py' in s:
        print('  diffusion1.f: すでに高速化済み')
        return
    i0 = s.index('       do 100 i=2, nx-2')
    i1 = s.index(' 100   continue') + len(' 100   continue')
    blk = s[i0:i1]
    assert blk.count('       do 100 j=2, ny-2') == 1
    assert blk.count(' 110   continue') == 1
    blk = blk.replace('       do 100 j=2, ny-2', '       do 105 j=2, ny-2')
    blk = re.sub(r' 110   continue\s*\n\s*\n?\s*100   continue',
                 ' 110   continue\n 105   continue\n 100   continue', blk)
    assert ' 105   continue' in blk
    head = ('c     --- 高速化(2026-10-22, patch_fast.py): 格子ごとに独立なので\n'
            'c     OpenMP で並列化（計算の中身は旧版と同じ）---\n'
            'c$omp parallel do private(i,j,cg0,cg2,cg4,cg6,cg8,cgn0)\n')
    blk = head + blk + '\nc$omp end parallel do'
    shutil.copy(p, p + '.orig')
    open(p, 'w', encoding='utf-8').write(s[:i0] + blk + s[i1:])
    print('  diffusion1.f: 並列化')


def patch_makefile():
    p = os.path.join(D, 'Makefile')
    s = open(p, encoding='utf-8', errors='replace').read()
    if '-march=native' in s:
        return
    s = s.replace('gfortran -c -O3 -fopenmp', 'gfortran -c -O3 -march=native -fopenmp')
    open(p, 'w', encoding='utf-8').write(s)
    print('  Makefile: -march=native を追加')


print(D)
patch_intraction()
patch_diffusion1()
patch_makefile()
print('完了。make をやり直してください（make clean がなければ rm -f *.o; make）')
