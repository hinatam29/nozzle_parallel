# -*- coding: utf-8 -*-
"""
並列ノズル v3 版を作るパッチ（2026-10-22）
  python3 patch_v3.py <fieldFIX2 の計算フォルダのコピー>

物理・数値の変更（10/22 ゼミの指摘への対応、第1段階）:
  1. 液滴の出し方（initial.f, main.f）
     1回に 5 個を、ノズルの内径 0.1 mm の幅に横に並べて同時に出す
     （x = 0, ±0.025, ±0.05 mm、下向き 3 m/s）。出す回数は今まで通り
     7 ステップに1回／ノズル。鏡面（ノズルの軸）を使っているので、計算する
     液滴は軸の上 1 個と、片側の 2 個（x = 0.025, 0.05 mm）の計 3 個で、
     反対側の 2 個は鏡像として入る。→ 実際の液滴の数（流量）は 5 倍
     液滴番号 k の並び: (k-1) を 6 で割った余り 0,1,2 = 左のノズル
     （軸から 0, 0.025, 0.05 mm）、3,4,5 = 右のノズル（同じ）
     集計するときの重み: 軸の上で出た液滴 1、横に出た液滴 2
  2. 液滴どうしの力（intraction.f）: 自分自身の鏡像の力も入れる＋高速化
  3. 液滴の位置の電場（eforce.f）: 4つの格子点からの補間を、距離の逆数の
     重みから「双一次補間」に変更。軸の上（x = 0）では横向きの電場が
     ちょうど 0 になる（旧版は隣の格子点の外向きの電場が混ざっていた）
  4. 電位（potential.f, ej.f）: 鏡面をちょうど格子点の列（x = 0, x = 間隔）に
     置く（旧版は φ1=φ2 で、鏡面が半マスずれていた）。横向きの電場は鏡面で 0。
     解き直しを変化がなくなるまで繰り返す
  5. 形（geo.f）: ノズルと電極の端を一番近い格子点に合わせる
     （旧版はノズル先端 98.8 mm、電極 91.8〜92.6 mm と 0.2 mm 下にずれていた
      → ノズル先端 99.0 mm、電極 92.0〜93.0 mm）
  6. 気体（diffusion1.f, diffusion_gas.f）: 鏡面の列も拡散を計算し、
     気体の減り方＝液滴の吸収量になるようにした
"""
import os, re, sys, shutil

D = sys.argv[1]
SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'v3src')


def rd(f):
    return open(os.path.join(D, f), encoding='utf-8', errors='replace').read()


def wr(f, s):
    p = os.path.join(D, f)
    if not os.path.exists(p + '.fix2'):
        shutil.copy(p, p + '.fix2')
    open(p, 'w', encoding='utf-8').write(s)


def rep(s, old, new, f):
    n = s.count(old)
    if n != 1:
        sys.exit('%s: 置き換える文字列が %d 回見つかりました: %r' % (f, n, old[:60]))
    return s.replace(old, new)


if 'patch_v3' in rd('intraction.f'):
    sys.exit('すでに v3 になっています')

xlenline = re.search(r'^      xlen = .*$', rd('initial.f'), re.M).group(0)
print(D, ' ', xlenline.strip())

# ---- 1. 液滴の出し方 ----
s = rd('initial.f')
s = rep(s, """      if( mod(k,2) .eq. 0 )then
         xp(k) = xlen
      else
         xp(k) = 0.
      end if""", """c     --- v3 (2026-10-22): 1回に 5 個を内径 0.1 mm の幅に並べて同時に出す ---
c     (k-1) mod 6 = 0,1,2 : 左のノズル（軸から 0, 0.025, 0.05 mm）
c                   3,4,5 : 右のノズル（軸から 0, 0.025, 0.05 mm）
c     反対側の 0.025, 0.05 mm の2個は鏡像として入る（計 5 個）
      kk6 = mod(k-1,6)
      if( kk6 .ge. 3 )then
         xp(k) = xlen - 0.025e-3*real(kk6-3)
      else
         xp(k) = 0.025e-3*real(kk6)
      end if""", 'initial.f')
wr('initial.f', s)

s = rd('main.f')
s = rep(s, """      if( ip .lt. ipmx-1 )then
         if( iph .eq. ns1 ) ip = ip + 1
         if( iph .eq. ns2 ) ip = ip + 1
      end if""", """c     --- v3 (2026-10-22): 1回に 3 個（＋鏡像 2 個 ＝ 5 個）を同時に出す ---
      if( ip .lt. ipmx-6 )then
         if( iph .eq. ns1 ) ip = ip + 3
         if( iph .eq. ns2 ) ip = ip + 3
      end if""", 'main.f')
# 最初は左のノズルの 3 個が出た状態から始める（液滴番号の並びをそろえるため）
s = rep(s, """      ip = 1.
""", """c     v3: 最初は左のノズルの 3 個（＋鏡像 2 個）が出た状態から始める
      ip = 3
""", 'main.f')
wr('main.f', s)

# ---- 2. 液滴どうしの力 ----
t = open(os.path.join(SRC, 'intraction_v3.f'), encoding='utf-8').read()
t = t.replace('XLENLINE', xlenline)
t = t.replace('c     =================================================================\n      real xa',
              'c     （patch_v3）\nc     =================================================================\n      real xa', 1)
wr('intraction.f', t)

# ---- 3. 電場の補間 ----
s = rd('eforce.f')
i0 = s.index('          Ernex(k) = ( Er(ijk)*l1(k)')
i1 = s.index('/ ( l1(k)+l2(k)+l3(k)+l4(k) ) ', i0 + 10)
i1 = s.index('\n', i1) + 1
s = s[:i0] + """c         --- v3 (2026-10-22): 双一次補間（旧: 距離の逆数の重み）---
c         軸の上 x=0 では横向きの電場がちょうど 0 になる
          fxi = ( xp(k) - x(inex(k),jnex(k)) )
     &        / ( x(inex(k)+1,jnex(k)) - x(inex(k),jnex(k)) )
          fyi = ( yp(k) - y(inex(k),jnex(k)) )
     &        / ( y(inex(k),jnex(k)+1) - y(inex(k),jnex(k)) )
          fxi = min( max( fxi, 0. ), 1. )
          fyi = min( max( fyi, 0. ), 1. )
          Ernex(k) = (1.-fxi)*(1.-fyi)*Er(ijk)
     &             + fxi*(1.-fyi)*Er(ijk+(ny-1))
     &             + fxi*fyi*Er(ijk+ny)
     &             + (1.-fxi)*fyi*Er(ijk+1)
          Eznex(k) = (1.-fxi)*(1.-fyi)*Ez(ijk)
     &             + fxi*(1.-fyi)*Ez(ijk+(ny-1))
     &             + fxi*fyi*Ez(ijk+ny)
     &             + (1.-fxi)*fyi*Ez(ijk+1)
""" + s[i1:]
wr('eforce.f', s)

# ---- 4. 電位：鏡面を格子点の列に ----
s = rd('potential.f')
s = rep(s, """  700 continue
""", """  700 continue

c     --- v3 (2026-10-22): 鏡面（x=0 の列 i=1 と x=間隔 の列 i=nx-1）を
c     ちょうど格子点の上に置く。外側の値＝内側の隣の値（鏡）として
c     ラプラス方程式を解く（旧: φ1=φ2 で鏡面が半マスずれていた）---
      do 705 k=2, ny-2
         dzs = 1./( 0.5*(y(1,k+1) - y(1,k-1)) )
         drs = 1./( x(2,k) - x(1,k) )
         ijk = jj(1,k)
         A(ijk,2) = 0.
         A(ijk,8) = 2.*drs*F1(ijk)
         A(ijk,4) = dzs*F3(ijk-1)
         A(ijk,6) = dzs*F3(ijk)
         A(ijk,5) = -2.*drs*F1(ijk) + dzs*(-F3(ijk)-F3(ijk-1))
         B(ijk)   = 0.
         ijk = jj(nx-1,k)
         A(ijk,8) = 0.
         A(ijk,2) = 2.*drs*F1(ijk-(ny-1))
         A(ijk,4) = dzs*F3(ijk-1)
         A(ijk,6) = dzs*F3(ijk)
         A(ijk,5) = -2.*drs*F1(ijk-(ny-1)) + dzs*(-F3(ijk)-F3(ijk-1))
         B(ijk)   = 0.
  705 continue
""", 'potential.f')
s = rep(s, """      call DBICGST(A, N3, N1, N2, M3, M6, 
     &             B, R0, DD, PP, RR,
     &             S, T, V, XX, EPS, ITR, IER, SGM)
      write(*,*) 'POTENTIAL solved: iterations=',ITR,' res=',EPS,
     &           ' ier=',IER
""", """c     --- v3: 今の解から変化がなくなるまで解き直す（円筒版と同じ）---
      do 995 irep=1, 40
      EPS = 1.0E-14
      ITR = 50000
      call DBICGST(A, N3, N1, N2, M3, M6, 
     &             B, R0, DD, PP, RR,
     &             S, T, V, XX, EPS, ITR, IER, SGM)
      write(*,*) 'POTENTIAL solved: pass=',irep,' iterations=',ITR,
     &           ' res=',EPS,' ier=',IER
      if( ITR .le. 2 ) go to 996
 995  continue
 996  continue
""", 'potential.f')
wr('potential.f', s)

s = rd('ej.f')
s = rep(s, """  210 continue
""", """  210 continue

c     --- v3 (2026-10-22): 鏡面（i=1, i=nx-1）の上では横向きの電場は 0 ---
      do 220 k=1, ny-1
         Er(jj(1,k)) = 0.
         Er(jj(nx-1,k)) = 0.
  220 continue
""", 'ej.f')
wr('ej.f', s)

# ---- 5. 形：一番近い格子点に合わせる ----
s = rd('geo.f')
s = rep(s, """c     --- 2nd nozzle (mirror image at right boundary, x=xlen) ---""",
"""c     --- v3 (2026-10-22): ノズルと電極の端を一番近い格子点に合わせる ---
c     旧: ノズル先端 98.8 mm、電極 91.8〜92.6 mm と 0.2 mm 下にずれていた
c     格子点 i の座標は x = (i-1)*dx、j の座標は y = (j-1)*dy
      inoz  = nint( xnoz/dx ) + 1
      jnoz  = nint( (ylen-ynoz)/dy ) + 1
      ihol  = nint( xhol/dx ) + 1
      jhol1 = nint( yhol1/dy ) + 1
      jhol2 = nint( yhol2/dy ) + 1
      write(*,*) 'v3 nozzle tip y=', y(1,jnoz),
     &           ' half-width x=', x(inoz,1)
      write(*,*) 'v3 electrode y=', y(1,jhol2), '-', y(1,jhol1),
     &           ' x from ', x(ihol,1)

c     --- 2nd nozzle (mirror image at right boundary, x=xlen) ---""", 'geo.f')
wr('geo.f', s)

# ---- 6. 気体 ----
for f in ('diffusion1.f', 'diffusion_gas.f'):
    wr(f, open(os.path.join(SRC, f), encoding='utf-8').read())

# Makefile
s = rd('Makefile')
if '-march=native' not in s:
    s = s.replace('gfortran -c -O3 -fopenmp', 'gfortran -c -O3 -march=native -fopenmp')
    wr('Makefile', s)
print('v3 にしました。rm -f *.o go.out; make でビルドしてください')
