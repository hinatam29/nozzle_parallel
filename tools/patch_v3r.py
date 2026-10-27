# -*- coding: utf-8 -*-
"""
v3 版（fieldFIX3_multi5）から「出す位置をランダムにした版」を作るパッチ（2026-10-27）
  python3 patch_v3r.py <fieldFIX3_multi5 の計算フォルダのコピー> [--q 0.8] [--d 500]

  1. 液滴の出す位置（initial.f）
     旧 v3: 毎回 x = 0, 0.025, 0.05 mm の決まった位置 → 噴霧が3本の線に分かれた
     v3r : 毎回、計算する3個とも横の位置を 0〜0.05 mm（＝内径 0.1 mm の幅の半分）
           の中でランダムに選ぶ（一様分布、決まった種で再現できる）。
           (k-1) mod 3 = 0 の液滴は1個だけの液滴（重み 1、自分の鏡像なし）、
           それ以外は鏡像と対の液滴（重み 2、自分の鏡像あり）。合計 5 個は同じ
  2. --q : 帯電量を限界帯電量の何倍にするか（既定 1.0）。
           eforce.f（電場からの力）、pcount.f（クーロン力の加速度）、
           intraction.f（クーロン力の源）の3か所すべてに掛ける
  3. --d : 液滴径の中央値 [nm]（既定 10.5）。対数正規分布の広がり（0.3985）は同じ
"""
import os, re, sys, shutil

args = sys.argv[1:]
D = args[0]
q = 1.0
dmed = None
if '--q' in args:
    q = float(args[args.index('--q') + 1])
if '--d' in args:
    dmed = float(args[args.index('--d') + 1])


def rd(f):
    return open(os.path.join(D, f), encoding='utf-8', errors='replace').read()


def wr(f, s):
    p = os.path.join(D, f)
    if not os.path.exists(p + '.v3'):
        shutil.copy(p, p + '.v3')
    open(p, 'w', encoding='utf-8').write(s)


def rep(s, old, new, f):
    n = s.count(old)
    if n != 1:
        sys.exit('%s: 置き換える文字列が %d 回見つかりました: %r' % (f, n, old[:60]))
    return s.replace(old, new)


if 'v3r' in rd('initial.f'):
    sys.exit('すでに v3r になっています')
print(D, ' 帯電量 x', q, ' 液滴径の中央値', dmed if dmed else 10.5, 'nm')

# 1. ランダムな位置
s = rd('initial.f')
s = rep(s, "      real vp(100000)\n", "      real vp(100000)\n      integer*8 isr\n", 'initial.f')
s = rep(s, "      do 10 k=1,ipmx\n", """c     --- v3r (2026-10-27): 出す位置を毎回ランダムにする（乱数の種は固定）---
      isr = 20261027
      do 10 k=1,ipmx
""", 'initial.f')
s = rep(s, """      kk6 = mod(k-1,6)
      if( kk6 .ge. 3 )then
         xp(k) = xlen - 0.025e-3*real(kk6-3)
      else
         xp(k) = 0.025e-3*real(kk6)
      end if""", """c     --- v3r (2026-10-27): 横の位置は 0〜0.05 mm の中でランダム（一様）---
c     (k-1) mod 6 = 0,1,2 : 左のノズル、3,4,5 : 右のノズル
c     (k-1) mod 3 = 0 は1個だけの液滴（重み1）、1,2 は鏡像と対（重み2）
      isr = mod( isr*16807_8, 2147483647_8 )
      xoff = 0.05e-3*real(isr)/2147483647.
      kk6 = mod(k-1,6)
      if( kk6 .ge. 3 )then
         xp(k) = xlen - xoff
      else
         xp(k) = xoff
      end if""", 'initial.f')
wr('initial.f', s)

# 2. 帯電量
if q != 1.0:
    s = rd('eforce.f')
    s = rep(s, """        vq(k) = pi*sqrt( vq(k) )
""", """        vq(k) = pi*sqrt( vq(k) )
c       --- v3r: 帯電量を限界帯電量の %g 倍にする ---
        vq(k) = %g * vq(k)
""" % (q, q), 'eforce.f')
    wr('eforce.f', s)
    s = rd('pcount.f')
    s = rep(s, """       vq(k) = pi*sqrt( vq(k) )
""", """       vq(k) = pi*sqrt( vq(k) )
c      --- v3r: 帯電量を限界帯電量の %g 倍にする ---
       vq(k) = %g * vq(k)
""" % (q, q), 'pcount.f')
    wr('pcount.f', s)
    s = rd('intraction.f')
    s = rep(s, """            qa(na) = pi*sqrt( 8.*eps0*fs*dp(k)*dp(k)*dp(k) )
""", """            qa(na) = pi*sqrt( 8.*eps0*fs*dp(k)*dp(k)*dp(k) )
c           --- v3r: 帯電量を限界帯電量の %g 倍にする ---
            qa(na) = %g * qa(na)
""" % (q, q), 'intraction.f')
    wr('intraction.f', s)

# 3. 液滴径
if dmed:
    s = rd('rancal.f')
    s = rep(s, "      dmed = 10.50e-9\n", "c     v3r: 液滴径の中央値を変更\n      dmed = %ge-9\n" % dmed, 'rancal.f')
    s = rep(s, "      dp(1) = 10.50e-9", "      dp(1) = %ge-9" % dmed, 'rancal.f')
    wr('rancal.f', s)
    p = os.path.join(D, 'visualize.py')
    if os.path.exists(p):
        s = open(p, encoding='utf-8').read()
        lo = round(dmed * 0.4, -1) or 2
        st = max(round(dmed * 0.1, -1), 1)
        s = s.replace("D_LEVELS = np.arange(4.0, 30.0, 2.0)",
                      "D_LEVELS = np.arange(%g, %g, %g)   # v3r: 液滴径の中央値 %g nm 用" % (lo, dmed * 2.6, st, dmed))
        open(p, 'w', encoding='utf-8').write(s)
print('v3r にしました。rm -f *.o go.out; make でビルドしてください')
