# -*- coding: utf-8 -*-
"""
Single（ノズル1本）を xy（平面2次元）で計算する版を作るパッチ（2026-10-27）
  python3 patch_xysingle.py <並列 fieldFIX2 の計算フォルダのコピー>

円筒 rz の Single（1_single/fieldFIX/single_7kV_d11nm_k1e-4）と同じ条件にする：
  - 領域：横 100 mm × 縦 100 mm（rz の半径 100 mm と同じ）、格子 約 0.2 mm
  - ノズルは x = 0 の1本だけ（x = 0 は鏡面＝ノズルの中心）。右のノズルは置かない
  - 対向電極は x = 5 mm から右端（100 mm）までの板（rz の「穴の半径 5 mm の大きな板」と同じ）
  - 液滴は 7 ステップに 1 個、11 nm、吸収係数 1e-4、すべり補正あり（fieldFIX2 と同じ）
  - 電位の左右対称化（並列用）は外す
  - 右端 x = 100 mm は電場 0（rz の外側の境界と同じ）。液滴の鏡像の右側（200 mm 先）は
    遠いので影響はほぼない
違い：xy では奥行き方向に無限に長い「板のノズル」と「板の電極」を計算していることになる
"""
import os, re, sys, shutil
D = sys.argv[1]


def rd(f):
    return open(os.path.join(D, f), encoding='utf-8', errors='replace').read()


def wr(f, s):
    p = os.path.join(D, f)
    if not os.path.exists(p + '.par'):
        shutil.copy(p, p + '.par')
    open(p, 'w', encoding='utf-8').write(s)


def rep(s, old, new, f):
    n = s.count(old)
    if n != 1:
        sys.exit('%s: 置き換える文字列が %d 回見つかりました: %r' % (f, n, old[:60]))
    return s.replace(old, new)


def setxlen(f, indent):
    s = rd(f)
    s, n = re.subn(r'^ +xlen = 1\.1e-2.*$', indent + 'xlen = 1.e-1      ! xy の Single: 領域の横幅 100 mm', s, flags=re.M)
    if n != 1:
        sys.exit(f + ': xlen の行が見つかりません')
    wr(f, s)

setxlen('geo.f', '      ')
setxlen('initial.f', '      ')
setxlen('intraction.f', '      ')
setxlen('pcount.f', '       ')

s = rd('geo.f')
s = rep(s, "      ihol2 = nx - ihol\n", """      ihol2 = nx - ihol
c     --- xy の Single: 電極は x = 5 mm から右端までの板（右の穴なし）---
      ihol2 = nx - 1
""", 'geo.f')
wr('geo.f', s)

s = rd('potential.f')
s = rep(s, "c     ------ nozzle 2 (mirror image, right boundary) ------\n",
        "c     ------ nozzle 2 (mirror image, right boundary) ------\nc     --- xy の Single: 右のノズルは置かない ---\n      go to 745\n", 'potential.f')
s = rep(s, "c     ------ counter electrode (holes mirrored on both sides) ------\n",
        "  745 continue\nc     ------ counter electrode (holes mirrored on both sides) ------\n", 'potential.f')
s = rep(s, "         im = nx - i\n         if( im .gt. i )then\n",
        "         im = nx - i\nc        xy の Single: 左右対称ではないので対称化しない\n         if( .false. )then\n", 'potential.f')
wr('potential.f', s)

s = rd('initial.f')
i0 = s.index("      if( mod(k,2) .eq. 0 )then")
i1 = s.index("      end if", i0) + len("      end if")
s = s[:i0] + "c     xy の Single: すべての液滴をノズル（x = 0）から出す\n      xp(k) = 0." + s[i1:]
wr('initial.f', s)

s = rd('main.f')
s = rep(s, "         if( iph .eq. ns2 ) ip = ip + 1\n",
        "c        xy の Single: ノズル1本なので 7 ステップに 1 個だけ\nc        if( iph .eq. ns2 ) ip = ip + 1\n", 'main.f')
wr('main.f', s)

s = rd('pcount.f')
s = rep(s, "       if( xp(k) .ge. xhol  .and.  xp(k) .le. xlen-xhol\n",
        "c      xy の Single: 電極は x = 5 mm から右端まで\n       if( xp(k) .ge. xhol  .and.  xp(k) .le. xlen\n", 'pcount.f')
wr('pcount.f', s)
print('xy の Single にしました')
