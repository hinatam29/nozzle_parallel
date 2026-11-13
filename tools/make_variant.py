"""11/13 教授との相談で決めた試験版を作る（元：multi5random の計算フォルダ）。
使い方:
  python3 tools/make_variant.py <元フォルダ> <新フォルダ> [--ugas U] [--hole R_mm] [--sheet Y_mm V_rel]
  --ugas U       : 噴霧と同じ向きの一様なガスの流れ [m/s]（tools/make_gas.py と同じ。0 なら流れなし）
  --hole R_mm    : 電極の穴の半径 [mm]（今は 5 mm）。geo.f・pcount.f・visualize.py の xhol を変える
  --sheet Y V    : y = Y [mm] より下を全部シート電極にする（電位 V [V]、計算の電位。リング = 0 V、ノズル = +7 kV）。
                   シートより下に入った液滴は pout = 6（シートに付着）として計算から外す
"""
import sys, os, re, subprocess, argparse
ap = argparse.ArgumentParser()
ap.add_argument('src'); ap.add_argument('dst')
ap.add_argument('--ugas', type=float, default=0.0)
ap.add_argument('--hole', type=float, default=None)
ap.add_argument('--sheet', type=float, nargs=2, default=None)
a = ap.parse_args()
here = os.path.dirname(os.path.abspath(__file__))
subprocess.check_call([sys.executable, os.path.join(here, 'make_gas.py'), a.src, a.dst, str(a.ugas)])
# flow.dat は最初と 11000 ステップごと（make_gas.py の後に、11/10 の gasflow と同じ変更）
os.chdir(a.dst)
m = open('main.f').read()
if "mod(istp,11000)" not in m:
    L = m.split('\n')
    i = [k for k, l in enumerate(L) if "write(100,*) 'zone t=\"n= ',istp,'\", i=',nx-1" in l and not l.startswith('c')]
    assert len(i) == 1
    L.insert(i[0], "c     --- flow.dat は最初と 11000 ステップごとだけ書く（ディスク節約）---\n      if( istp .eq. 1 .or. mod(istp,11000) .eq. 0 )then")
    m = '\n'.join(L); old = " 150      continue\n\nc        --- spray (particle) output"; assert m.count(old) == 1
    m = m.replace(old, " 150      continue\n      end if\n\nc        --- spray (particle) output")
    open('main.f', 'w').write(m)
if a.hole is not None:
    r = a.hole * 1e-3
    for fn, pat in [('geo.f', r'^(      xhol = )5\.e-3'), ('pcount.f', r'^(       xhol = )5\.e-3')]:
        s = open(fn).read(); s2, n = re.subn(pat, lambda mo: mo.group(1) + '%g' % r, s, flags=re.M)
        assert n == 1, fn; open(fn, 'w').write(s2)
    if os.path.exists('visualize.py'):
        s = open('visualize.py').read(); s2, n = re.subn(r'^XHOL = 5\.0e-3', 'XHOL = %g' % r, s, flags=re.M)
        if n == 1: open('visualize.py', 'w').write(s2)
if a.sheet is not None:
    ys, vs = a.sheet[0] * 1e-3, a.sheet[1]
    s = open('potential.f').read()
    old = "  780 continue\n"
    assert s.count(old) == 1
    s = s.replace(old, old + """
c     --- 2026-11-13: シート電極（y <= ysht より下を全部電極、電位 vsht）---
      ysht = %g
      vsht = %g
      jsht = nint( ysht/(y(1,2)-y(1,1)) ) + 1
      do 790 i=1, nx-1
      do 790 k=1, jsht
         ijk = jj(i,k)
         A(ijk,5) = 1.0
         B(ijk)   = vsht
         A(ijk,2) = 0.
         A(ijk,4) = 0.
         A(ijk,6) = 0.
         A(ijk,8) = 0.
  790 continue
""" % (ys, vs))
    open('potential.f', 'w').write(s)
    s = open('pcount.f').read()
    old = "       do 20 k=1,ip\n"
    assert s.count(old) == 1
    s = s.replace(old, """c      --- 2026-11-13: シート電極（y <= %g m）に届いた液滴は付着（pout = 6）---
       do 15 k=1, ip
       if( pout(k) .eq. 0.  .and.  yp(k) .le. %g )then
         pout(k) = 6.
       end if
 15    continue
""" % (ys, ys) + old)
    open('pcount.f', 'w').write(s)
print('done', a.dst)
