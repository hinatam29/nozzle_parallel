"""一様な空気の流れ（噴霧の向き = -y、速さ ugas [m/s]）を入れた版を作る（感度を調べるための試験版）。
使い方: python3 tools/make_gas.py <元フォルダ> <新フォルダ> <ugas m/s>
- ソース（*.f, subcom.inc, Makefile, run.sh, *.py）だけをコピーし、0 ms から計算する（restart.dat は持ってこない）
- subcom.inc に parameter( ugas = ... ) を追加
- viscous.f：抗力を「液滴と空気の相対速度」で計算する
- pcount.f：抗力の陰解法の式に空気の速度を入れる（u_new = (u + dt*F - dt*beta*ugas) / (1 + dt*beta)）
- 局所配列の大きさ 100000 を m（subcom.inc の配列上限）にそろえる（上限を増やしたときのはみ出し防止）
空気の流れは領域全体で一様（電極やノズルのまわりの流れは考えない）。"""
import sys,os,shutil,re,glob
src,dst,ug=sys.argv[1],sys.argv[2],float(sys.argv[3])
os.makedirs(dst,exist_ok=True)
for p in glob.glob(os.path.join(src,'*')):
    b=os.path.basename(p)
    if os.path.isdir(p): continue
    if b.endswith(('.f','.inc','.py','.sh')) or b=='Makefile':
        shutil.copy2(p,os.path.join(dst,b))
os.chdir(dst)
def chk(fn,s):
    for i,l in enumerate(s.split('\n')):
        if l and l[0] not in 'cC*!' and len(l)>72: raise SystemExit(f'{fn}:{i+1} 72桁を超える: {l}')
# 局所配列 100000 -> m
for fn in glob.glob('*.f'):
    s=open(fn).read()
    s2=re.sub(r'^(\s+(?:real|integer)\s.*)$', lambda mo: mo.group(1).replace('(100000)','(m)'), s, flags=re.M)
    if s2!=s: open(fn,'w').write(s2)
# subcom.inc
s=open('subcom.inc').read(); L=s.split('\n')
hit=[k for k,l in enumerate(L) if l.startswith('      parameter(') and 'm=' in l]
assert len(hit)==1
L.insert(hit[0]+1,"c     2026-11-10: 一様な空気の流れ（噴霧の向き -y）の速さ [m/s]\n      parameter( ugas = %g )"%ug)
open('subcom.inc','w').write('\n'.join(L))
# viscous.f
s=open('viscous.f').read()
n1=s.count('uyp(k)*uyp(k)'); assert n1==3, n1
s=s.replace('uyp(k)*uyp(k)','(uyp(k)+ugas)**2')
old="     &          * uyp(k)\n"; assert s.count(old)==1
s=s.replace(old,"     &          * (uyp(k)+ugas)\n")
s=s.replace("      subroutine viscous\n","      subroutine viscous\nc     2026-11-10: 抗力は空気との相対速度 (uxp, uyp+ugas) で計算\n",1)
chk('viscous.f',s); open('viscous.f','w').write(s)
# pcount.f
s=open('pcount.f').read()
old="       umag = sqrt( uxp(k)*uxp(k) + uyp(k)*uyp(k) )\n"; assert s.count(old)==1
s=s.replace(old,"c      2026-11-10: 空気との相対速度の大きさ\n       umag = sqrt( uxp(k)*uxp(k) + (uyp(k)+ugas)**2 )\n")
old="     &          + dt*gv\n"; assert s.count(old)==1
s=s.replace(old,old+"     &          - dt*beta*ugas\n")
chk('pcount.f',s); open('pcount.f','w').write(s)
print('done',dst,'ugas=',ug)
