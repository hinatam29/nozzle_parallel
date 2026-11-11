"""11 ms で終わった fieldFIX3_multi5r 系の計算を、restart.dat の続きから長く計算するフォルダを作る。
使い方: python3 tools/make_long.py <元フォルダ> <新フォルダ> <終了時刻ms>
変更: main.f (ftime, ipmx=300000, flow.dat は 11000 ステップごと), subcom.inc (m=300000),
      restart.f (上限が増えても古い restart.dat を読める)。元のファイルは *.11ms。物理は変えない。"""
import sys,os,shutil,re
src,dst,tms=sys.argv[1],sys.argv[2],float(sys.argv[3])
os.makedirs(dst,exist_ok=True)
for f in os.listdir(src):
    p=os.path.join(src,f)
    if f in('flow.dat',) or os.path.isdir(p): continue
    shutil.copy2(p,os.path.join(dst,f))
# flow.dat は最初のコマだけ
with open(os.path.join(src,'flow.dat')) as fi, open(os.path.join(dst,'flow.dat'),'w') as fo:
    nz=0
    for i,l in enumerate(fi):
        if 'zone t=' in l:
            nz+=1
            if nz>=2: break
        fo.write(l)
os.chdir(dst)
for f in('main.f','subcom.inc','restart.f'): shutil.copy2(f,f+'.11ms')
s=open('subcom.inc').read()
L=s.split('\n'); hit=[k for k,l in enumerate(L) if l.startswith('      parameter(') and 'm=100000' in l]
assert len(hit)==1, 'subcom.inc の形が想定と違う'
L[hit[0]]=L[hit[0]].replace('m=100000','m=300000')
open('subcom.inc','w').write('\n'.join(L))
m=open('main.f').read()
m=re.sub(r"\n      ftime = [0-9.e-]+\n", "\nc     2026-11-08: 11 ms の続きから %g ms まで延長\n      ftime = %g.e-3\n"%(tms,tms),m,count=1)
a="      ipmx = 100000.\n"; assert m.count(a)==1
m=m.replace(a,"c     2026-11-08: 液滴配列の上限を 10万 -> 30万（約 70 ms 分）\n      ipmx = 300000.\n")
if "mod(istp,11000)" not in m:
    L=m.split('\n')
    i=[k for k,l in enumerate(L) if "write(100,*) 'zone t=\"n= ',istp,'\", i=',nx-1" in l and not l.startswith('c')]
    assert len(i)==1; L.insert(i[0],"c     --- 2026-11-08: flow.dat は 11000 ステップごとだけ書く（ディスク節約）---\n      if( istp .eq. 1 .or. mod(istp,11000) .eq. 0 )then")
    m='\n'.join(L)
    old=" 150      continue\n\nc        --- spray (particle) output"; assert m.count(old)==1
    m=m.replace(old," 150      continue\n      end if\n\nc        --- spray (particle) output")
open('main.f','w').write(m)
r=open('restart.f').read()
old="      if( nx0 .ne. nx .or. ny0 .ne. ny .or. ipmx0 .ne. ipmx )then"; assert r.count(old)==1
r=r.replace(old,"c     2026-11-08: 上限 ipmx を増やした場合は、古い restart.dat の\nc     ipmx0 個分を読み、残り（ipmx0+1〜ipmx）は initial の値のままにする\n      if( nx0 .ne. nx .or. ny0 .ne. ny .or. ipmx0 .gt. ipmx )then")
a,b=r.split('subroutine rrestart')
open('restart.f','w').write(a+'subroutine rrestart'+b.replace("k=1,ipmx)","k=1,ipmx0)"))
# 局所配列の大きさ 100000 を m にそろえる（ip が 10 万を超えたときのはみ出し防止、2026-11-10）
import glob
for fn in glob.glob('*.f'):
    t=open(fn).read()
    t2=re.sub(r'^(\s+(?:real|integer)\s.*)$', lambda mo: mo.group(1).replace('(100000)','(m)'), t, flags=re.M)
    if t2!=t: open(fn,'w').write(t2)
print('done',dst)
