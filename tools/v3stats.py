import numpy as np,re,sys
xl=float(sys.argv[2])*1e-3
fr=[];cur=None;rows=[]
for l in open(sys.argv[1]):
    if 'zone' in l:
        if rows: fr.append((cur,np.array(rows)))
        cur=int(re.search(r"n=\s*(\d+)",l).group(1)); rows=[]
    elif len(l.split())==5 and 'var' not in l and 'title' not in l:
        rows.append([float(re.sub(r"(\d)([+-]\d{3})$",r"\1E\2",t)) for t in l.split()])
fr.append((cur,np.array(rows)))
log={}
if len(sys.argv)>3:
    for l in open(sys.argv[3]):
        p=l.split()
        if len(p)==4 and p[0].isdigit(): log[int(p[0])]=float(re.sub(r"(\d)([+-]\d{3,4})$",r"\1E\2",p[3]).replace('E-0','E-'))
for istp,a in fr[1:]:
    k=np.arange(len(a)); m6=k%6; left=m6<3; off=(m6%3)
    x,y,d,c,po=a.T; n=np.pi/6*d**3*c; fly=po==0
    rL=x[left&fly]*1e3; rR=(xl-x[~left&fly])*1e3
    on=(off==0)&fly; dax=np.where(left,x,xl-x)[on]*1e3
    print(f"t={istp*0.2e-3:.2f}ms N={len(a)} dep={(po==1).sum()} | 軸からの距離 中央値/90% 左 {np.median(rL):.3f}/{np.percentile(rL,90):.3f} 右 {np.median(rR):.3f}/{np.percentile(rR,90):.3f} mm | 高さ5%/中央 {np.percentile(y[fly],5)*1e3:.2f}/{np.median(y[fly])*1e3:.2f} | 軸上で出た液滴の軸からの距離 中央値 {np.median(dax)*1e3:.2f} um | 液滴の吸収 {n.sum():.4e} 気体 {log.get(istp,float('nan')):.4e}")
