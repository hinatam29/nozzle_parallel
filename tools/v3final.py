import numpy as np,re,sys
f=sys.argv[1]; xl=float(sys.argv[2])*1e-3; N=int(sys.argv[3])
import subprocess
L=subprocess.run(['tail','-n',str(N),f],capture_output=True,text=True).stdout.splitlines()
a=np.array([[float(re.sub(r"(\d)([+-]\d{3})$",r"\1E\2",t)) for t in l.split()] for l in L if len(l.split())==5])
x,y,d,c,po=a.T; k=np.arange(len(a)); m6=k%6; Lf=m6<3; w=np.where(m6%3==0,1.,2.)
fly=po==0; dep=po==1; dist=np.where(Lf,x,xl-x)*1e3
def wq(v,ww,q):
    i=np.argsort(v); cw=np.cumsum(ww[i]); return v[i][np.searchsorted(cw,q*cw[-1])]
n=w*np.pi/6*d**3*c
up=fly&(y>=88e-3); dn=fly&(y<88e-3)
print(f"付着 {w[dep].sum()/w.sum()*100:.1f}% | 軸から 中央値/90% {wq(dist[fly],w[fly],.5):.2f}/{wq(dist[fly],w[fly],.9):.2f} mm | 高さ 下5%/中央 {wq(y[fly],w[fly],.05)*1e3:.1f}/{wq(y[fly],w[fly],.5)*1e3:.1f} mm | 吸収 合計 {n.sum():.3e} mol（飛行中 z>=88mm {n[up].sum()/n.sum()*100:.0f}%, z<88mm {n[dn].sum()/n.sum()*100:.0f}%, 付着 {n[dep].sum()/n.sum()*100:.0f}%）| 液滴中の濃度 飛行中 {np.average(c[fly],weights=w[fly]):.1f} 付着 {np.average(c[dep],weights=w[dep]) if dep.any() else 0:.1f} mol/m3 | 液滴径 付着 {np.median(d[dep])*1e9 if dep.any() else 0:.1f} 飛行中 {np.median(d[fly])*1e9:.1f} nm")
