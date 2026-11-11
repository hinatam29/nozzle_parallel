import numpy as np,re,sys
f=sys.argv[1]; xl=float(sys.argv[2])*1e-3
txt=open(f).read()
parts=re.split(r'zone t="n=\s*(\d+)\s*", i=\s*(\d+)\s*,f=point',txt)
# parts: [header, istp, n, body, istp, n, body,...]
seen={}
for i in range(1,len(parts)-2,3):
    istp=int(parts[i]); n=int(parts[i+1]); body=parts[i+2]
    body=re.sub(r'(\d)([+-]\d{3})(?=\s)',r'\1E\2',body)
    a=np.array(body.split(),dtype=float)
    if a.size!=n*5: continue
    seen[istp]=a.reshape(n,5)
prev_dep=0;prev_t=0
print("t[ms] N 飛行中 | 付着%(累計) 付着速度[実際の個/0.2ms] | 電極より上(z>=92)の飛行中: 個数 横90%[mm] | z88-92: 個数 | 下から5%[mm] | 上の割合")
for istp in sorted(seen):
    if istp<1000: continue
    a=seen[istp]; x,y,d,c,po=a.T; k=np.arange(len(a)); w=np.where(k%3==0,1.,2.)
    L=(k%6)<3; dist=np.where(L,x,xl-x)*1e3; fly=po==0; dep=po==1
    up=fly&(y>=92e-3); mid=fly&(y<92e-3)&(y>=88e-3)
    D=w[dep].sum(); t=istp*0.2e-3
    if istp%2000==0 or istp in (1000,55000,110000):
        print(f"{t:5.1f} {w.sum():7.0f} {w[fly].sum():7.0f} | {D/w.sum()*100:5.1f} {(D-prev_dep)/((t-prev_t)/0.2):7.1f} | {w[up].sum():6.0f} {np.percentile(dist[up],90):5.2f} | {w[mid].sum():6.0f} | {np.percentile(y[fly],5)*1e3:5.1f} | {w[up].sum()/w[fly].sum()*100:4.1f}%")
        prev_dep=D; prev_t=t
