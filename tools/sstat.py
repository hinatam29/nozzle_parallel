import numpy as np,re,sys
want=int(sys.argv[2])
cur=None;rows=[];keep=False
for l in open(sys.argv[1]):
    if 'zone' in l:
        if keep: break
        cur=int(re.search(r"n=\s*(\d+)",l).group(1)); keep=(cur==want); rows=[]
    elif keep and len(l.split())==5:
        rows.append([float(re.sub(r"(\d)([+-]\d{3})$",r"\1E\2",t)) for t in l.split()])
a=np.array(rows); x,y,d,c,po=a.T; f=po==0
r=x[f]*1e3
print(f"t={want*0.2e-3:.1f}ms N={len(a)} 付着={(po==1).sum()} 半径 中央値/90% {np.median(r):.2f}/{np.percentile(r,90):.2f} mm 高さ 中央値/下5% {np.median(y[f])*1e3:.1f}/{np.percentile(y[f],5)*1e3:.1f} mm 液滴径 中央値 {np.median(d)*1e9:.1f} nm")
