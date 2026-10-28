import numpy as np,re,sys
fr=[];rows=[];cur=None
for l in open(sys.argv[1]):
    if 'zone' in l:
        if rows: fr.append((cur,np.array(rows)))
        cur=int(re.search(r"n=\s*(\d+)",l).group(1)); rows=[]
        if cur>int(sys.argv[2]): break
    elif len(l.split())==5 and 'title' not in l and 'var' not in l:
        rows.append([float(re.sub(r"(\d)([+-]\d{3})$",r"\1E\2",t)) for t in l.split()])
else:
    fr.append((cur,np.array(rows)))
for istp,a in fr[1:]:
    if istp%5000 and istp!=1000: continue
    x,y,d,c,po=a.T; f=po==0
    print(istp*0.2e-3,"ms N",len(a),"dep",(po==1).sum(),"r med/90 %.2f/%.2f"%(np.median(x[f])*1e3,np.percentile(x[f],90)*1e3),"z med/5%% %.1f/%.1f"%(np.median(y[f])*1e3,np.percentile(y[f],5)*1e3))
