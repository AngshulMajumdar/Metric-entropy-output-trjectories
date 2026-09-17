"""Convergence/bias study of the particle estimator (Algorithm 1)."""
import numpy as np, pickle, sys
from particles import run, M, E
part=sys.argv[1]; out={}
try: out=pickle.load(open("pstudy.pkl","rb"))
except Exception: pass
cfgs={"N":[(5000,500,100),(20000,500,100),(50000,500,100)],"T":[(20000,300,100),(20000,1000,100),(20000,1000,300)]}[part]
for ch,name in ((M,"M"),(E,"E")):
    for N,T,B in cfgs:
        est=[run(ch['weight'],ch['step'],-0.5,0.5,N=N,T=T,burn=B,seed=100+s) for s in range(5)]
        out[(name,N,T,B)]=(float(np.mean(est)),float(np.std(est,ddof=1)/np.sqrt(5)))
        print(name,N,T,B,out[(name,N,T,B)],flush=True)
pickle.dump(out,open("pstudy.pkl","wb"))
