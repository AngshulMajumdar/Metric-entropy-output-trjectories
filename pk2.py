import numpy as np, pickle
def run(N=20000,T=500,B=100,seed=0):
    r=np.random.default_rng(seed); x1=np.zeros(N); x2=np.zeros(N); acc=0.0
    for t in range(T):
        x=r.uniform(-.5,.5,N); w=2/np.cosh(2*x+0.2*x1+0.1*x2)**2; m=w.mean()
        if t>=B: acc+=np.log2(m)
        x2=x1; x1=x; idx=r.choice(N,N,p=w/w.sum()); x1=x1[idx]; x2=x2[idx]
    return acc/(T-B)
est=[run(seed=s) for s in range(5)]
out=(float(np.mean(est)),float(np.std(est,ddof=1)/np.sqrt(5))); print(out)
pickle.dump(out,open("pk2.pkl","wb"))
