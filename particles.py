"""Feynman-Kac particle estimate of v = lim T^-1 log2 int prod_t d_x f (not certified)."""
import numpy as np, pickle
def run(weight, step, lo, hi, N=50000, T=500, burn=100, seed=0):
    rng=np.random.default_rng(seed); xi=np.zeros(N); y=np.zeros(N); acc=0.0
    for t in range(T):
        x=rng.uniform(lo,hi,N); w=weight(x,xi,y); m=w.mean()
        if t>=burn: acc+=np.log2(m)
        y=step(x,xi,y); xi=x
        idx=rng.choice(N,N,p=w/w.sum()); xi=xi[idx]; y=y[idx]
    return acc/(T-burn)+np.log2(hi-lo)
M=dict(weight=lambda x,xi,y:2/np.cosh(2*x+0.3*xi)**2, step=lambda x,xi,y:0.6*y+np.tanh(2*x+0.3*xi))
E=dict(weight=lambda x,xi,y:(1+0.3*np.cos(y))/np.cosh(x+0.1*xi)**2,
       step=lambda x,xi,y:0.5*y+0.2*np.sin(y)+(1+0.3*np.cos(y))*np.tanh(x+0.1*xi))
if __name__=="__main__":
    import sys, os
    name=sys.argv[1]; ch={"M":M,"E":E}[name]
    out=pickle.load(open("particles.pkl","rb")) if os.path.exists("particles.pkl") else {}
    if True:
        est=[run(ch['weight'],ch['step'],-0.5,0.5,seed=s) for s in range(10)]
        out[name]=(float(np.mean(est)),float(np.std(est,ddof=1)/np.sqrt(len(est))))
        print(name,out[name],flush=True)
    pickle.dump(out,open("particles.pkl","wb"))
    