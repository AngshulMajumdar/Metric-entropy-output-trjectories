"""Channel E: independent check of the volume identity Leb(S_T) = int J_T (Lemma 2) for small T.
Hit-or-miss: sample y uniformly in a box, invert the recursion in closed form, test x in X^T."""
import numpy as np, pickle
rng=np.random.default_rng(1)
g=lambda y:0.5*y+0.2*np.sin(y); A=lambda y:1+0.3*np.cos(y)
def member(Y):
    n,T=Y.shape; xi=np.zeros(n); yp=np.zeros(n); ok=np.ones(n,bool)
    for t in range(T):
        z=(Y[:,t]-g(yp))/A(yp); ok&=np.abs(z)<1
        x=np.arctanh(np.clip(z,-0.999999,0.999999))-0.1*xi; ok&=np.abs(x)<=0.5
        xi=x; yp=Y[:,t]
    return ok
def intJ(T,n):
    xi=np.zeros(n); y=np.zeros(n); P=np.ones(n)
    for t in range(T):
        x=rng.uniform(-0.5,0.5,n); u=x+0.1*xi; P*=A(y)/np.cosh(u)**2
        y=g(y)+A(y)*np.tanh(u); xi=x
    return P.mean(), P.std()/np.sqrt(n)
out={}
b1=1.3*np.tanh(0.5); a=1.4375
for T in (1,2,3,4):
    lo=np.array([-b1]+[-a]*(T-1)); hi=-lo; vol=np.prod(hi-lo)
    hits=0; tot=0
    for _ in range(20):
        Y=rng.uniform(lo,hi,(1_000_000,T)); hits+=member(Y).sum(); tot+=len(Y)
    p=hits/tot; hm=(vol*p, vol*np.sqrt(p*(1-p)/tot))
    ij=intJ(T,4_000_000)
    out[T]=dict(hitmiss=hm,intJ=ij); print(T,hm,ij,flush=True)
pickle.dump(out,open("volcheck.pkl","wb"))
