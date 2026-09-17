"""Channel E: certified enclosure of v via monotone step operators on the 2-D state (xi, y).
Enclosures are vectorized double-precision evaluations of monotone pieces, widened outward by 1e-12
(far above the <=1 ulp error of the platform's tanh/sin/cos); sums are checked with a 1e-12 relative
safety factor, which dominates floating-point summation error for n<=512 terms."""
import numpy as np, pickle, math, sys, hashlib
from fractions import Fraction as Fr
W=1e-12
def build(n,m,a):
    ex=np.linspace(-0.5,0.5,n+1); ey=np.linspace(-a,a,m+1)
    # u = x + 0.1 xi over box (i: xi cell, j: x cell) -> shape (n,1,n)
    ulo=(ex[None,:-1]+0.1*ex[:-1,None])-W; uhi=(ex[None,1:]+0.1*ex[1:,None])+W
    tlo=np.tanh(ulo)-W; thi=np.tanh(uhi)+W
    amin=np.where((ulo<=0)&(uhi>=0),0.0,np.minimum(abs(ulo),abs(uhi))); amax=np.maximum(abs(ulo),abs(uhi))
    s2lo=1-np.tanh(amax+W)**2-W; s2hi=1-np.tanh(np.maximum(amin-W,0))**2+W
    cy0=np.cos(ey[:-1]); cy1=np.cos(ey[1:])
    Alo=1+0.3*np.minimum(cy0,cy1)-W
    Ahi=1+0.3*np.where((ey[:-1]<=0)&(ey[1:]>=0),1.0,np.maximum(cy0,cy1))+W
    glo=0.5*ey[:-1]+0.2*np.sin(ey[:-1])-W; ghi=0.5*ey[1:]+0.2*np.sin(ey[1:])+W
    T=lambda v:v[:,None,:]; Y=lambda v:v[None,:,None]
    Jlo=Y(Alo)*T(s2lo)/n*(1-W); Jhi=Y(Ahi)*T(s2hi)/n*(1+W)
    plo=np.minimum(Y(Alo)*T(tlo),Y(Ahi)*T(tlo)); phi=np.maximum(Y(Alo)*T(thi),Y(Ahi)*T(thi))
    flo=Y(glo)+plo-W; fhi=Y(ghi)+phi+W
    k0=np.clip(np.searchsorted(ey[1:],flo,'left'),0,m-1); k1=np.clip(np.searchsorted(ey[:-1],fhi,'right')-1,0,m-1)
    return Jlo,Jhi,k0,k1
class RMQ:
    def __init__(s,k0,k1,n):
        L=k1-k0+1; s.l=np.floor(np.log2(L)).astype(int); s.k0=k0; s.k2=k1-(1<<s.l)+1
        s.jj=np.broadcast_to(np.arange(n)[None,None,:],k0.shape); s.lmax=int(s.l.max())
    def q(s,g,op):
        tabs=[g]
        for l in range(1,s.lmax+1):
            p=tabs[-1]; h=1<<(l-1); t=p.copy(); t[:,:-h]=op(p[:,:-h],p[:,h:]); tabs.append(t)
        out=np.empty(s.k0.shape)
        for l in range(s.lmax+1):
            msk=s.l==l
            if msk.any():
                jj=s.jj[msk]; out[msk]=op(tabs[l][jj,s.k0[msk]],tabs[l][jj,s.k2[msk]])
        return out
def run(n,m,a,iters=400):
    Jlo,Jhi,k0,k1=build(n,m,a); R=RMQ(k0,k1,n)
    U=lambda g:(Jhi*R.q(g,np.maximum)).sum(axis=2)
    L=lambda g:(Jlo*R.q(g,np.minimum)).sum(axis=2)
    w=np.ones((n,m))
    for _ in range(iters): w=U(w); w/=w.max()
    rup=float(np.max(U(w)/w))*(1+1e-9); assert w.min()>0
    z=np.ones((n,m))
    for _ in range(iters): z=L(z); z/=z.max()
    Lz=L(z); pos=z>0; rlo=float(np.min(Lz[pos]/z[pos]))*(1-1e-9)
    h=hashlib.sha256(np.ascontiguousarray(w).tobytes()+np.ascontiguousarray(z).tobytes()).hexdigest()
    return math.floor(math.log2(rlo)*1e4)/1e4, math.ceil(math.log2(rup)*1e4)/1e4, h
if __name__=="__main__":
    out={}
    for n,m in [(int(sys.argv[1]),int(sys.argv[2]))]:
        lo,up,h=run(n,m,1.4375); out[(n,m)]=(lo,up); print(n,m,lo,up,h[:12],flush=True)
    try: old=pickle.load(open("E2.pkl","rb"))
    except Exception: old={}
    old.update(out); pickle.dump(old,open("E2.pkl","wb"))
