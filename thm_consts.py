"""Explicit Theorem-level brackets for Channel M (rho=1/2) and S, using the O(eps) constants."""
import math, pickle
ln2=math.log(2)
def eta_p(D,l,k1): return math.log(1+2*D/l)+math.log(1+k1*D)
def eta_m(D,l,k1): return -math.log(1-2*D/l)+math.log(1+k1*D)
out={}
# Channel M, separable (Remark): kappa1=L2/c, L_sigma enters only through D
c=2/math.cosh(1.15)**2; L2=2*(4/(3*math.sqrt(3)))*(2+0.3); Lx=0.3; Ls=0.6; l=1.0; k1=L2/c
vM=pickle.load(open("M.pkl","rb"))[2]
D=lambda r:(1+Ls)*r/(c-Lx)
for e in (1e-2,1e-3,1e-4):
    lam=2*e*(1+e)
    up=eta_p(D(2*e),l,k1)/ln2; lo=(eta_m(D(lam),l,k1)+math.log(1+e))/ln2
    out[("M",e)]=dict(base=math.log2(1/(2*e)),vlo=vM['v_low'],vup=vM['v_up'],eta_up=up,eta_lo=lo)
    print("M",e,round(up,4),round(lo,4))
out["M_consts"]=dict(c=c,L2=L2,k1=k1,Dcoef=(1+Ls)/(c-Lx))
pickle.dump(out,open("thm_consts.pkl","wb"))
