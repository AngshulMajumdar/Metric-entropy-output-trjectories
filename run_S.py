"""Channel S: y_t = 0.6 y + 0.5(1+0.5cos y) sin(3x), |x|<=1 (k_x=0): reach-interval system."""
import pickle, math, sys
from core import *
fm=lambda Y:0.6*Y-0.5*(1+0.5*iv.cos(Y)); fp=lambda Y:0.6*Y+0.5*(1+0.5*iv.cos(Y))
fpf=lambda y:0.6*y+0.5*(1+0.5*math.cos(y))
L=Fr(17,20)
out={}
if sys.argv[1]=="v":
    for dk in (7,8,9,10):
        a=invariant_Y(fpf,fp,I(0),dk)
        A,d=outer_kernel(fm,fp,I(0),-a,a,dk); up=log2_up(Fr(d)*cw_upper(A,f"S_vup_{dk}"))
        B,d=inner_kernel(fm,fp,-a,a,dk); r,_=cw_lower(B,f"S_vlow_{dk}"); lo=log2_down(Fr(d)*r)
        out[dk]=(lo,up,a); print(dk,lo,up,a,flush=True)
    pickle.dump(out,open("S_v.pkl","wb"))
else:
    for q in (20,50,100,200,500,1000):
        eps=1/q; rho=Irat(37,20*q)
        a=invariant_Y(fpf,fp,I(0),12)
        A,v,ep=lattice_graph(fm,fp,eps,-a,a)
        r,W=cw_lower(A,f"S_lat_{q}")
        F0=(v>=-0.75)&(v<=0.75); ok=bool((W[F0]>0).any())
        Y=invariant_Y(fpf,fp,rho,10)
        B,d=outer_kernel(fm,fp,rho,-Y,Y,10); R=cw_upper(B,f"S_vol_{q}")
        ub=log2_up(Fr(d)*R*Fr(q,2))
        rho2=Irat(37,10*q)
        Y2=invariant_Y(fpf,fp,rho2,10)
        B2,d=outer_kernel(fm,fp,rho2,-Y2,Y2,10); R2=cw_upper(B2,f"S_ent_{q}")
        eub=log2_up(Fr(d)*R2*Fr(q,2))
        out[q]=dict(eps=eps,lattice=log2_down(r),start_ok=ok,upper=ub,ent_upper=eub,N=len(v),asym=math.log2(q/2))
        print(q,out[q],flush=True)
    pickle.dump(out,open("S_eps.pkl","wb"))
