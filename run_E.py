"""Channel E: y_t = 0.5 y + 0.2 sin y + (1+0.3 cos y) tanh(x_t + 0.1 x_{t-1}), |x|<=1/2 (coupled, k_x=k_y=1).
Certified volume exponents of the union / intersection reach-interval relaxations."""
import pickle, math
from core import *
out={}
for name,t in (("union",0.55),("inter",0.45)):
    tI=TANH(I(t))
    fm=lambda Y,tI=tI:0.5*Y+0.2*iv.sin(Y)-(1+0.3*iv.cos(Y))*tI
    fp=lambda Y,tI=tI:0.5*Y+0.2*iv.sin(Y)+(1+0.3*iv.cos(Y))*tI
    fpf=lambda y,t=t:0.5*y+0.2*math.sin(y)+(1+0.3*math.cos(y))*math.tanh(t)
    a=invariant_Y(fpf,fp,I(0),10)
    A,d=outer_kernel(fm,fp,I(0),-a,a,10); up=log2_up(Fr(d)*cw_upper(A,f"E_{name}_up"))
    B,d=inner_kernel(fm,fp,-a,a,10); r,_=cw_lower(B,f"E_{name}_low"); lo=log2_down(Fr(d)*r)
    out[name]=(lo,up,a); print(name,lo,up,a,flush=True)
pickle.dump(out,open("E_rel.pkl","wb"))
