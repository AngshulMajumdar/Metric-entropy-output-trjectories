"""Channel M: y_t = 0.6 y_{t-1} + tanh(2 x_t + 0.3 x_{t-1}), |x|<=rho (k_x=k_y=1, decoupled).
Transfer operator on the past input: (T phi)(xi) = int_{-rho}^{rho} 2 sech^2(2x+0.3 xi) phi(x) dx."""
import pickle, math, sys, numpy as np, scipy.sparse as sp
from core import *
def step_mats(rho_p,rho_q,dk):
    rho=rho_p/rho_q; d=2.0**-dk; N=int(round(2*rho/d)); assert -rho+N*d==rho
    e=[-rho+d*k for k in range(N+1)]
    tl=[None]*(0)
    L=np.zeros((N,N)); U=np.zeros((N,N))
    for i in range(N):
        Xi=I([e[i],e[i+1]])
        for j in range(N):
            Z=2*I([e[j],e[j+1]])+0.3*Xi
            S=2*SECH2(Z)
            L[i,j]=float((S*d).a); U[i,j]=float((S*d).b)
    return L,U,N
out={}
for (p,q,dk) in [(1,2,8),(1,4,9),(1,8,10),(1,16,11),(1,32,12)]:
    rho=p/q
    Lm,Um,N=step_mats(p,q,dk)
    lo,_=cw_lower(sp.csr_matrix(Lm),f"M_low_{q}"); up=cw_upper(sp.csr_matrix(Um),f"M_up_{q}")
    out[q]=dict(rho=rho,N=N,v_low=log2_down(lo,5),v_up=log2_up(up,5),
               union=math.log2(2*math.tanh(2.3*rho)),inter=math.log2(2*math.tanh(1.7*rho)),
               local=math.log2(4*rho),c=2/math.cosh(2.3*rho)**2)
    print(q,out[q],flush=True)
pickle.dump(out,open("M.pkl","wb"))
