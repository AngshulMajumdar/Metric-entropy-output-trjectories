"""Additional certified families: (a) Channel M_beta near the input-dominance boundary;
(b) Channel S_a near the contraction boundary; (c) separable system with k_x=2;
(d) closed-form check u = a x + b x xi + c x^2."""
import pickle, math, sys, numpy as np, scipy.sparse as sp
from core import *
def sep1(kern,rho,N,tag):
    d=2*rho/N; e=[-rho+d*k for k in range(N+1)]
    L=np.zeros((N,N)); U=np.zeros((N,N))
    for i in range(N):
        Xi=I([e[i],e[i+1]])
        for j in range(N):
            S=kern(I([e[j],e[j+1]]),Xi)*d
            L[i,j]=float(S.a); U[i,j]=float(S.b)
    lo,_=cw_lower(sp.csr_matrix(L),tag+"_low"); up=cw_upper(sp.csr_matrix(U),tag+"_up")
    return log2_down(lo,5),log2_up(up,5)
part=sys.argv[1]; out={}
if part=="beta":
    for bn in (0,15,30,45,53):
        b=bn/100; bI=Irat(bn,100)
        lo,up=sep1(lambda X,Xi,bI=bI:2*SECH2(2*X+bI*Xi),0.5,int(sys.argv[2]),f"Mb_{bn}")
        out[bn]=dict(beta=b,c=2/math.cosh(1+b/2)**2,v=(lo,up),vu=math.log2(2*math.tanh(1+b/2)),vi=math.log2(2*math.tanh(1-b/2)))
        print(out[bn],flush=True)
elif part=="a":
    for an in (30,45,60,70,74):
        aI=Irat(an,100); af=an/100
        fm=lambda Y,aI=aI:aI*Y-0.5*(1+0.5*iv.cos(Y)); fp=lambda Y,aI=aI:aI*Y+0.5*(1+0.5*iv.cos(Y))
        fpf=lambda y,af=af:af*y+0.5*(1+0.5*math.cos(y))
        Y=invariant_Y(fpf,fp,I(0),9)
        A,d=outer_kernel(fm,fp,I(0),-Y,Y,9); up=log2_up(Fr(d)*cw_upper(A,f"Sa_{an}_up"))
        B,d=inner_kernel(fm,fp,-Y,Y,9); r,_=cw_lower(B,f"Sa_{an}_low"); lo=log2_down(Fr(d)*r)
        out[an]=dict(a=af,Ls=af+0.25,Y=Y,v=(lo,up)); print(out[an],flush=True)
elif part=="kx2":
    N=int(sys.argv[2]); rho=0.5; d=2*rho/N; e=[-rho+d*k for k in range(N+1)]
    rows=[];cols=[];lo_=[];hi_=[]
    X=[I([e[k],e[k+1]]) for k in range(N)]
    for i1 in range(N):
        for i2 in range(N):
            base=Irat(1,5)*X[i1]+Irat(1,10)*X[i2]
            for j in range(N):
                S=2*SECH2(2*X[j]+base)*d
                rows.append(i1*N+i2); cols.append(j*N+i1); lo_.append(float(S.a)); hi_.append(float(S.b))
    L=sp.csr_matrix((lo_,(rows,cols)),shape=(N*N,N*N)); U=sp.csr_matrix((hi_,(rows,cols)),shape=(N*N,N*N))
    lo,_=cw_lower(L,f"K2_{N}_low"); up=cw_upper(U,f"K2_{N}_up")
    out[N]=(log2_down(lo,4),log2_up(up,4)); print(N,out[N],flush=True)
elif part=="quad":
    a,b,c=1,Fr(2,5),Fr(3,10)
    for q in (2,4,8):
        rho=1/q
        lo,up=sep1(lambda X,Xi:1+Irat(2,5)*Xi+Irat(3,5)*X,rho,128,f"Q_{q}")
        R=Fr(1,q); exact=math.log2(float(R)*a+math.sqrt(float(R*R*a*a+8*b*c*R**4/3)))
        out[q]=dict(rho=rho,v=(lo,up),exact=exact,second=math.log2(2*rho)+float(2*b*c/3)*rho**2/math.log(2))
        print(out[q],flush=True)
elif part=="Mfine":
    lo,up=sep1(lambda X,Xi:2*SECH2(2*X+Irat(3,10)*Xi),0.5,512,"Mfine_512")
    out[512]=(lo,up); print(out,flush=True)
fn=f"fam_{part}.pkl"; pickle.dump(out,open(fn,"wb"))
