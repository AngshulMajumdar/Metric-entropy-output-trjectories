"""Re-derive every certified number. With VERIFY=1 the integer certificates are loaded from certs/
(no eigensolver); Channel E 2-D enclosures are recomputed deterministically.
Usage: python3 verify.py            (all jobs, ~20 min)
       python3 verify.py 0 3 5      (selected jobs)"""
import os, sys, pickle, subprocess
J=[(["run_S.py","v"],"S_v.pkl","ref/S_v.pkl",None),
   (["run_S.py","eps"],"S_eps.pkl","ref/S_eps.pkl",None),
   (["run_M.py"],"M.pkl","ref/M.pkl",None),
   (["run_E.py"],"E_rel.pkl","ref/E_rel.pkl",None),
   (["run_fam.py","quad"],"fam_quad.pkl","ref/fam_quad.pkl",None),
   (["run_fam.py","a"],"fam_a.pkl","ref/fam_a.pkl",None),
   (["run_fam.py","kx2","32"],"fam_kx2.pkl","ref/fam_kx2_32.pkl",None),
   (["run_fam.py","kx2","48"],"fam_kx2.pkl","ref/fam_kx2_48.pkl",None),
   (["run_fam.py","beta","256"],"fam_beta.pkl","ref/fam_beta.pkl",None),
   (["run_fam.py","Mfine"],"fam_Mfine.pkl","ref/fam_Mfine.pkl",None),
   (["run_E2.py","32","128"],"E2.pkl","ref/E2.pkl",(32,128)),
   (["run_E2.py","64","256"],"E2.pkl","ref/E2.pkl",(64,256)),
   (["run_E2.py","96","384"],"E2.pkl","ref/E2.pkl",(96,384)),
   (["run_E2.py","128","512"],"E2.pkl","ref/E2.pkl",(128,512))]
sel=[int(a) for a in sys.argv[1:]] or range(len(J))
env=dict(os.environ,VERIFY="1"); bad=0
for k in sel:
    cmd,out,ref,key=J[k]
    subprocess.run(["python3"]+cmd,env=env,check=True,stdout=subprocess.DEVNULL)
    a=pickle.load(open(out,"rb")); b=pickle.load(open(ref,"rb"))
    ok=(a[key]==b[key]) if key else (a==b); bad+=not ok
    print(("PASS " if ok else "FAIL ")+f"[{k}] "+" ".join(cmd),flush=True)
print("ALL SELECTED PASS" if bad==0 else f"{bad} FAILURES")
