from mpmath import iv, mpf
iv.prec=53
def tanh_pt(x):            # enclosure of tanh at a point (x may be interval point)
    e=iv.exp(2*iv.mpf(x)); return 1-2/(e+1)
def TANH(X):               # monotone: [lo(tanh(a)), hi(tanh(b))]
    X=iv.mpf(X); return iv.mpf([tanh_pt(X.a).a, tanh_pt(X.b).b])
def ABS(X):
    X=iv.mpf(X)
    if X.a>=0: return X
    if X.b<=0: return -X
    return iv.mpf([0, max(-X.a, X.b)])
def lo(X): return float(X.a)
def hi(X): return float(X.b)
def SECH2(X):
    """enclosure of sech^2 over interval X (even, decreasing in |x|)."""
    X=iv.mpf(X); a,b=X.a,X.b
    lo_abs = 0 if (a<=0<=b) else min(abs(a),abs(b))
    hi_abs = max(abs(a),abs(b))
    t_hi=tanh_pt(hi_abs); t_lo=tanh_pt(lo_abs)
    return iv.mpf([(1-t_hi*t_hi).a, (1-t_lo*t_lo).b])
