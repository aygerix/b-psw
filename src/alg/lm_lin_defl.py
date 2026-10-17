# Real LINEAR Williams families with MIXED-sign slopes (not covered by Lemma 3 / Lemma 6).
# f_1 = t, f_i = a_i t + b_i.  Conditions: for each j, P_j = prod_{i!=j} f_i equals 1 at both roots of f_j^2-1.
import numpy as np, sys
k=int(sys.argv[1]); seed=int(sys.argv[2]); trials=int(sys.argv[3])
rng=np.random.default_rng(seed)
def coeffs(x):
    a=np.concatenate([[1.0],x[:k-1]]); b=np.concatenate([[0.0],x[k-1:]]); return a,b
def res(x):
    a,b=coeffs(x); out=[]
    for j in range(k):
        for e in (1.0,-1.0):
            r=(e-b[j])/a[j]; P=1.0
            for i in range(k):
                if i!=j: P*=a[i]*r+b[i]
            out.append(P-1.0)
    F=np.stack([a,b],1); pen=1.0+sum(0.1/(np.sum((F[i]-F[j])**2)+1e-300) for i in range(k) for j in range(i+1,k))
    return np.array(out)*pen
def lm(x,it=400):
    lam=1e-3; r=res(x); c=r@r
    for _ in range(it):
        J=np.empty((len(r),len(x))); h=1e-8
        for m in range(len(x)):
            e=np.zeros(len(x)); e[m]=h; J[:,m]=(res(x+e)-r)/h
        A=J.T@J; g=J.T@r
        while True:
            try: dx=np.linalg.solve(A+lam*np.diag(np.diag(A)+1e-12),-g)
            except np.linalg.LinAlgError: return x,c
            xn=x+dx; rn=res(xn); cn=rn@rn
            if np.isfinite(cn) and cn<c: x,r,c=xn,rn,cn; lam=max(lam/3,1e-14); break
            lam*=4
            if lam>1e12: return x,c
        if c<1e-28: break
    return x,c
best=[]; sols=[]
for trial in range(trials):
    x=rng.normal(size=2*(k-1))*rng.choice([0.5,1,2,5])
    x,c=lm(x); best.append(c)
    a,b=coeffs(x)
    if c<1e-20 and min(abs(a))>1e-6:
        F=np.stack([a,b],1); d=min(np.linalg.norm(F[i]-F[j]) for i in range(k) for j in range(i+1,k))
        if d>1e-5:
            key=tuple(np.round(sorted(np.abs(a)/np.max(np.abs(a))),5))
            if key not in sols:
                sols.append(key); print("SOLUTION c=%.1e signs=%s a=%s b=%s"%(c,''.join('+' if s>0 else '-' for s in a),np.round(a,6),np.round(b,6)),flush=True)
print("k=%d trials=%d min residual %.3e  #<1e-20: %d  distinct: %d"%(k,trials,min(best),sum(c<1e-20 for c in best),len(sols)))
