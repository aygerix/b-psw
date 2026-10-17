# Real cubic families, normalised f_i = 1 + t g_i, f_i(-1) = -1 (k odd).  TRUE remainders (monic division),
# barrier against a_i -> 0 and against repeated members.  Any sign of a_i allowed.
import numpy as np, sys
k=int(sys.argv[1]); seed=int(sys.argv[2]); trials=int(sys.argv[3])
rng=np.random.default_rng(seed)
def fam(x):
    a=x[:k]; p=x[k:]; return [np.array([a[i],p[i],2-a[i]+p[i],1.0]) for i in range(k)]
def raw(x):
    f=fam(x); out=[]
    for j in range(k):
        P=np.array([1.0])
        for i in range(k):
            if i!=j: P=np.polymul(P,f[i])
        P=np.polysub(P,[1.0])
        g=f[j][:3]; h=np.polydiv(np.polyadd(f[j],[1.0]),[1.0,1.0])[0]
        for m in (g,h):
            r=np.polydiv(P,m/m[0])[1]; r=np.concatenate([np.zeros(2-len(r)),r])[-2:]
            out.append(r/max(1.0,np.max(np.abs(P))))
    return np.concatenate(out)
def res(x):
    V=x.reshape(2,k).T
    pen=1+sum(0.05/(np.sum((V[i]-V[j])**2)+1e-300) for i in range(k) for j in range(i+1,k))+sum(1e-4/(x[i]**2+1e-300) for i in range(k))
    return raw(x)*pen
def lm(x,it=400):
    lam=1e-3; r=res(x); c=r@r
    for _ in range(it):
        J=np.empty((len(r),len(x)))
        for m in range(len(x)):
            e=np.zeros(len(x)); e[m]=1e-7; J[:,m]=(res(x+e)-r)/1e-7
        A=J.T@J; g=J.T@r
        while True:
            try: dx=np.linalg.solve(A+lam*np.diag(np.diag(A)+1e-12),-g)
            except np.linalg.LinAlgError: return x,c
            xn=x+dx; rn=res(xn); cn=rn@rn
            if np.isfinite(cn) and cn<c: x,r,c=xn,rn,cn; lam=max(lam/3,1e-14); break
            lam*=4
            if lam>1e12: return x,c
        if c<1e-26: break
    return x,c
best=[]
for tr in range(trials):
    x=rng.normal(size=2*k)*rng.choice([0.3,1,3]); x,c=lm(x); best.append(c)
    if c<1e-16:
        rr=raw(x); print("ZERO c=%.1e raw=%.1e"%(c,np.max(np.abs(rr))), repr(x.tolist()),flush=True)
print("k=%d trials=%d min residual %.3e  #<1e-16: %d"%(k,trials,min(best),sum(b<1e-16 for b in best)))
