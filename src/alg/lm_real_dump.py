# Numerical search (over C) for quadratic Williams families: f_i^2-1 | P_i - 1, P_i = prod_{j!=i} f_j.
# f_1 = t^2 + c1 (affine normalisation), f_i = a_i t^2 + b_i t + c_i.  Levenberg-Marquardt, random starts.
import numpy as np, sys
rng=np.random.default_rng(int(sys.argv[2]) if len(sys.argv)>2 else 1)
k=int(sys.argv[1]) if len(sys.argv)>1 else 4
nv=1+3*(k-1)
REAL=True
def polys(z):
    f=[np.array([1,0,z[0]])]
    for i in range(k-1): f.append(np.array([z[1+3*i],z[2+3*i],z[3+3*i]]))
    return f
def res(z):
    f=polys(z); out=[]
    for j in range(k):
        P=np.array([1.0+0j])
        for i in range(k):
            if i!=j: P=np.polymul(P,f[i])
        P=np.polysub(P,[1])
        M=np.polysub(np.polymul(f[j],f[j]),[1])
        M=M/M[0]
        q,r=np.polydiv(P,M)
        r=np.concatenate([np.zeros(4-len(r)),r]) if len(r)<4 else r[-4:]
        out.append(r)
    return np.concatenate(out)
def realres(x):
    if REAL: r=res(x[:nv]+0j); return np.concatenate([r.real,r.imag])
    z=(x+0j) if REAL else x[:nv]+1j*x[nv:]; r=res(z); return np.concatenate([r.real,r.imag])
def lm(x,it=300):
    lam=1e-3
    r=realres(x); c=r@r
    for _ in range(it):
        J=np.empty((len(r),len(x))); h=1e-7
        for m in range(len(x)):
            e=np.zeros(len(x)); e[m]=h; J[:,m]=(realres(x+e)-r)/h
        A=J.T@J; g=J.T@r
        while True:
            try: dx=np.linalg.solve(A+lam*np.diag(np.diag(A)+1e-12),-g)
            except np.linalg.LinAlgError: return x,c
            xn=x+dx; rn=realres(xn); cn=rn@rn
            if cn<c: x,r,c=xn,rn,cn; lam=max(lam/3,1e-12); break
            lam*=4
            if lam>1e12: return x,c
        if c<1e-26: break
    return x,c
best=[]
for trial in range(int(sys.argv[3]) if len(sys.argv)>3 else 200):
    x=rng.normal(size=(nv if REAL else 2*nv))*rng.choice([0.5,1,2,4])
    x,c=lm(x)
    z=(x+0j) if REAL else x[:nv]+1j*x[nv:]
    f=polys(z)
    lead=min(abs(fi[0]) for fi in f)
    # distinctness
    dist=min(np.linalg.norm(f[i]-f[j]) for i in range(k) for j in range(i+1,k))
    if c<1e-18:
        print("ZERO c=%.1e lead=%.2e dist=%.2e"%(c,lead,dist), np.round(z.real,5)); sys.stdout.flush()
    best.append(c)
print("k=%d trials=%d min residual %.3e, #<1e-10: %d"%(k,len(best),min(best),sum(b<1e-10 for b in best)))
