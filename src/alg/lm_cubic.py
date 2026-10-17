# LM search on the normalised real cubic system (gen_cubic.py, k odd), deflated against repeated members.
import sympy as sp, numpy as np, re, sys
k=int(sys.argv[1]); seed=int(sys.argv[2]); trials=int(sys.argv[3])
s=open(f'cub_{k}.sing').read()
body=re.search(r"ideal I=\n(.*?);\n",s,re.S).group(1)
A=sp.symbols(f'a1:{k+1}'); P=sp.symbols(f'p1:{k+1}')
eqs=[sp.sympify(e.replace('^','**')) for e in body.split(',\n') if 'w' not in e]
F=sp.lambdify(list(A)+list(P),eqs,'numpy')
rng=np.random.default_rng(seed)
def res(x):
    r=np.array(F(*x),dtype=float); V=x.reshape(2,k).T
    pen=1+sum(0.05/(np.sum((V[i]-V[j])**2)+1e-300) for i in range(k) for j in range(i+1,k))
    return r*pen
def lm(x,it=300):
    lam=1e-3; r=res(x); c=r@r
    for _ in range(it):
        J=np.empty((len(r),len(x)))
        for m in range(len(x)):
            e=np.zeros(len(x)); e[m]=1e-7; J[:,m]=(res(x+e)-r)/1e-7
        A_=J.T@J; g=J.T@r
        while True:
            try: dx=np.linalg.solve(A_+lam*np.diag(np.diag(A_)+1e-12),-g)
            except np.linalg.LinAlgError: return x,c
            xn=x+dx; rn=res(xn); cn=rn@rn
            if np.isfinite(cn) and cn<c: x,r,c=xn,rn,cn; lam=max(lam/3,1e-14); break
            lam*=4
            if lam>1e12: return x,c
        if c<1e-26: break
    return x,c
best=[]
for tr in range(trials):
    x=rng.normal(size=2*k)*rng.choice([0.3,1,3])
    x,c=lm(x); best.append(c)
    if c<1e-18: print("ZERO",c,np.round(x,6),flush=True)
print("k=%d trials=%d min residual %.3e"%(k,trials,min(best)))
