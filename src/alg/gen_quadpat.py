# Lemma 8 reduction for REAL quadratic Williams families (positive leading coeffs):
# all special points non-real (Lemma 6 + Lemma 3), and each side needs a proportional pair.
# Pattern A1: f2-1 = lam (f1-1), f3+1 = mu (f1+1);  A2: f2-1 = lam (f1-1), f4+1 = mu (f3+1).
# Other members free quadratics.  Normalisation f1 = t^2 + c.  Output Singular ideal; dim -1 => no family.
import sympy as sp, sys
k=int(sys.argv[1]); pat=sys.argv[2]
t,c,lam,mu=sp.symbols('t c lam mu')
f=[None]*k; free=[]; inv=[]
f[0]=t**2+c
f[1]=1+lam*(f[0]-1)
if pat=='A1':
    f[2]=-1+mu*(f[0]+1); start=3
else:
    a3,b3,e3=sp.symbols('a3 b3 e3'); f[2]=a3*t**2+b3*t+e3; free+= [a3,b3,e3]; inv.append(a3)
    f[3]=-1+mu*(f[2]+1); start=4
for j in range(start,k):
    a,b,e=sp.symbols(f'a{j+1} b{j+1} e{j+1}'); f[j]=a*t**2+b*t+e; free+=[a,b,e]; inv.append(a)
eqs=[]
for j in range(k):
    P=sp.expand(sp.prod([f[i] for i in range(k) if i!=j]))-1
    for s in (1,-1):
        r=sp.prem(P, sp.expand(f[j]-s), t)
        eqs+= [x for x in sp.Poly(r,t).all_coeffs() if x!=0]
vars_=[c,lam,mu]+free
ws=sp.symbols(f'w0:{3+len(inv)+1}')
# nondegeneracy: lam, mu, a_j nonzero; lam != 1, mu != 1 (else repeated members)
nd=[lam,mu,lam-1,mu-1]+inv
eqs+= [sp.expand(nd[i]*ws[i]-1) for i in range(len(nd))]
allv=vars_+list(ws[:len(nd)])
with open(f'qp_{k}_{pat}.sing','w') as F:
    F.write(f"ring R=0,({','.join(map(str,allv))}),dp;\noption(redSB);\nideal I=\n")
    F.write(",\n".join(str(sp.expand(e)).replace('**','^') for e in eqs)+";\n")
    F.write('ideal G=std(I);\nprint("dim:"); print(dim(G)); print("size:"); print(size(G));\nif(size(G)<40){print(G);}\nquit;\n')
print(k,pat,len(eqs),"eqs",len(allv),"vars")
