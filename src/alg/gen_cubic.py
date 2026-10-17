# Real CUBIC Williams families, a_j > 0.  Lemma 6: every real special point is shared by all members;
# a cubic f-1 (f+1) has a real root, so all f_i - 1 share one simple real root r and all f_i + 1 one root s,
# with s < r (else f_i - 1 would have 3 real roots -> all proportional -> Lemma 3).  Normalise r = 0, s = -1:
#   f_i = 1 + t*g_i,  g_i = a_i t^2 + p_i t + q_i,  q_i = 2 - a_i + p_i  (so f_i(-1) = -1).
# P_i(-1) = (-1)^(k-1) forces k odd.  Remaining conditions: P_i == 1 mod g_i and mod h_i = (f_i+1)/(t+1).
import sympy as sp, sys
k=int(sys.argv[1]); t=sp.symbols('t')
A=sp.symbols(f'a1:{k+1}'); Pp=sp.symbols(f'p1:{k+1}'); W=sp.symbols(f'w1:{k+1}')
f=[sp.expand(1+t*(A[i]*t**2+Pp[i]*t+2-A[i]+Pp[i])) for i in range(k)]
eqs=[]
for i in range(k):
    g=sp.expand((f[i]-1)/t); h=sp.quo(f[i]+1,t+1,t); assert sp.rem(f[i]+1,t+1,t)==0
    P=sp.expand(sp.prod([f[l] for l in range(k) if l!=i]))-1
    for m in (g,h):
        r=sp.prem(P,m,t); eqs+=[x for x in sp.Poly(r,t).all_coeffs() if x!=0]
eqs+=[A[i]*W[i]-1 for i in range(k)]
allv=list(A)+list(Pp)+list(W)
with open(f'cub_{k}.sing','w') as F:
    F.write(f"ring R=0,({','.join(map(str,allv))}),dp;\noption(redSB);\nideal I=\n")
    F.write(",\n".join(str(sp.expand(e)).replace('**','^') for e in eqs)+";\n")
    F.write('ideal G=std(I);\nprint("dim:"); print(dim(G)); print("size:"); print(size(G));\nquit;\n')
print(k,len(eqs),"eqs",len(allv),"vars")
