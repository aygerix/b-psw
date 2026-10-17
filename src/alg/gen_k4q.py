# Generate Singular input: k=4 quadratic Williams family f_i^2-1 | P_i - 1 in C[t].
# Normalise f1 = t^2 + c1 (affine reparametrisation over C). a_i != 0 via inverse vars u_i.
import sympy as sp
t=sp.Symbol('t')
c1=sp.Symbol('c1'); a=sp.symbols('a2:5'); b=sp.symbols('b2:5'); c=sp.symbols('c2:5'); u=sp.symbols('u2:5')
f=[t**2+c1]+[a[i]*t**2+b[i]*t+c[i] for i in range(3)]
eqs=[]
for j in range(4):
    P=sp.expand(sp.prod([f[i] for i in range(4) if i!=j]))-1
    M=sp.expand(f[j]**2-1)
    lc=1 if j==0 else a[j-1]**2
    iv=1 if j==0 else u[j-1]**2
    Mm=sp.expand(M*iv)   # monic modulo a*u=1
    # reduce P modulo Mm treating leading coeff as 1 (valid since a*u=1 is in the ideal)
    R=P
    for deg in range(6,3,-1):
        co=sp.Poly(R,t).coeff_monomial(t**deg)
        R=sp.expand(R-co*t**(deg-4)*Mm)
    for d in range(4):
        eqs.append(sp.Poly(R,t).coeff_monomial(t**d))
eqs+= [a[i]*u[i]-1 for i in range(3)]
vars_=[c1]+list(a)+list(b)+list(c)+list(u)
s="ring R=0,("+",".join(map(str,vars_))+"),dp;\nideal I=\n"+",\n".join(str(sp.expand(e)).replace('**','^') for e in eqs)+";\n"
s+="option(redSB);\nideal G=std(I);\nprint(\"dim:\"); print(dim(G)); print(\"size:\"); print(size(G));\nif(size(G)<20){print(G);}\nquit;\n"
open('k4q.sing','w').write(s); print(len(eqs),"equations,",len(vars_),"vars")
