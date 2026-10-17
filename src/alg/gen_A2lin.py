# k=5 real quadratic, pattern A2.  n is LINEAR in the free member f5 = x t^2 + y t + z, so the conditions
# n == eps mod M for the six moduli not involving f5 are linear in (x,y,z).  Eliminate (x,y,z) first.
import sympy as sp
t,c,lam,mu,a,b,e,x,y,z=sp.symbols('t c lam mu a b e x y z')
f1=t**2+c; f2=1+lam*(f1-1); f3=a*t**2+b*t+e; f4=-1+mu*(f3+1); f5=x*t**2+y*t+z
n=sp.expand(f1*f2*f3*f4*f5)
mods=[(f1-1,1),(f3-1,1),(f4-1,1),(f1+1,-1),(f2+1,-1),(f3+1,-1)]
eqs=[]
for M,eps in mods:
    r=sp.prem(n-eps,sp.expand(M),t); eqs+=[q for q in sp.Poly(r,t).all_coeffs() if q!=0]
w=sp.symbols('w1:6')
nd=[a,lam,mu,lam-1,mu-1]
eqs+=[sp.expand(nd[i]*w[i]-1) for i in range(5)]
allv=[x,y,z]+[c,lam,mu,a,b,e]+list(w)
with open('A2lin.sing','w') as F:
    F.write(f"ring R=0,({','.join(map(str,allv))}),dp;\noption(redSB);\nideal I=\n")
    F.write(",\n".join(str(sp.expand(q)).replace('**','^') for q in eqs)+";\n")
    F.write('ideal J=eliminate(I,x*y*z);\nideal G=std(J);\nprint("dimJ:"); print(dim(G)); print("sizeJ:"); print(size(G));\nif(size(G)<15){print(G);}\nquit;\n')
print(len(eqs),"eqs")
