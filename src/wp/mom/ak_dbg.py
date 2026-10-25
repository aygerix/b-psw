import sys, mpmath as mp, ast
mp.mp.dps = 50
import numpy as np
m = int(sys.argv[1]); z = eval(open(sys.argv[2]).read().strip(), {"np": np})
z = [mp.mpc(complex(c)) for c in z]; U, W = z[:m], z[m:]; k = m; L = 3*k
def xser(u, v):
    c = [mp.mpc(0)]*L; c[0] = mp.mpc(1); c[1] = -(1+u)
    for j in range(2, L): c[j] = -(1+u)*c[j-1] - (u+2*v)*c[j-2]
    x = [mp.mpc(0)]*L
    for j in range(3, L): x[j] = v*c[j-3]
    return x
G = [mp.mpc(0)]*L; G[0] = mp.mpc(1)
for i in range(m):
    x = xser(U[i], W[i]); x[0] += 1
    G = [sum(G[a]*x[n-a] for a in range(n+1)) for n in range(L)]
print("resid max:", mp.nstr(max(abs(g) for g in G[k:]), 5))
print("q_i(-1) =", [mp.nstr(1 - (1+U[i]) + (U[i]+2*W[i]), 5) for i in range(m)])
print("u+2v =", [mp.nstr(U[i]+2*W[i], 5) for i in range(m)], " v=", [mp.nstr(W[i],5) for i in range(m)], " u=", [mp.nstr(U[i],5) for i in range(m)])
P1 = sum(G[n]*(-1)**n for n in range(k)); print("P(-1) =", mp.nstr(P1, 10), " expected", 2.0**(1-k))
