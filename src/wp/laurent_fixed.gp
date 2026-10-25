\\ Check (k = 5): with V = y^4 + R(y), R = 15/4 y + 45/16 FIXED, the "-" half Laurent identity
\\   sum_{r<k} y^(r-1) e_r(a) - e_k(a) R(y) == 0
\\ has vanishing z^0 and z^1 coefficients (z = y + 1) identically in the member parameters,
\\ i.e. V(-1) = 2^(1-k), V'(-1) = (2-2k)/2^k hold for every family (expansion at the shared "+" point).
\\ Members: a_i = 2 + z d_i + z^2 b_i (psi_i = 1 + y a_i, psi_i(0) = 1, psi_i(-1) = -1).
k = 5;
bv = vector(k, i, eval(Str("b", i))); dv = vector(k, i, eval(Str("d", i)));
a = vector(k, i, 2 + (y+1)*dv[i] + (y+1)^2*bv[i]);
G = prod(i = 1, k, 1 + X*a[i]);
e = vector(k+1, r, polcoef(G, r-1, X));   \\ e[r+1] = e_r
R = 15/4*y + 45/16;
S = sum(r = 1, k-1, y^(r-1)*e[r+1]) - e[k+1]*R;
Sz = subst(S, y, z - 1);
for(m = 0, 3, print("z^", m, " coefficient: ", factor(polcoef(Sz, m, z))));
print("deg_y S = ", poldegree(S, y));
print("top y coefficient: ", factor(polcoef(S, poldegree(S, y), y)));
print("y^0 coefficient: ", factor(polcoef(S, 0, y)));
