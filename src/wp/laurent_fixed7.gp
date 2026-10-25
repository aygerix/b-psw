\\ k = 7 check of the universal constraints: with R = A + B z + z^2 (q0 + q1 y), A = 2^(1-k) - 1,
\\ B = (k-1)(1 - 2^(1-k)), the z^0 and z^1 coefficients of sum_{r<k} y^(r-1) e_r(a) - e_k(a) R vanish identically.
k = 7;
a = vector(k, i, 2 + (y+1)*eval(Str("d", i)) + (y+1)^2*eval(Str("b", i)));
G = prod(i = 1, k, 1 + X*a[i]); e = vector(k+1, r, polcoef(G, r-1, X));
A = 2^(1-k) - 1; B = (k-1)*(1 - 2^(1-k));
R = A + B*(y+1) + (y+1)^2*(q0 + q1*y);
S = sum(r = 1, k-1, y^(r-1)*e[r+1]) - e[k+1]*R;
Sz = subst(S, y, z - 1);
print("z^0: ", polcoef(Sz, 0, z) == 0, "   z^1: ", polcoef(Sz, 1, z) == 0, "   deg_y S = ", poldegree(S, y));
