\\ Line OC library: per-prime orders, 2-adic classes, Selfridge D list.
\\ SD = Selfridge D values that can be selected (prime discriminants or -15; Lemma 12(c)), in Selfridge order.
SD = [5,-7,-11,13,-15,17,-19,-23,29,-31,37,41,-43,-47,53];
QofD(D) = (1-D)/4;
\\ x-coefficient of alpha^k mod p equals U_k (alpha^k = U_k alpha - Q U_{k-1})
Uk(k, D, p) = my(Q = QofD(D), a = Mod(Mod(1,p)*x, x^2 - x + Q)^k); polcoeff(lift(lift(a)),1);
\\ rank of apparition z(p) for D, given the factorisation of N = p - (D/p)
zrank_f(p, D, N, fa) = { my(z = N, r, e);
  for(i = 1, #fa~, r = fa[i,1]; e = fa[i,2];
    for(j = 1, e, if(Uk(z/r, D, p) == 0, z /= r, break)));
  z; }
zrank(p, D) = { my(eps = kronecker(D, p), N = p - eps); zrank_f(p, D, N, factor(N)); }
\\ multiplicative order of 2 given factorisation of p-1
ord2_f(p, fa) = znorder(Mod(2,p), [p-1, fa]);
oddpart(n) = n >> valuation(n, 2);
\\ for p with p+1 = 2^a * c (c small, known), inert for D: compute z quickly
\\ z(p) | p+1; returns z using the explicit factorisation of p+1
zinert_known(p, D, fa) = zrank_f(p, D, p + 1, fa);
