\\ E2 doubly-primitive census for a general Selfridge D (P = 1, Q = (1-D)/4).
\\ All primes p with ord_p(2) <= X and z_D(p) <= X (z_D = entry point of U_n(1,Q)).
\\ Usage:  echo 'X=10000; DD=-7; OUT="dso_Dm7_X1e4.txt"; read("dso_census_D.gp")' | gp -q -D parisizemax=3000000000
\\ Same method as dso_census.gp: B = prod_{e<=X} Phi_e(alpha,beta), Phi_e = prod_{k|e} U_k^mu(e/k) (an integer for a
\\ nondegenerate Lucas sequence).  p | B iff z_D(p) <= X (p | Q never divides U_k when P = 1).
tt() = getabstime();
QQ = (1 - DD)/4;
T0 = tt();
A = vecprod(vector(X, d, polcyclo(d, 2)));
UU = vector(X); UU[1] = 1; if(X >= 2, UU[2] = 1); for(k = 3, X, UU[k] = UU[k-1] - QQ*UU[k-2]);
if(#select(u -> u == 0, UU), error("degenerate Lucas sequence for D = ", DD));
Phiu(e) = {my(num = 1, den = 1); fordiv(e, k, my(m = moebius(e/k)); if(m == 1, num *= UU[k], m == -1, den *= UU[k])); num/den};
B = vecprod(vector(X, e, abs(Phiu(e))));
TB = tt(); print("D = ", DD, ": A ", exponent(A)+1, " bits, B ", exponent(B)+1, " bits, ", round((TB-T0)/100.)/10., " s");
G = gcd(A, B);
TG = tt(); print("G = gcd(A,B): ", exponent(G)+1, " bits, ", round((TG-TB)/100.)/10., " s");
A = 0; B = 0;
LX = lcm(vector(X, i, i));
lucmod(k, p) = lift((Mod([1, -QQ; 1, 0], p)^k)[2,1]);       \\ U_k(1,Q) mod p
zdesc(p, m) = {my(f = factor(m)); for(i = 1, #f~, my(r = f[i,1]); for(j = 1, f[i,2], if(lucmod(m/r, p) == 0, m /= r, break))); m};
PR = Map(); R = G;
{for(d = 1, X,
  my(g = gcd(R, polcyclo(d, 2)));
  if(g > 1,
    my(f = factor(g)[,1]);
    for(i = 1, #f, my(p = f[i]);
      if(!isprime(p), error("non-prime factor ", p));
      if(!mapisdefined(PR, p), mapput(PR, p, d));
      while(R % p == 0, R /= p))));
}
if(R != 1, error("G not exhausted, cofactor ", R));
rows = List();
{foreach(Vec(PR), p,
  if(p == 2, next);
  my(o, z, kr = kronecker(DD, p));
  o = znorder(Mod(2, p), mapget(PR, p));
  z = if(kr == 0, p, zdesc(p, gcd(p - kr, LX)));
  if(lucmod(z, p) != 0 || o > X || z > X, error("order check failed at ", p));
  listput(rows, [p, o, z, kr]));
}
rows = vecsort(Vec(rows), 1);
TZ = tt(); print("rows ", #rows, "; max p ", rows[#rows][1], "; total ", round((TZ-T0)/100.)/10., " s");
write(OUT, rows);
