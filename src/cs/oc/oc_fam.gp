\\ Line OC: generic per-member analysis for algebraic prime families.
\\ For a prime p and Selfridge D (p not | D): eps = (D/p), Lucas group order NL = p - eps, Fermat order NF = p - 1.
\\ class: j1 = v2(ord_p 2), j2 = v2(z(p)); computed exactly from the odd parts (no factorisation needed).
\\ usable at B: odd(ord_p 2) is B-smooth AND odd(z(p)) is B-smooth, tested by one exponentiation each:
\\   odd(ord) | S_F (the B-smooth part of odd(p-1))  <=>  2^(2^b S_F) == 1;  same for z with U.
\\ rough bits: log2 of the B-rough part of odd(group order) NOT absorbed (upper bound on new group bits when not usable).
read("oc_lib.gp");
smoothpart(N, B) = { my(S = 1, r); forprime(r = 3, B, while(N % r == 0, N /= r; S *= r)); [S, N]; } \\ N odd in
alpha_pow(k, D, p) = Mod(Mod(1,p)*x, x^2 - x + QofD(D))^k;
isUzero(a) = (polcoeff(lift(lift(a)), 1) == 0);
\\ j2 = min r with U_{m 2^r} == 0, m = odd(NL)
j2of(p, D, NL) = { my(m = oddpart(NL), a = alpha_pow(m, D, p), r = 0);
  while(!isUzero(a), a = a^2; r++; if(r > 4000, return(-1))); r; }
j1of(p) = { my(m = oddpart(p - 1), y = Mod(2, p)^m, r = 0); while(y != 1, y = y^2; r++); r; }
\\ odd z exactly when its B-smooth test passes (reduce S); else return 0
oddz_exact(p, D, NL, B) = { my(m = oddpart(NL), v = valuation(NL, 2), sp = smoothpart(m, B), S = sp[1], fa, r);
  if(!isUzero(alpha_pow(S << v, D, p)), return(0));
  fa = factor(S); for(i = 1, #fa~, r = fa[i,1]; for(j = 1, fa[i,2], if(isUzero(alpha_pow((S/r) << v, D, p)), S /= r, break))); S; }
oddord_exact(p, B) = { my(m = oddpart(p - 1), v = valuation(p - 1, 2), sp = smoothpart(m, B), S = sp[1], fa, r);
  if(Mod(2, p)^(S << v) != 1, return(0));
  fa = factor(S); for(i = 1, #fa~, r = fa[i,1]; for(j = 1, fa[i,2], if(Mod(2, p)^((S/r) << v) == 1, S /= r, break))); S; }
\\ unabsorbed rough part (bits) of odd order on each side, at bound B: test whether the rough part R is fully in the index
roughF(p, B) = { my(m = oddpart(p - 1), sp = smoothpart(m, B), R = sp[2]);
  if(R == 1, return(0.)); if(Mod(2,p)^((p-1)/R) == 1, return(0.)); log(R)/log(2); }
roughL(p, D, NL, B) = { my(m = oddpart(NL), sp = smoothpart(m, B), R = sp[2]);
  if(R == 1, return(0.)); if(isUzero(alpha_pow(NL/R, D, p)), return(0.)); log(R)/log(2); }
\\ one output row per (member, D): fam par bits D eps j1 j2 oz200 oo200 roughF1500 roughL1500 oz1500 oo1500
analyze(p, fam, par, Dl) = { my(D, e, NL, j1, j2, s, oz, oo, oz2, oo2);
  j1 = j1of(p);
  for(k = 1, #Dl, D = Dl[k]; e = kronecker(D, p); if(e == 0, next); NL = p - e;
    j2 = j2of(p, D, NL);
    oo = oddord_exact(p, 200); oz = oddz_exact(p, D, NL, 200);
    oo2 = oddord_exact(p, 1500); oz2 = oddz_exact(p, D, NL, 1500);
    printf("%s %d %d %d %d %d %d %d %d %d %d %.1f %.1f\n", fam, par, #binary(p), D, e, j1, j2, oo, oz, oo2, oz2,
          roughF(p, 1500), roughL(p, D, NL, 1500))); }
