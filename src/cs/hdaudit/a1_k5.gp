\\ AUDIT (new code) of the k = 5 HD counterexample. Linear common-root family f_j = 1 + s_j x.
\\ Korselt holds automatically (F(0) = 1).  Member j Lucas-good  <=>  F(-2/s_j) = -1  <=>  prod_{i != j}(1 - 2 s_i/s_j) = 1.
\\ Independent identity: for a root y of G with other roots y', y'': y' + y'' = e1 - y, y' y'' = e3/y, hence
\\ y*[(y-2)(y-2M)(y-2y')(y-2y'') - y^4] = (y-2)(y-2M)(y^3 - 2y^2(e1-y) + 4e3) - y^5.  Must be 0 mod G.
default(realprecision, 120);
G(a,M) = (7*M^3+57*M^2+57*M+7)*a^3 - (98*M^3+308*M^2+98*M)*a^2 + (392*M^3+392*M^2)*a - 448*M^3;
{
  my(g = G('y,'M), lc = pollead(g,'y), e1 = -polcoeff(g,2,'y)/lc, e3 = polcoeff(g,0,'y)/lc*(-1));
  \\ for monic cubic y^3 - e1 y^2 + e2 y - e3: e3 = -c0/lc
  my(E = ('y-2)*('y-2*'M)*('y^3 - 2*'y^2*(e1-'y) + 4*e3) - 'y^5);
  print("[sym] remainder of y*E(y) mod G in Q(M)[y]: ", lift(Mod(E, g)) );
  \\ numeric check at several M, incl. near M0
  my(dM = poldisc(g,'y), fs = factor(dM));
  print("[disc] factorisation: ", fs);
  my(sx = fs[#fs[,1],1]);
  print("[disc] real roots of sextic: ", polroots(sx));
  foreach([20, 24.8, 24.86, 24.87, 25, 30, 135, 256, 1000, 10^6, 1/30, 0.03, -2, -119/54], M,
     my(p = G('y, M), n = polsturm(p));
     my(r = [real(z) | z <- polroots(p), abs(imag(z)) < 1e-60]);
     my(s = concat([1], concat(r, [M])), R = vector(#s, j, prod(i=1,#s, if(i==j,1, 1-2*s[i]/s[j]))));
     printf("[M=%s] #real roots=%d  slopes=%s\n     R_j-1=%s\n", M, n, s, apply(t->if(abs(t-1)<1e-80,0,t-1), R)));
  \\ G(1;M), G(M;M): roots must avoid M > M0
  print("[avoid] G(1;M) = ", factor(G(1,'M)), "  G(M;M) = ", factor(G('M,'M)));
  print("[avoid] real roots of G(1;M): ", polroots(G(1,'M)));
  print("[avoid] real roots of G(M;M)/M^3: ", polroots(G('M,'M)/'M^3));
}
\\ The rational mixed-sign orbit point and the integer example
{
  my(chk(s) = vector(#s, j, prod(i=1,#s, if(i==j,1, 1-2*s[i]/s[j]))));
  my(s1 = [1, 238/41, -476/95, 952/481, -119/54]);
  print("[rat] R_j for ", s1, " : ", chk(s1));
  my(s2 = [101168730, 587272140, -506908584, 200234160, -222945905]);
  print("[int] ratios s2/s2[1]: ", s2/s2[1]);
  print("[int] R_j: ", chk(s2));
  my(F = prod(j=1,5, 1 + s2[j]*'x));
  print("[int] Korselt rem: ", vector(5,j, (F-1) % (s2[j]*'x)), "  deg gcd(F+1,f_j+1): ", vector(5,j, poldegree(gcd(F+1, 2 + s2[j]*'x))));
}
\\ G5 over Q at M = 135, 256: irreducible?
print("[irr] M=135: ", polisirreducible(G('y,135)), "  M=256: ", polisirreducible(G('y,256)));
