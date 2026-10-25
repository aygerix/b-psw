\\ AUDIT: A-reduction lemma. f_j = 1 + c_j A B_j, L_j = (c_jB_j)^(k-1) - prod_{i != j}(c_jB_j - 2c_iB_i).
\\ Claim: gcd(F+1, f_j+1) = gcd(f_j+1, L_j) (as polynomials, with multiplicity). Test on random families where a
\\ common root is FORCED: pick B's, c's, a root theta of L_j, then choose A with c_j A(theta) B_j(theta) = -2
\\ (and optionally a double root: also derivative), so both gcds are nontrivial.
setrand(12345);
rp(d) = 'x^d + sum(i=0, d-1, (random(11)-5)*'x^i);
test(k, a, bvec, dbl) = {
  my(B = vector(k, j, if(bvec[j]==0, 1, rp(bvec[j]))), c = vector(k, j, random(7)+1), j = random(k)+1);
  my(L = (c[j]*B[j])^(k-1) - prod(i=1,k, if(i==j,1, c[j]*B[j] - 2*c[i]*B[i])));
  if (L == 0 || poldegree(L) < 1, return([-1]));
  my(fl = factor(L)[,1], th = fl[1]);   \\ an irreducible factor theta-polynomial of L
  \\ want A with A == -2/(c_j B_j) mod th (mod th^2 if dbl): solve by CRT in Q[x]/(th^e)
  my(e = if(dbl, 2, 1), m = th^e);
  if (poldegree(gcd(c[j]*B[j], m)) > 0, return([-1]));
  my(inv = lift(Mod(c[j]*B[j], m)^(-1)), A0 = lift(Mod(-2*inv, m)));
  \\ A = A0 + m * random of degree a - deg m (if a > deg m), monic not required
  my(A = A0 + m*if(a > poldegree(m), rp(a - poldegree(m)), 0));
  my(f = vector(k, i, 1 + c[i]*A*B[i]), F = prod(i=1,k,f[i]));
  my(g1 = gcd(F+1, f[j]+1), g2 = gcd(f[j]+1, L));
  [j, poldegree(g1), poldegree(g2), g1/pollead(g1) == g2/pollead(g2), poldegree(A)];
}
{
  my(bad = 0, n = 0, nontriv = 0);
  foreach([[3,[0,1,1]],[3,[1,1,2]],[3,[0,2,1]],[4,[1,0,2,1]],[5,[1,1,0,2,1]],[5,[0,0,1,1,2]]], kb,
    for(t = 1, 40, my(r = test(kb[1], 6, kb[2], t%2));
      if (r[1] == -1, next); n++;
      if (r[2] > 0, nontriv++);
      if (!r[4], bad++; print("MISMATCH ", kb, " ", r))));
  print("A-reduction tests: ", n, " families, ", nontriv, " with nontrivial gcd, mismatches: ", bad);
}
\\ L_j == 0 at k = 3 forces a pencil under Korselt (hand proof: X | Y+Z and yz | (y+z)^2). Check the unique
\\ non-pencil solution form X = 2h y z, Y = h y (y+z), Z = h z (y+z): Korselt at member j fails unless y, z constant.
{
  my(y = 'x + 2, z = 'x - 3, h = 'x^2 + 1, X = 2*h*y*z, Y = h*y*(y+z), Z = h*z*(y+z));
  print("harmonic relation L = X(Y+Z)*2 - 4YZ = ", 2*X*(Y+Z) - 4*Y*Z, ";  (Y+Z) mod X = ", (Y+Z) % X);
}
