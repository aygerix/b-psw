\\ general shapes: write the 8 boundary expressions in evaluation form (at r1, r2 for V; at +-1 for W) and in normalised
\\ variables, to look for a positive combination generalising the symmetric slice.
read("q22_moments.gp");
{
  my(k = 7, B = bnd(k, e1, e2, lam, mu));
  \\ remainder coefficients a1 t + a0 -> evaluations: value at root r is a1 r + a0. Pairs (B1,B2)=R1 coeffs (t^0,t^1), etc.
  my(ev(b0, b1, r) = b1*r + b0);
  my(R1a = ev(B[1], B[2], r1), R1b = ev(B[1], B[2], r2), R2a = ev(B[3], B[4], r1), R2b = ev(B[3], B[4], r2));
  my(R3p = ev(B[5], B[6], 1), R3m = ev(B[5], B[6], -1), R4p = ev(B[7], B[8], 1), R4m = ev(B[7], B[8], -1));
  my(sb(P) = my(Q = subst(subst(P, e1, r1 + r2), e2, r1*r2), Dl = (1 - r1^2)*(1 - r2^2));
     Q = subst(subst(Q, lam, 2*(r1 + r2)/Dl), mu, 2*(1 - r1^2 - r2^2 - r1*r2)/Dl); Q);
  my(E = [R1a, R1b, R2a, R2b, R3p, R3m, R4p, R4m], nm = ["V(r1)","V(r2)","V'(r1)","V'(r2)","W(1)","W(-1)","W'(1)","W'(-1)"]);
  for(i = 1, 8, my(Q = sb(E[i])); print(nm[i], ": ", factor(denominator(Q)), "  numerator terms: ", #Str(numerator(Q))));
  \\ print one of them fully in a compact form: W'(1) - c0 (H^6)'(1)
  print("W(1) numerator = ", numerator(sb(R3p)));
}
