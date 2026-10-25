\\ pc7_cert_dbl.gp -- (1) the B-encoding of q22_moments.gp at a DOUBLE root r1 = r2 = r: identity B3 = B2 - r B4 for all
\\ lam, mu, w (so only 3 independent V-conditions there, versus 4 for (V - c0 G^6) mod H^2), checked symbolically in r with
\\ 3 symbolic members; (2) the Cauchy-Schwarz inequality used at the corner.
read("pc7_cert_dec.gp");
r;
{
  my(ws = [w, v, u], B = Bof(ws));
  my(Bd = apply(P -> substvec(P, [e1, e2], [2*r, r^2]), B));
  print("double root, 3 symbolic members, symbolic r, lam, mu: B3 - (B2 - r*B4) = ", Bd[3] - (Bd[2] - r*Bd[4]));
  \\ the D-encoding keeps rank 4: its t^3-coefficient (coefficient of (V - c0 G^6) mod H^2) is not a combination of B1..B4
  my(V = Vof(ws), R = substvec((V - c0*G^6) % H^2, [e1, e2], [2*r, r^2]));
  my(B1 = Bd[1], B2 = Bd[2], B4 = Bd[4]);
  \\ Taylor data at r: R = d0 + d1 (t-r) + d2 (t-r)^2 + d3 (t-r)^3 ; B-data determine d0, d1, d2 only
  my(Rt = subst(R, t, t + r), d = vector(4, k, polcoef(Rt, k - 1, t)));
  print("B1 + r B2 - d0 = ", simplify(B1 + r*B2 - d[1]), ",  B2 - d1 = ", simplify(B2 - d[2]), ",  B4 - 2 d2 = ", simplify(B4 - 2*d[3]));
  print("d3 (the condition lost by the B-encoding) is nonzero: ", d[4] != 0);
}
{
  my(f = 140*(1 + x)^2*(3 + 2*x + 3*x^2) - (17 + 26*x + 17*x^2)^2);
  print("corner C-S polynomial: ", f, "   (all coefficients positive: ", vecmin(Vec(f)) > 0, ")");
}
