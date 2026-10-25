\\ m9v (adversarial verifier): INDEPENDENT exact check of the m9_forced claims (k = 9, 11), not using
\\ m9_lib / m9_moments / cubic_moments.  V is computed as the Euclidean quotient
\\   V = quo( t^(k-1) prod_i v_i f_i [* (2t+1)],  prod_i v_i (f_i - 1) [* 2t] )
\\ (= polynomial part of t^(k-1) prod(1 + x_i)); v f = v(2t+1) + t(t+1)(t+u) is monic, so the division is exact over Q[u,v].
\\ Structure argument: mod s^k the product prod(1 + x_i) only contains products of <= r cubic x's
\\ (r = 2 for k = 9, r = 3 for k = 11, since x_i = O(s^3)), so every boundary quantity is
\\ C + sum_i A(z_i) + sum_{i<j} B(z_i,z_j) [+ sum_{i<j<l} C3(...)] with m-independent pieces; likewise any
\\ polynomial of degree <= r in the moments M(a,b) = sum v^a w^b.  Hence an identity holding for m = 0..r
\\ symbolic members (k fixed) holds for every m, in particular m = k - lin.
mem(i) = [eval(Str("uu", i)), eval(Str("vv", i))];
Vq(k, lin, L, refl) = {my(N = t^(k-1), D = 1);
  for(i = 1, #L, my(z = L[i], u = if(refl, 1 - z[1], z[1]), v = z[2]);
    N *= v*(2*t + 1) + t*(t + 1)*(t + u); D *= v*(2*t) + t*(t + 1)*(t + u));
  if(lin, N *= 2*t + 1; D *= 2*t);
  my(qr = divrem(N, D, t)); qr[1]};
bndL(k, lin, L) = {my(c0 = 2^(1-k), c1 = (1-k)*2^(1-k), V = Vq(k, lin, L, 0), Vr = Vq(k, lin, L, 1));
  [subst(V, t, -1) - c0, subst(deriv(V, t), t, -1) - c1, subst(Vr, t, -1) - c0, subst(deriv(Vr, t), t, -1) - c1]};
bnd(k, lin, m) = bndL(k, lin, vector(m, i, mem(i)));
\\ moment substitution: M(a,b) -> sum_i v_i^a (u_i - 1/2)^b over the m formal members
msubst(P, m) = {my(R = P, vl = variables(P));
  foreach(vl, z, my(st = Str(z)); if(Vec(st)[1] != "M", next);
    my(sp = strsplit(st, "_"), a = eval(strjoin(Vec(sp[1])[2..#Vec(sp[1])])), b = eval(sp[2]));
    R = subst(R, z, sum(i = 1, m, mem(i)[2]^a * (mem(i)[1] - 1/2)^b)));
  R};
Zm0(E) = (E[1] - E[3])/2;
rq() = (random(2001) - 1000)/(1 + random(97));
Mv(a, b) = sum(i = 1, #PT, PT[i][2]^a * (PT[i][1] - 1/2)^b);
Zm1(E, k) = -(E[2] - E[4])/2 - (k-1)*Zm0(E);
chk(name, k, lin, rmax, f) = {my(ok = 1);
  for(m = 0, rmax, my(E = bnd(k, lin, m), d = simplify(f(E, m)));
    if(d != 0, ok = 0; print("   ", name, ": m = ", m, " residual NONZERO: ", d)));
  print(name, " (k = ", k, ", lin = ", lin, ", m = 0..", rmax, "): ", if(ok, "IDENTITY HOLDS", "FAILS"))};
{
  \\ ---- k = 9 claims, transcribed from the claim text ----
  my(S = M1_0);
  \\ 9 cubics: odd part of (k-1)Z0 + Z1 = S M11 - 75/4 M11 - 5 M13 + 21 M21 ; E2 = -((k-1)Z0 + Z1)
  chk("k9 m9 odd V'", 9, 0, 2, (E, m) -> (E[2] - E[4])/2 + msubst(S*M1_1 - 75/4*M1_1 - 5*M1_3 + 21*M2_1, m));
  chk("k9 m9 even V'", 9, 0, 2, (E, m) -> (E[2] + E[4])/2 + msubst(5/2*S^2 - 351/16*S - 12*M1_2 + 27*M2_0 - M1_4 + 6*M2_2 - 4*M3_0 + 255/32, m));
  chk("k9 8+lin odd V'", 9, 1, 2, (E, m) -> (E[2] - E[4])/2 + msubst(S*M1_1 - 111/8*M1_1 - 9/2*M1_3 + 19*M2_1, m));
  chk("k9 8+lin even V'", 9, 1, 2, (E, m) -> (E[2] + E[4])/2 + msubst(9/4*S^2 - 15*S - 39/4*M1_2 + 89/4*M2_0 - M1_4 + 6*M2_2 - 4*M3_0 + 143/32, m));
  \\ Phi: claimed consequence of the odd equations.  Z0- = (E1 - E1r)/2, Z1 = -E2 - (k-1) E1  =>  Z1- = -(E2-E2r)/2 - (k-1) Z0-
  my(Ph9 = 1149/16*M1_1 + 13*M1_3 - 51*M2_1 - M1_5 + 8*M2_3 - 12*M3_1, Ph8 = 759/16*M1_1 + 21/2*M1_3 - 41*M2_1 - M1_5 + 8*M2_3 - 12*M3_1);
  chk("k9 m9 Phi = -35 Z0- - 9/2 Z1-", 9, 0, 2, (E, m) -> msubst(Ph9, m) - (-35*Zm0(E) - 9/2*Zm1(E, 9)));
  chk("k9 8+lin Phi = -31 Z0- - 4 Z1-", 9, 1, 2, (E, m) -> msubst(Ph8, m) - (-31*Zm0(E) - 4*Zm1(E, 9)));
  \\ also: the odd parts Z0-, Z1- themselves have the M1_1*S coefficient ratio claimed in m9_run9 (9/2 : -35, resp 4 : -31)
  chk("k9 m9 Z0- formula", 9, 0, 2, (E, m) -> (E[1]-E[3])/2 - msubst(9/2*M1_1*M1_0 - 201/16*M1_1 - 19/2*M1_3 + 87/2*M2_1 - M1_5 + 8*M2_3 - 12*M3_1, m));
  chk("k9 8+lin Z0- formula", 9, 1, 2, (E, m) -> (E[1]-E[3])/2 - msubst(4*M1_1*M1_0 - 129/16*M1_1 - 15/2*M1_3 + 35*M2_1 - M1_5 + 8*M2_3 - 12*M3_1, m));
  \\ ---- k = 11 analogues, transcribed from m9_forced.out ----
  my(e11 = -1/6*M1_0^3 + 67/4*M1_0^2 + (6*M1_2 + (-25/2*M2_0 - 2815/64))*M1_0 + (3*M1_1^2 + (-699/16*M1_2 + (903/8*M2_0 + (-69/4*M1_4 + (114*M2_2 + (-250/3*M3_0 + (-M1_6 + (10*M2_4 + (-24*M3_2 + (8*M4_0 + 5115/512)))))))))));
  my(o11 = (35/2*M1_1 + (M1_3 - 4*M2_1))*M1_0 + ((M1_2 + (-2*M2_0 - 381/8))*M1_1 + (-32*M1_3 + (303/2*M2_1 + (-6*M1_5 + (50*M2_3 - 78*M3_1))))));
  my(e10 = -1/6*M1_0^3 + 209/16*M1_0^2 + (11/2*M1_2 + (-23/2*M2_0 - 28))*M1_0 + (11/4*M1_1^2 + (-501/16*M1_2 + (1329/16*M2_0 + (-29/2*M1_4 + (387/4*M2_2 + (-214/3*M3_0 + (-M1_6 + (10*M2_4 + (-24*M3_2 + (8*M4_0 + 2811/512)))))))))));
  my(o10 = (59/4*M1_1 + (M1_3 - 4*M2_1))*M1_0 + ((M1_2 + (-2*M2_0 - 1023/32))*M1_1 + (-99/4*M1_3 + (477/4*M2_1 + (-11/2*M1_5 + (46*M2_3 - 72*M3_1))))));
  chk("k11 m11 even V'", 11, 0, 3, (E, m) -> (E[2] + E[4])/2 + msubst(e11, m));
  chk("k11 m11 odd V'", 11, 0, 3, (E, m) -> (E[2] - E[4])/2 + msubst(o11, m));
  chk("k11 10+lin even V'", 11, 1, 3, (E, m) -> (E[2] + E[4])/2 + msubst(e10, m));
  chk("k11 10+lin odd V'", 11, 1, 3, (E, m) -> (E[2] - E[4])/2 + msubst(o10, m));
  \\ negative control: perturb one coefficient; must FAIL
  chk("NEGATIVE CONTROL (75/4 -> 76/4)", 9, 0, 2, (E, m) -> (E[2] - E[4])/2 + msubst(S*M1_1 - 76/4*M1_1 - 5*M1_3 + 21*M2_1, m));
  \\ ---- direct check at the true size m = k - lin, at random exact rational points (no structure argument) ----
  setrand(424242);
  foreach([[9, 0], [9, 1]], c, my(k = c[1], lin = c[2], m = k - lin, bad = 0);
    for(trial = 1, 6,
      my(pt = vector(m, i, [rq(), abs(rq()) + 1/7]), Ev = bndL(k, lin, pt));
      PT = pt; my(Sv = Mv(1, 0));
      my(odd = if(lin, Sv*Mv(1,1) - 111/8*Mv(1,1) - 9/2*Mv(1,3) + 19*Mv(2,1), Sv*Mv(1,1) - 75/4*Mv(1,1) - 5*Mv(1,3) + 21*Mv(2,1)));
      my(evn = if(lin, 9/4*Sv^2 - 15*Sv - 39/4*Mv(1,2) + 89/4*Mv(2,0) - Mv(1,4) + 6*Mv(2,2) - 4*Mv(3,0) + 143/32,
                       5/2*Sv^2 - 351/16*Sv - 12*Mv(1,2) + 27*Mv(2,0) - Mv(1,4) + 6*Mv(2,2) - 4*Mv(3,0) + 255/32));
      my(ph = if(lin, 759/16*Mv(1,1) + 21/2*Mv(1,3) - 41*Mv(2,1), 1149/16*Mv(1,1) + 13*Mv(1,3) - 51*Mv(2,1)) - Mv(1,5) + 8*Mv(2,3) - 12*Mv(3,1));
      my(z0 = (Ev[1] - Ev[3])/2, z1 = -(Ev[2] - Ev[4])/2 - 8*z0, phc = if(lin, -31*z0 - 4*z1, -35*z0 - 9/2*z1));
      if((Ev[2] - Ev[4])/2 + odd != 0 || (Ev[2] + Ev[4])/2 + evn != 0 || ph - phc != 0, bad++));
    print("direct random-point check k = ", k, ", m = ", m, ", lin = ", lin, ": ", if(bad, Str(bad, " FAILURES"), "6/6 exact agreements")));
}
