\\ pc7_adv_check.gp -- independent adversarial re-derivation of the identities used in the Pencil Lemma proof.
\\ Methods deliberately differ from pc7_ext_verify.gp: root sums are traces of P(M), M the companion matrix of the
\\ monic polynomial (no polsym); the (2,2) pencil is parametrised by r1, r2 directly (not e1, e2); the full quartic
\\ f -+ 1 is used (not only the cofactors K-+).
t; u; x; y;
ok = 1;
chk(name, v) = print(name, ": ", v); if(v != 0, ok = 0);
\\ trace of P(M) for M = companion matrix of g / lead(g)
trP(P, g) = {
  my(gm = g / pollead(g, t), n = poldegree(gm, t), M, A, S);
  M = matrix(n, n, i, j, 0);
  for(i = 2, n, M[i, i - 1] = 1);
  for(i = 1, n, M[i, n] = -polcoef(gm, i - 1, t));
  A = matid(n); S = 0 * matid(n);
  for(m = 0, poldegree(P, t), S += polcoef(P, m, t) * A; A = A * M);
  trace(S);
}
\\ sanity of trP on a numeric example: g = (t-2)(t-3)(t+5), P = t^3 + t: sum = 8+2 + 27+3 + (-125-5) = -90
chk("(0) trP sanity (expect 0)", trP(t^3 + t, (t - 2)*(t - 3)*(t + 5)) + 90);
\\ (I) root-sum identity, generic symbolic f and monic P, d = 1..6
{for(d = 1, 6,
  my(f = c*t^d + sum(m = 0, d - 1, eval(Str("bb", m))*t^m), P = t^d + sum(m = 0, d - 1, eval(Str("qq", m))*t^m));
  chk(Str("(I) d = ", d, " trP(f-1) - trP(f+1) - 2d/c"), simplify(trP(P, f - 1) - trP(P, f + 1) - 2*d/c)))};
\\ (I') complex leading coefficient is NOT allowed to be dropped: identity is 2d/c for any c != 0 (complex too) -- symbolic anyway.
\\ (II) explicit (2,2) pencil parametrised by r1, r2
{
  my(e1 = r1 + r2, e2 = r1*r2, Dl = (1 - r1^2)*(1 - r2^2), lam, mu, G, H, f0, fc, P);
  lam = 2*e1/Dl; mu = 2*(1 - e1^2 + e2)/Dl;
  G = t^2 - 1; H = (t - r1)*(t - r2); f0 = 1 + G*(lam*t + mu); P = G*H;
  chk("(II.a) f0(1) - 1", simplify(subst(f0, t, 1) - 1));
  chk("(II.b) f0(-1) - 1", simplify(subst(f0, t, -1) - 1));
  chk("(II.c) f0(r1) + 1", simplify(subst(f0, t, r1) + 1));
  chk("(II.d) f0(r2) + 1", simplify(subst(f0, t, r2) + 1));
  chk("(II.e) deg f0 <= 3 (expect 0 = coefficient of t^4)", polcoef(f0, 4, t));
  fc = f0 + c*P;
  chk("(II.f) full quartic: trP(P, fc-1) - trP(P, fc+1) - 8/c", simplify(trP(P, fc - 1) - trP(P, fc + 1) - 8/c));
  \\ K-_c and its values at r1, r2
  my(Km = (fc - 1) / G, Kp = (fc + 1) / H);
  chk("(II.g) (fc-1) mod G", simplify((fc - 1) % G)); chk("(II.h) (fc+1) mod H", simplify((fc + 1) % H));
  chk("(II.i) K-(r1) - 2/(1-r1^2)", simplify(subst(Km, t, r1) - 2/(1 - r1^2)));
  chk("(II.j) K-(r2) - 2/(1-r2^2)", simplify(subst(Km, t, r2) - 2/(1 - r2^2)));
  chk("(II.k) lead(K-) - c", simplify(pollead(Km, t) - c));
  chk("(II.l) deg K- - 2", poldegree(Km, t) - 2);
  \\ K+ values at +-1 (mirror sign forcing on the other side): K+(+-1) = 2/H(+-1)
  chk("(II.m) K+(1) - 2/H(1)", simplify(subst(Kp, t, 1) - 2/subst(H, t, 1)));
  chk("(II.n) K+(-1) - 2/H(-1)", simplify(subst(Kp, t, -1) - 2/subst(H, t, -1)));
  chk("(II.o) lead(K+) - c", simplify(pollead(Kp, t) - c));
  \\ step 1: f_i - f_j = (c_i - c_j) P exactly
  chk("(II.p) (f0 + ci P) - (f0 + cj P) - (ci - cj) P", simplify((f0 + ci*P) - (f0 + cj*P) - (ci - cj)*P));
  \\ symmetric slice e1 = 0: lam = 0, mu = 2/(1 - r^2)
  my(lam0 = subst(subst(lam, r1, -r), r2, r), mu0 = subst(subst(mu, r1, -r), r2, r));
  chk("(II.q) lam at r1 = -r2", simplify(lam0)); chk("(II.r) mu - 2/(1 - r^2) at r1 = -r2", simplify(mu0 - 2/(1 - r^2)));
}
\\ (III) half-plane lemma identity |1 + a w|^2 = 1 + 2 a Re w + a^2 |w|^2, w = x + I y
chk("(III.a)", simplify((1 + a*x)^2 + (a*y)^2 - (1 + 2*a*x + a^2*(x^2 + y^2))));
chk("(III.b)", simplify((1 - a*x)^2 + (a*y)^2 - (1 - 2*a*x + a^2*(x^2 + y^2))));
\\ (IV) step 5 with aa = 2/(1 - r^2), H = t^2 - r^2
{
  my(aa = 2/(1 - r^2), G = t^2 - 1, H = t^2 - r^2, f0 = 1 + aa*G, fj = f0 + cj*G*H, phij = 1 + aa*u + cj*u*(u + 2/aa));
  chk("(IV.a) f0 + 1 - aa H", simplify(f0 + 1 - aa*H));
  chk("(IV.b) phi_j(G) - f_j", simplify(subst(phij, u, G) - fj));
  chk("(IV.c) phi_j - 1 - u(cj u + aa + 2cj/aa)", simplify(phij - 1 - u*(cj*u + aa + 2*cj/aa)));
  chk("(IV.d) phi_j + 1 - (u + 2/aa)(cj u + aa)", simplify(phij + 1 - (u + 2/aa)*(cj*u + aa)));
  my(us = -(aa^2 + 2*cj)/(aa*cj));
  chk("(IV.e) phi_j(u*) - 1", simplify(subst(phij, u, us) - 1));
  chk("(IV.f) u* - (-aa/cj - 2/aa)", simplify(us - (-aa/cj - 2/aa)));
  chk("(IV.g) phi_j(-aa/cj) + 1", simplify(subst(phij, u, -aa/cj) + 1));
  chk("(IV.h) -2/aa - (r^2 - 1)", simplify(-2/aa - (r^2 - 1)));
  \\ K-_j in t-variable from phi: (f_j - 1)/G = cj G + aa + 2cj/aa, and at t = r: G(r) = r^2 - 1, value 2/(1 - r^2) = aa
  chk("(IV.i) (cj G + aa + 2cj/aa)(t = r) - aa", simplify(subst(cj*G + aa + 2*cj/aa, t, r) - aa));
}
print(if(ok, "ADV: ALL CHECKS PASSED", "ADV: SOME CHECK FAILED")); print("ok = ", ok);
