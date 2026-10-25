\\ ak_ident.gp -- exact verification of the identities used by the all-k line (real cubic Williams families).
\\ Member (a, c): f = 2t+1 + t(t+1)(a t + c), g = (f-1)/t, h = (f+1)/(t+1), p = t(t+1).
\\ For a real quadratic q with non-real roots z, zbar:  Tr_q(phi) = phi(z) + phi(zbar) = 2 Re phi(z),
\\ computed exactly from phi mod q.  Moment variables V = 1/a, u = c/a, ut = u - 1/2, X = -c/(2a).
tr(phi, q) = {my(r = phi % q); 2*polcoeff(r, 0, t) + polcoeff(r, 1, t)*(-polcoeff(q, 1, t)/polcoeff(q, 2, t))};
nrm(phi, q) = {polresultant(q, phi, t)/polcoeff(q, 2, t)^poldegree(phi, t)};
fm(a, c) = 2*t + 1 + t*(t+1)*(a*t + c);
gm(a, c) = a*t^2 + (a + c)*t + c + 2;
hm(a, c) = a*t^2 + c*t + 2;
p = t*(t+1);
{
  my(g = gm(a, c), h = hm(a, c));
  print("I1  Tr_g(p) - Tr_h(p)                  = ", simplify(tr(p, g) - tr(p, h)));
  print("I2  Tr_g(t p) - Tr_h(t p) - 6/a        = ", simplify(tr(t*p, g) - tr(t*p, h) - 6/a));
  print("I5  Tr_g(t) - (2X - 1), Tr_h(t) - 2X   = ", simplify(tr(t, g) - (-c/a - 1)), ", ", simplify(tr(t, h) - (-c/a)));
  print("I6  Tr_g(p)/2 - R, R = -2V - u(1-u)/2   = ", simplify(tr(p, g)/2 - (-2/a - (c/a)*(1 - c/a)/2)));
  \\ I3/I4: relative formulas, member l = (a2, c2) evaluated at the special points of member i = (a, c)
  my(f2 = fm(a2, c2), dV = 1/a2 - 1/a, du = c2/a2 - c/a);
  print("I3  [f_l - 1 - a_l(2 dV t + du p)] mod g_i   = ", simplify((f2 - 1 - a2*(2*dV*t + du*p)) % g));
  print("I4  [-f_l - 1 + a_l(2 dV (t+1) + du p)] mod h_i = ", simplify((-f2 - 1 + a2*(2*dV*(t+1) + du*p)) % h));
  \\ I8: linear member 2t+1:  |2w+1|^2 = 8V + 2ut,  |2v+1|^2 = 8V - 2ut
  print("I8  N_g(2t+1) - (8V + 2ut), N_h(2t+1) - (8V - 2ut) = ", simplify(nrm(2*t+1, g) - (8/a + 2*(c/a - 1/2))), ", ", simplify(nrm(2*t+1, h) - (8/a - 2*(c/a - 1/2))));
  \\ K: for i = l + (0, gamma) (same a), N_g(f_l) N_h(f_l) = 1 + K gamma^2 + O(gamma^3), K = 8V ut^2 - 8V^2 - (1/4 - ut^2)^2
  my(Pg = nrm(fm(a, c + gam), g)*nrm(-fm(a, c + gam), h), K = 8/a*(c/a - 1/2)^2 - 8/a^2 - (1/4 - (c/a - 1/2)^2)^2);
  print("I7  [gam^1] = ", simplify(polcoeff(Pg, 1, gam)), ";  [gam^2] - K = ", simplify(polcoeff(Pg, 2, gam) - K));
  \\ symmetric members c = a/2: Re w = -3/4, Re v = -1/4, and f_l(w_i) = 1 + 2(1 - a_l/a_i) w_i exactly
  my(gs = gm(a, a/2), hs = hm(a, a/2));
  print("S1  Tr_g(t) + 3/2, Tr_h(t) + 1/2 (sym)  = ", simplify(tr(t, gs) + 3/2), ", ", simplify(tr(t, hs) + 1/2));
  print("S2  [f_l - 1 - 2(1 - a_l/a) t] mod g (sym) = ", simplify((fm(a2, a2/2) - 1 - 2*(1 - a2/a)*t) % gs));
  print("S3  [-f_l - 1 + 2(1 - a_l/a)(t+1)] mod h (sym) = ", simplify((-fm(a2, a2/2) - 1 + 2*(1 - a2/a)*(t+1)) % hs));
  print("S4  N_g(1 + 2 d t) (sym) = ", simplify(nrm(1 + 2*d*t, gs)), "   [should be 1 - 3d + d^2(2 + 8/a)]: diff ", simplify(nrm(1 + 2*d*t, gs) - (1 - 3*d + d^2*(2 + 8/a))));
}
