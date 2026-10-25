\\ audC_ext.gp -- exact identities used in the auditor's partial extension of C2 to collinear cubics + 2t+1.
\\ Member X = (a, c) (V = 1/a, u = c/a), other point Q = (a2, c2) (any; Q = (0,0) is 2x+1). Main variable x.
f(A,C) = 2*x + 1 + x*(x+1)*(A*x + C);
htr(phi, q) = {my(r = phi % q, q2 = polcoeff(q,2), q1 = polcoeff(q,1)); polcoeff(r,0) + polcoeff(r,1)*(-q1/q2)/2};
nrm(phi, q) = polresultant(q, phi)/polcoeff(q, 2)^poldegree(phi);
g = (f(a,c) - 1)/x; h = (f(a,c) + 1)/(x+1);
Da = a2 - a; Dc = c2 - c;
H = 4*a^2 - 6*a*c - c^2*(a - c); K = 4*a^2 + a*c*(a - c); G = 2*a^2 + 6*a*c + c^2*(a - c);
print("(1) E^v_X(Q) - (H Da + K Dc)/(2a^3) = ", simplify(htr(-f(a2,c2) - 1, h) - (H*Da + K*Dc)/(2*a^3)));
print("(2) E^w_X(Q) - (G Da - K Dc)/(2a^3) = ", simplify(htr(f(a2,c2) - 1, g) - (G*Da - K*Dc)/(2*a^3)));
print("(3) G + H - 6a^2 = ", simplify(G + H - 6*a^2), ";  R + K/(2a^3) = ", simplify((-2/a - (c/a)*(1-c/a)/2) + K/(2*a^3)));
\\ relative chart: lam = 1 - a2/a, kap = c2 - (a2/a) c
V = 1/a; u = c/a; R = -2*V - u*(1-u)/2; N = u + 2*V; Np = 1 - u + 2*V;
lam = 1 - a2/a; kap = c2 - a2*c/a;
Qw = 4*lam^2 + 2*lam*kap*(1-u) + 2*kap^2*V; Qv = 4*lam^2 - 2*lam*kap*u + 2*kap^2*V;
Mw = 2*(-lam*(1+u) + kap*R) + N*Qw; Mv = 2*(-lam*(2-u) - kap*R) + Np*Qv;
print("(4) N_g(f_Q) - 1 - Mw = ", simplify(nrm(f(a2,c2), g) - 1 - Mw));
print("(5) N_h(f_Q) - 1 - Mv = ", simplify(nrm(f(a2,c2), h) - 1 - Mv));
print("(6) htr_g(f_Q - 1) - (-lam(1+u) + kap R) = ", simplify(htr(f(a2,c2) - 1, g) - (-lam*(1+u) + kap*R)));
print("(7) htr_h(-f_Q - 1) - (-lam(2-u) - kap R) = ", simplify(htr(-f(a2,c2) - 1, h) - (-lam*(2-u) - kap*R)));
print("(8) disc of Qw, Qv in (lam:kap): ", simplify(poldisc(subst(subst(Qw, a2, a*(1 - L)), c2, 0)*0 + (4*y^2 + 2*y*(1-u) + 2*V), y)), " ; ", simplify(poldisc(4*y^2 - 2*y*u + 2*V, y)));
\\ (9) H(c+4, c) and dH/da
print("(9) H(c+4,c) = ", factor(subst(H, a, c+4)), ";  dH/da at a = c+4: ", simplify(subst(deriv(H, a), a, c+4)));
\\ (10) for c in D: a - sqrt(8a) >= -2  <=>  (sqrt a - sqrt 2)^2 >= 0 ; check identity a - 2 sqrt(2a) + 2 = (sqrt a - sqrt 2)^2 with a = s^2
print("(10) s^2 - 2*sqrt2*s + 2 - (s - sqrt2)^2 = ", simplify(s^2 - 2*quadgen(8)*s + 2 - (s - quadgen(8))^2), "  [quadgen(8) = sqrt 2]");
\\ (11) vertical slice: R < 0 <=> c^2 - a c - 4a < 0; roots vs D-slice ends for a > 2
print("(11) 2a^3 R + K = 0 identity checked in (3); D-slice ends: (sqrt(8a))^2 - a sqrt(8a) - 4a and (a - sqrt(8a)) end, with a = s^2/8:");
print("     top end   : ", factor(subst(subst(c^2 - a*c - 4*a, c, s), a, s^2/8)));
print("     bottom end: ", factor(subst(subst(c^2 - a*c - 4*a, c, s^2/8 - s), a, s^2/8)));
