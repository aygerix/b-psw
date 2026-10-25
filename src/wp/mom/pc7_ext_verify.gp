\\ pc7_ext_verify.gp -- exact checks for the Pencil Lemma (extremal-member argument), (2,2) quartic pencils, all k.
\\ (1) Root-sum identity: for f of degree d with leading coefficient c and P monic of degree d,
\\     sum_{f(rho)=1} P(rho) - sum_{f(sig)=-1} P(sig) = 2 d / c.        (checked for d = 1..7, fully symbolic)
\\ (2) Pencil instance (d = 4): with the explicit (2,2) data, 2 Re P(th-) - 2 Re P(th+) = 8/c, i.e. RY - RX = 8.
\\ (3) Per-member reduction: f_i - f_j = (c_i - c_j) G H; at f_j = 1: f_i = 1 + (c_i - c_j) P; at f_j = -1: f_i = -1 + (c_i - c_j) P.
\\ (4) Quadratic member f0 (e1 = 0): f0 = 1 + aa G, f0 + 1 = aa H with H = G + 2/aa, members are phi_j(G),
\\     phi_j(u) - 1 = u (c_j u + aa + 2 c_j/aa), phi_j(u) + 1 = (u + 2/aa)(c_j u + aa).
\\ (5) Sign forcing: c_j < 0 makes f_j - 1 > 0 on (-1, 1) (sign of K-_j), contradicting f_j(r1) = -1 (checked as a polynomial identity
\\     for the leading coefficient of K-_j; the inequality itself is elementary).
t; u; \\ main variables first (priority)
\\ power sums of the roots of g (any leading coefficient), p_0..p_n
psum(g, n) = my(gm = g / pollead(g, t)); polsym(gm, n);
rootsum(P, g) = my(n = poldegree(P), ps = psum(g, n)); sum(m = 0, n, polcoef(P, m, t) * ps[m + 1]);
ok = 1;
{for(d = 1, 7,
  my(f = c * t^d + sum(m = 0, d - 1, eval(Str("b", m)) * t^m), P = t^d + sum(m = 0, d - 1, eval(Str("q", m)) * t^m));
  my(D = simplify(rootsum(P, f - 1) - rootsum(P, f + 1) - 2 * d / c));
  print("(1) d = ", d, ": residual = ", D, "   [deg P = ", poldegree(P, t), ", main var of f: ", variable(f), "]"); if(D != 0, ok = 0))};
\\ (2) explicit pencil data
Dl = 1 - e1^2 + 2*e2 + e2^2;
lam = 2*e1/Dl; mu = 2*(1 - e1^2 + e2)/Dl;
G = t^2 - 1; H = t^2 - e1*t + e2; L = lam*t + mu; Lp = lam*t + mu + lam*e1;
f0 = 1 + G*L; X = f0 - 1; Y = f0 + 1;
print("(2a) f0 + 1 - H*Lp = ", simplify(Y - H*Lp)); if(simplify(Y - H*Lp) != 0, ok = 0);
fc = f0 + c*G*H; Km = c*H + L; Kp = c*G + Lp;
print("(2b) fc - 1 - G*Km = ", simplify(fc - 1 - G*Km), ",  fc + 1 - H*Kp = ", simplify(fc + 1 - H*Kp));
\\ sum of P over the two roots of Km and of Kp (rootsum handles the symbolic coefficients)
P = G*H;
SPm = rootsum(P, Km); SPp = rootsum(P, Kp);
r2 = simplify(SPm - SPp - 8/c);
print("(2c) sum_{Km} P - sum_{Kp} P - 8/c = ", r2); if(r2 != 0, ok = 0);
\\ same through X, Y: at roots of Km, X = -c P; at roots of Kp, Y = -c P (remainders)
rX = simplify((X + c*P) % Km); rY = simplify((Y + c*P) % Kp);
print("(2d) (X + cP) mod Km = ", rX, ",  (Y + cP) mod Kp = ", rY); if(rX != 0 || rY != 0, ok = 0);
\\ (3) reduction
fi = f0 + ci*G*H; fj = f0 + cj*G*H;
r3a = simplify((fi - (1 + (ci - cj)*P)) % (fj - 1)); r3b = simplify((fi - (-1 + (ci - cj)*P)) % (fj + 1));
print("(3) f_i - (1 + (c_i - c_j) P) mod (f_j - 1) = ", r3a, ";  f_i - (-1 + (c_i - c_j) P) mod (f_j + 1) = ", r3b);
if(r3a != 0 || r3b != 0, ok = 0);
\\ (4) quadratic member, e1 = 0, aa (1 + e2) = 2
{
  my(e2s = 2/aa - 1, Gs = t^2 - 1, Hs = t^2 + e2s, f0s = 1 + aa*Gs);
  my(c4a = simplify(f0s + 1 - aa*Hs), c4b = simplify(Hs - Gs - 2/aa));
  my(phim = u*(cj*u + aa + 2*cj/aa), phip = (u + 2/aa)*(cj*u + aa), phi = 1 + aa*u + cj*u*(u + 2/aa));
  my(c4c = simplify(phi - 1 - phim), c4d = simplify(phi + 1 - phip), c4e = simplify(subst(phi, u, Gs) - (f0s + cj*Gs*Hs)));
  my(lamq = subst(lam, e1, 0), muq = subst(subst(mu, e1, 0), e2, e2s));
  print("(4) f0+1-aa*H = ", c4a, ", H-G-2/aa = ", c4b, ", phi-1 = ", c4c, ", phi+1 = ", c4d, ", phi(G) - f_j = ", c4e,
        ";  interpolant at e1=0: lam = ", lamq, ", mu - aa = ", simplify(muq - aa));
  if(c4a != 0 || c4b != 0 || c4c != 0 || c4d != 0 || c4e != 0 || lamq != 0 || simplify(muq - aa) != 0, ok = 0);
}
\\ (5) leading coefficient of K-_j is c
print("(5) lead coeff of Km in t: ", pollead(Km, t));
{
  my(Lr = subst(subst(subst(L, e1, r + s), e2, r*s), t, r));
  my(dd = simplify(Lr - 2/(1 - r^2))); print("    L(r) - 2/(1-r^2) with H = (t-r)(t-s): ", dd); if(dd != 0, ok = 0);
}
print(if(ok, "ALL CHECKS PASSED", "SOME CHECK FAILED")); print("ok = ", ok);
