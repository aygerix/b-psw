\\ audC_ident.gp -- independent exact check of identity C1 (audit of line ak), new code.
\\ Main variable x (highest gp priority). Members: f = 2x+1 + x(x+1)(A x + C), any A != 0 (sign-free).
\\ g = (f-1)/x, h = (f+1)/(x+1). p = x(x+1).
f(A,C) = 2*x + 1 + x*(x+1)*(A*x + C);
p = x*(x+1);
\\ symbolic members i = (ai, ci), l = (al, cl)
gi = (f(ai,ci) - 1)/x; hi = (f(ai,ci) + 1)/(x+1);
print("g_i = ", gi, "   h_i = ", hi);
print("check g_i - [ai(x+1)(x+ci/ai) + 2] = ", simplify(gi - (ai*(x+1)*x + ci*(x+1) + 2)));
print("check h_i - [x(ai x + ci) + 2] = ", simplify(hi - (x*(ai*x+ci) + 2)));
rho = al/ai;
\\ (1) f_l - 1 mod g_i = 2(1-rho) x + (cl - rho ci) p   (reduced mod g_i, both sides)
r1 = (f(al,cl) - 1) % gi;  s1 = (2*(1-rho)*x + (cl - rho*ci)*p) % gi;
print("(1) f_l - 1 - [2(1-rho)x + (cl-rho ci)p] mod g_i = ", simplify(r1 - s1));
\\ (2) f_l + 1 mod h_i = 2(1-rho)(x+1) + (cl - rho ci) p
r2 = (f(al,cl) + 1) % hi;  s2 = (2*(1-rho)*(x+1) + (cl - rho*ci)*p) % hi;
print("(2) f_l + 1 - [2(1-rho)(x+1) + (cl-rho ci)p] mod h_i = ", simplify(r2 - s2));
\\ (3) the a_l (2 dV x + du p) form of the line: a_l dV = 1 - rho, a_l du = cl - rho ci
dV = 1/al - 1/ai; du = cl/al - ci/ai;
print("(3) al*dV - (1-rho) = ", simplify(al*dV - (1-rho)), ";  al*du - (cl - rho ci) = ", simplify(al*du - (cl - rho*ci)));
\\ (4) half-traces: for quadratic q = q2 x^2 + q1 x + q0 and phi, (1/2)(phi(z1)+phi(z2)) computed from phi mod q
htr(phi, q) = {my(r = phi % q, q2 = polcoeff(q,2), q1 = polcoeff(q,1)); polcoeff(r,0) + polcoeff(r,1)*(-q1/q2)/2};
u = ci/ai; V = 1/ai; R = -2*V - u*(1-u)/2;
print("(4a) htr_g(x) + (1+u)/2 = ", simplify(htr(x, gi) + (1+u)/2), ";  htr_g(p) - R = ", simplify(htr(p, gi) - R));
print("(4b) htr_h(x+1) - (2-u)/2 = ", simplify(htr(x+1, hi) - (2-u)/2), ";  htr_h(p) - R = ", simplify(htr(p, hi) - R));
\\ (5) the C1 statement in trace form: htr_g(f_l - 1) = al <n^w, dP>, htr_h(-f_l - 1) = -al <n^v, dP>
nw = [-u-1, R]; nv = [2-u, R]; dP = [dV, du];
print("(5a) htr_g(f_l-1) - al<n^w,dP> = ", simplify(htr(f(al,cl)-1, gi) - al*(nw[1]*dP[1] + nw[2]*dP[2])));
print("(5b) htr_h(-f_l-1) + al<n^v,dP> = ", simplify(htr(-f(al,cl)-1, hi) + al*(nv[1]*dP[1] + nv[2]*dP[2])));
print("(5c) n^v - n^w = ", simplify(nv - nw));
\\ (6) (a,c)-plane form: htr_g(f_l-1) = <m^w, dQ>, htr_h(-f_l-1) = <m^v, dQ>, dQ = (al-ai, cl-ci)
mw = [(1+u)*V - u*R, R]; mv = [(2-u)*V + u*R, -R]; dQ = [al-ai, cl-ci];
print("(6a) htr_g(f_l-1) - <m^w,dQ> = ", simplify(htr(f(al,cl)-1, gi) - (mw[1]*dQ[1] + mw[2]*dQ[2])));
print("(6b) htr_h(-f_l-1) - <m^v,dQ> = ", simplify(htr(-f(al,cl)-1, hi) - (mv[1]*dQ[1] + mv[2]*dQ[2])));
print("(6c) m^w + m^v - (3V, 0) = ", simplify(mw + mv - [3*V, 0]));
\\ (7) linear member 2x+1 (the point Q = (0,0)): same formulas with al = cl = 0 in the (a,c) form
lin = 2*x + 1;
print("(7a) htr_g(lin - 1) - <m^w, -Q_i> = ", simplify(htr(lin - 1, gi) - (mw[1]*(-ai) + mw[2]*(-ci))));
print("(7b) htr_h(-lin - 1) - <m^v, -Q_i> = ", simplify(htr(-lin - 1, hi) - (mv[1]*(-ai) + mv[2]*(-ci))));
\\ (8) norms of the linear member: N_q(phi) = phi(z1) phi(z2) = res(q, phi)/lc(q)^deg(phi)
nrm(phi, q) = polresultant(q, phi)/polcoeff(q, 2)^poldegree(phi);
print("(8a) N_g(2x+1) - (8V + 2u - 1) = ", simplify(nrm(lin, gi) - (8*V + 2*u - 1)));
print("(8b) N_h(2x+1) - (8V - 2u + 1) = ", simplify(nrm(lin, hi) - (8*V - 2*u + 1)));
\\ (9) norms: N_g(x) = |w|^2 = u + 2V, N_h(x+1) = |v+1|^2 = 1 - u + 2V
print("(9) N_g(x) - (u + 2V) = ", simplify(nrm(x, gi) - (u + 2*V)), ";  N_h(x+1) - (1 - u + 2V) = ", simplify(nrm(x+1, hi) - (1 - u + 2*V)));
\\ (10) discriminants: disc(g) = (a-c)^2 - 8a, disc(h) = c^2 - 8a
print("(10) disc g - ((ai-ci)^2 - 8ai) = ", simplify(poldisc(gi) - ((ai-ci)^2 - 8*ai)), ";  disc h - (ci^2 - 8ai) = ", simplify(poldisc(hi) - (ci^2 - 8*ai)));
\\ (11) ray families c = u0 a: f_l - 1 mod g_i = 2(1 - rho) x;  f_l + 1 mod h_i = 2(1-rho)(x+1)
gr = subst(gi, ci, u0*ai); hr = subst(hi, ci, u0*ai);
print("(11) ray: ", simplify((f(al, u0*al) - 1 - 2*(1-al/ai)*x) % gr), ", ", simplify((f(al, u0*al) + 1 - 2*(1-al/ai)*(x+1)) % hr));
\\ (12) negative control: perturb (5a) by a wrong sign of R
print("(12) control (must be nonzero): ", simplify(htr(f(al,cl)-1, gi) - al*(nw[1]*dP[1] - nw[2]*dP[2])) != 0);
