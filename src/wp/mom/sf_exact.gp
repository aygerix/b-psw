\\ sf_exact.gp -- exact identities used by the sign-free extension of Theorem 11 (step (iv)).
\\ Normal form f = 2t + 1 + t(t+1)(a t + c), a != 0 of EITHER sign.
f = 2*t + 1 + t*(t+1)*(a*t + c);
g = f - 1; if(polcoef(g, 0, t) != 0, error("f(0) != 1")); g = g / t;
h = f + 1; if(subst(h, t, -1) != 0, error("f(-1) != -1")); h = simplify(h / (t+1));
print("g = ", g, "   h = ", h);
print("g - (a t^2 + (a+c) t + c + 2) = ", g - (a*t^2 + (a+c)*t + c + 2));
print("h - (a t^2 + c t + 2) = ", h - (a*t^2 + c*t + 2));
print("disc(g) - ((a-c)^2 - 8a) = ", poldisc(g, t) - ((a-c)^2 - 8*a));
print("disc(h) - (c^2 - 8a) = ", poldisc(h, t) - (c^2 - 8*a));
print("g(-1) = ", subst(g, t, -1), "   h(0) = ", subst(h, t, 0));
\\ a < 0  =>  both discriminants are > 0 (sum of a square and -8a > 0): extra real special points.
\\ reflection f -> -f(-1-t) maps (a, c) -> (a, a - c)
fr = -subst(f, t, -1 - t);
print("reflection residual: ", fr - subst(subst(f, c, a - c), a, a));
\\ affine normalisation: any cubic phi with phi(r) = 1, phi(r2) = -1; psv(tau) = phi(r + (r - r2) tau)
\\ satisfies psv(0) = 1, psv(-1) = -1, lead(psv) = lead(phi) (r - r2)^3.
phi = A*(x - r)*(x^2 + p*x + q) + 1; \\ phi(r) = 1; impose phi(r2) = -1 by solving for A
Aval = -2 / ((r2 - r)*(r2^2 + p*r2 + q)); phi = subst(phi, A, Aval);
psv = subst(phi, x, r + (r - r2)*tau);
print("psv(0) = ", simplify(subst(psv, tau, 0)), "  psv(-1) = ", simplify(subst(psv, tau, -1)));
print("lead(psv)/(lead(phi)(r-r2)^3) = ", simplify(polcoef(psv, 3, tau) / (Aval*(r - r2)^3)));
\\ Lemma 6 sign rule: e_i in {+-1}, prod_{i != j} e_i = 1 for all j  =>  all e_i equal (and e = -1 only for k odd)
{
  for(k = 2, 9, my(cnt = 0, bad = 0);
    forvec(E = vector(k, i, [0, 1]), my(e = vector(k, i, (-1)^E[i]));
      if(prod(i = 1, k, e[i]) && vecmin(vector(k, j, prod(i = 1, k, if(i == j, 1, e[i])))) == 1
         && vecmax(vector(k, j, prod(i = 1, k, if(i == j, 1, e[i])))) == 1,
        cnt++; if(#Set(e) > 1, bad++)));
    print("k = ", k, ": sign vectors with all P_j = 1: ", cnt, ", non-constant among them: ", bad));
}
\\ Lemma 3 sign-free: a linear member is fixed by its "-" point r1 and "+" point r2: f = 1 + 2 (t - r1)/(r1 - r2)
lin = 1 + 2*(t - r1)/(r1 - r2);
print("linear member through (r1, 1), (r2, -1): f(r1) = ", subst(lin, t, r1), ", f(r2) = ", simplify(subst(lin, t, r2)));
