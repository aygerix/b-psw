\\ (2,2) real quartic Williams pencils: f_j = f0 + c_j G H, G = t^2 - 1 (shared "-" points +-1), H = t^2 - e1 t + e2
\\ (shared "+" points r1, r2 in (-1, 1)), f0 = 1 + G (lam t + mu), lam = 2 e1/Dl, mu = 2 (1 - e1^2 + e2)/Dl,
\\ Dl = (1 - r1^2)(1 - r2^2) = 1 - e1^2 + 2 e2 + e2^2. f_j - 1 = G K-_j, K-_j = c_j H + L, L = lam t + mu;
\\ f_j + 1 = H K+_j, K+_j = c_j G + L', L' = lam t + mu + lam e1.
\\ Cofactors: F - 1 = G prod K-_j V, V = G^(k-1) [prod(1 + x_j) - prod x_j], x_j = 1/(f_j - 1);
\\           F + 1 = H prod K+_j W, W = H^(k-1) [prod(1 - y_j) + prod y_j], y_j = 1/(f_j + 1).
\\ Boundary (k odd): at the roots of H, x_j = -1/2, so V == 2^(1-k) G^(k-1) and V' == 2^(1-k) (G^(k-1))' mod H;
\\                   at +-1, y_j = 1/2, so W == 2^(1-k) H^(k-1) and W' == 2^(1-k) (H^(k-1))' mod G.
\\ Members enter only through w_j = 1/c_j; the series are written in POWER SUMS p_a = sum_j w_j^a (formal variables).
s; t; w; trunc(P, n) = sum(j = 0, n - 1, polcoef(P, j, s) * s^j);
\\ single-member series x(w) up to s^(N-1): x = w s^4 / ((1 - s^2) (1 + (-e1) s + e2 s^2 + w s (lam + mu s)))
inv1(D, N) = {my(r = 1, e = 1 - D); for(i = 1, N, r = trunc(1 + e * r, N)); r};   \\ 1/D for D = 1 + O(s)
mkx(N, w, lam, mu, e1, e2) = trunc(w * s^4 * inv1(trunc((1 - s^2) * (1 - e1*s + e2*s^2 + w*s*(lam + mu*s)), N), N), N);
mky(N, w, lam, mu, e1, e2) = trunc(w * s^4 * inv1(trunc((1 - e1*s + e2*s^2) * (1 - s^2 + w*(lam*s + (mu + lam*e1)*s^2)), N), N), N);
\\ log(1 + X) truncated, X = O(s^4)
lg1(X, N) = {my(r = 0, P = 1); for(m = 1, (N - 1) \ 4, P = trunc(P * X, N); r += (-1)^(m+1) * P / m); r};
\\ replace w^a by p_a (power sums); P polynomial in w with no constant term
topow(P, w, pv) = sum(a = 1, poldegree(P, w), polcoef(P, a, w) * pv[a]);
ex1(Lg, N) = {my(r = 1, P = 1); for(m = 1, (N - 1) \ 4, P = trunc(P * Lg, N); r += P / m!); r};
\\ build the 8 boundary equations for k members (all quartic pencil members)
\\ bnd(k, e1, e2, lam, mu): numbers or formal symbols (pass lam, mu as symbols to avoid rational functions; the relations
\\ lam Dl = 2 e1, mu Dl = 2 (1 - e1^2 + e2) are then imposed separately)
\\ quad = 1: one extra member f0 itself (quadratic; needs e1 = 0, lam = 0, f0 - 1 = aa G, f0 + 1 = aa H, aa (1 + e2) = 2);
\\ then k = (number of quartic members) + 1, and the quartic members are k - 1 of them.
bnd(k, e1, e2, lam, mu, quad = 0) = {
  my(N = 2*k - 1, pv = vector(2*k, a, eval(Str("p", a))));
  my(X = mkx(N, w, lam, mu, e1, e2), Y = mky(N, w, lam, mu, e1, e2));
  \\ sum_j log(1 + x_j): coefficientwise power-sum substitution
  my(LX = lg1(X, N), LY = lg1(-Y, N));
  my(SX = sum(n = 0, N - 1, topow(polcoef(LX, n, s), w, pv) * s^n), SY = sum(n = 0, N - 1, topow(polcoef(LY, n, s), w, pv) * s^n));
  my(EX = ex1(SX, N), EY = ex1(SY, N));
  if(quad, my(x0 = trunc(s^2 / aa * inv1(1 - s^2, N), N), y0 = trunc(s^2 / aa * inv1(1 - e1*s + e2*s^2, N), N));
    EX = trunc(EX * (1 + x0), N); EY = trunc(EY * (1 - y0), N));
  my(G = t^2 - 1, H = t^2 - e1*t + e2, c0 = 2^(1 - k));
  my(BV = trunc((1 - s^2)^(k - 1) * EX, N), BW = trunc((1 - e1*s + e2*s^2)^(k - 1) * EY, N));
  my(V = sum(n = 0, N - 1, polcoef(BV, n, s) * t^(2*k - 2 - n)), W = sum(n = 0, N - 1, polcoef(BW, n, s) * t^(2*k - 2 - n)));
  my(R1 = (V - c0*G^(k-1)) % H, R2 = (deriv(V, t) - c0*deriv(G^(k-1), t)) % H);
  my(R3 = (W - c0*H^(k-1)) % G, R4 = (deriv(W, t) - c0*deriv(H^(k-1), t)) % G);
  [polcoef(R1, 0, t), polcoef(R1, 1, t), polcoef(R2, 0, t), polcoef(R2, 1, t), polcoef(R3, 0, t), polcoef(R3, 1, t), polcoef(R4, 0, t), polcoef(R4, 1, t)]};
