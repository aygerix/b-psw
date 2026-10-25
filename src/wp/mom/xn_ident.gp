\\ xn line (extra necessary conditions, real cubic Williams families): exact symbolic identity checks.
\\ Normal form f_i = 2t + 1 + t(t+1)(a_i t + c_i); g_i = (f_i - 1)/t, h_i = (f_i + 1)/(t+1); P = prod g, Q = prod h.
\\ V := (F - 1)/(t P), W := (F + 1)/((t + 1) Q) as RATIONAL functions (no family assumption).  All members symbolic.
\\ Checks (k = 3, 5, 7; series with polynomial coefficients, no rational-function simplification):
\\  (1) V(-1) = 2^(1-k), V'(-1) = (1-k) 2^(1-k)                                 [k odd]
\\  (2) V''(-1) = (k-1)(k-2) 2^(1-k) + 2^(-k) e_2(beta),  beta_i = 2 + a_i - c_i = h_i(-1) = f_i'(-1)
\\      [equivalently 2^(3-k)[(sum x')^2 - sum x'^2] with x_i' = x_i'(-1) = -(2 + a_i - c_i)/4]
\\  (3) V(0) = sum_i 1/prod_{j != i} (c_j + 2)
\\  (6) W(0) = 2^(1-k), W'(0) = (k-1) 2^(1-k), W''(0) = (k-1)(k-2) 2^(1-k) + 2^(-k) e_2(gamma), gamma_i = 2 + c_i = g_i(0)
\\  (7) even-subset identity  2^(k-1) t P V = sum_{|S| even} prod_{i in S} (t+1) h_i prod_{i not in S} t g_i  (polynomial identity)
\\      hence V = 2^(1-k) sum_j t^(k-1-2j) (t+1)^(2j) e_(2j)(h/g), every term >= 0 on R.
\\  (8) V'''(-1)/6 = -C(k-1,3) 2^(1-k) - (k-1) 2^(-1-k) e2(beta) + 2^(1-k) [p(e1) p(e2) - p(e1 e2)],
\\      e1_i = -beta_i/2, e2_i = beta_i + gamma_i/2 - 3 - beta_i^2/4  (Taylor coefficients of eps_i = (f_i+1)/(f_i-1) at t = -1)
t; s; z;
esym(v, j) = polcoef(prod(i = 1, #v, 1 + v[i]*z), j, z);
{
  foreach([3, 5, 7], k,
    my(A = vector(k, i, eval(Str("a", i))), C = vector(k, i, eval(Str("c", i))));
    my(f = vector(k, i, 2*t + 1 + t*(t+1)*(A[i]*t + C[i])), F = prod(i = 1, k, f[i]));
    my(g = vector(k, i, (f[i] - 1)/t), h = vector(k, i, (f[i] + 1)/(t + 1)), P = prod(i = 1, k, g[i]), Q = prod(i = 1, k, h[i]));
    my(Nm = (F - 1)/t, Np = (F + 1)/(t + 1));       \\ exact polynomial divisions
    my(c0 = 2^(1-k), ord = 4);
    \\ V around t = -1: tau = t + 1
    my(Nt = subst(Nm, t, t - 1), Pt = subst(P, t, t - 1));
    my(sV = (Nt + O(t^ord)) / (Pt + O(t^ord)));
    my(beta = vector(k, i, 2 + A[i] - C[i]), gam = vector(k, i, 2 + C[i]));
    print("k = ", k);
    print("  (1) V(-1) - 2^(1-k) = ", polcoef(sV, 0) - c0, ";  V'(-1) - (1-k)2^(1-k) = ", polcoef(sV, 1) - (1-k)*c0);
    print("  (2) V''(-1) - [(k-1)(k-2)2^(1-k) + 2^(-k) e2(beta)] = ", 2*polcoef(sV, 2) - ((k-1)*(k-2)*c0 + 2^(-k)*esym(beta, 2)));
    print("  (3) V(0) - sum 1/prod_(j!=i)(c_j+2) = ", simplify(subst(Nm, t, 0)/subst(P, t, 0) - sum(i = 1, k, 1/prod(j = 1, k, if(j == i, 1, C[j] + 2)))));
    my(E1 = vector(k, i, -beta[i]/2), E2 = vector(k, i, beta[i] + gam[i]/2 - 3 - beta[i]^2/4));
    my(rhs8 = -binomial(k-1, 3)*c0 - (k-1)*2^(-1-k)*esym(beta, 2) + c0*(sum(i = 1, k, E1[i])*sum(i = 1, k, E2[i]) - sum(i = 1, k, E1[i]*E2[i])));
    print("  (8) V'''(-1)/6 - formula = ", polcoef(sV, 3) - rhs8);
    my(sW = (Np + O(t^ord)) / (Q + O(t^ord)));
    print("  (6) W(0) - c0 = ", polcoef(sW, 0) - c0, ";  W'(0) - (k-1)c0 = ", polcoef(sW, 1) - (k-1)*c0, ";  W''(0) - formula = ", 2*polcoef(sW, 2) - ((k-1)*(k-2)*c0 + 2^(-k)*esym(gam, 2)));
    my(B = vector(k, i, (t + 1)*h[i]), Aa = vector(k, i, t*g[i]));
    my(ev = sum(r = 0, 2^k - 1, my(bits = binary(r + 2^k)); my(sz = sum(i = 2, k + 1, bits[i])); if(sz % 2, 0, prod(i = 1, k, if(bits[i + 1], B[i], Aa[i])))));
    print("  (7) even-subset identity residual = ", 2^(k-1)*(F - 1) - ev));
}
