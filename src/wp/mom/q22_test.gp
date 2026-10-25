\\ sanity check of q22_moments.gp against a direct construction: random rational pencil, k = 5 members, compute
\\ V, W directly as polynomial parts and compare the formal power-sum version evaluated at the actual power sums.
read("q22_moments.gp");
{
  my(k = 5, r1 = 1/3, r2 = -1/5, E1 = r1 + r2, E2 = r1*r2, cs = [2, 3, 7/2, 5, 11/3]);
  my(Dl0 = (1 - r1^2)*(1 - r2^2), B = bnd(k, E1, E2, 2*E1/Dl0, 2*(1 - E1^2 + E2)/Dl0));
  my(vals = vector(#B, i, my(b = B[i]);
     for(a = 1, 2*k, b = subst(b, eval(Str("p", a)), sum(j = 1, k, (1/cs[j])^a))); b));
  \\ direct: same remainders computed from exact polynomial division of F -+ 1
  my(G = t^2 - 1, H = (t - r1)*(t - r2), Dl = (1 - r1^2)*(1 - r2^2), lam = 2*E1/Dl, mu = 2*(1 - E1^2 + E2)/Dl);
  my(f0 = 1 + G*(lam*t + mu), f = vector(k, j, f0 + cs[j]*G*H), F = prod(j = 1, k, f[j]));
  my(Km = prod(j = 1, k, (f[j] - 1)/G), Kp = prod(j = 1, k, (f[j] + 1)/H));
  \\ polynomial parts: V = quotient of (F - 1) by G Km (as Euclidean quotient), W likewise
  my(V = (F - 1) \ (G*Km), W = (F + 1) \ (H*Kp), c0 = 2^(1-k));
  my(R1 = (V - c0*G^(k-1)) % H, R2 = (deriv(V, t) - c0*deriv(G^(k-1), t)) % H, R3 = (W - c0*H^(k-1)) % G, R4 = (deriv(W, t) - c0*deriv(H^(k-1), t)) % G);
  my(dir = [polcoef(R1,0,t), polcoef(R1,1,t), polcoef(R2,0,t), polcoef(R2,1,t), polcoef(R3,0,t), polcoef(R3,1,t), polcoef(R4,0,t), polcoef(R4,1,t)]);
  print("formal - direct: ", vals - dir);
}
