\\ pc7_cert_dec.gp -- 3-body decomposition of the k = 7 (2,2) pencil boundary functionals.
\\ V = [G^6 prod(1 + x_j)]_+ = G^6 + sum_j w_j (G^5 \ Q_j) + sum_{i<j} w_i w_j (G^4 \ (Q_i Q_j)) + e3(w),  Q_j = H + w_j L
\\ W = [H^6 prod(1 - y_j)]_+ = H^6 - sum_j w_j (H^5 \ P_j) + sum_{i<j} w_i w_j (H^4 \ (P_i P_j)) - e3(w),  P_j = G + w_j L'
\\ (x_j = w_j/(G Q_j), y_j = w_j/(H P_j); e_m(x) for m >= 4 has no polynomial part.)
t; w; v; u; e1; e2; lam; mu;
\\ shape data with lam, mu symbolic (relations lam Dl = 2 e1, mu Dl = 2(1 - e1^2 + e2) imposed afterwards)
G = t^2 - 1; H = t^2 - e1*t + e2; L = lam*t + mu; Lp = lam*t + mu + lam*e1;
c0 = 2^(-6);
Vof(ws) = {my(k = #ws, Q = vector(k, j, H + ws[j]*L), V = G^6);
  for(j = 1, k, V += ws[j] * (G^5 \ Q[j]));
  for(i = 1, k, for(j = i + 1, k, V += ws[i]*ws[j] * (G^4 \ (Q[i]*Q[j]))));
  for(i = 1, k, for(j = i + 1, k, for(l = j + 1, k, V += ws[i]*ws[j]*ws[l])));
  V};
Wof(ws) = {my(k = #ws, P = vector(k, j, G + ws[j]*Lp), W = H^6);
  for(j = 1, k, W -= ws[j] * (H^5 \ P[j]));
  for(i = 1, k, for(j = i + 1, k, W += ws[i]*ws[j] * (H^4 \ (P[i]*P[j]))));
  for(i = 1, k, for(j = i + 1, k, for(l = j + 1, k, W -= ws[i]*ws[j]*ws[l])));
  W};
Bof(ws) = {my(V = Vof(ws), W = Wof(ws));
  my(R1 = (V - c0*G^6) % H, R2 = (deriv(V, t) - c0*deriv(G^6, t)) % H);
  my(R3 = (W - c0*H^6) % G, R4 = (deriv(W, t) - c0*deriv(H^6, t)) % G);
  [polcoef(R1, 0, t), polcoef(R1, 1, t), polcoef(R2, 0, t), polcoef(R2, 1, t), polcoef(R3, 0, t), polcoef(R3, 1, t), polcoef(R4, 0, t), polcoef(R4, 1, t)]};
