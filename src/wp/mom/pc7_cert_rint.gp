\\ reality interval lemma (symbolic in r1, r2): with q1(w) = (w lam - e1)^2 - 4(e2 + w mu) [K-_j = c_j H + L root-free iff q1(w_j) < 0]
\\ and q2(w) = w^2 lam^2 - 4 w (mu + lam e1) + 4 [K+_j = c_j G + L' root-free iff q2(w_j) < 0]:
\\   (i) q2(w) = q1(w - Dl/2);  (ii) disc_w(q1) = 16 L(r1) L(r2) = 64/Dl;  hence I = {q1(w) < 0, q1(w - Dl/2) < 0} is the interval
\\   (u_- + Dl/2, u_+) between the roots u_- < u_+ of q1, nonempty iff u_+ - u_- = 8/(lam^2 sqrt(Dl)) > Dl/2 iff e1^4 < 16 Dl  (e1 != 0).
r1; r2; ww;
{
  my(E1 = r1 + r2, E2 = r1*r2, Dl = (1 - r1^2)*(1 - r2^2), LA = 2*E1/Dl, MU = 2*(1 - E1^2 + E2)/Dl);
  my(q1 = (ww*LA - E1)^2 - 4*(E2 + ww*MU), q2 = ww^2*LA^2 - 4*ww*(MU + LA*E1) + 4);
  print("(i)  q2(w) - q1(w - Dl/2) = ", simplify(q2 - subst(q1, ww, ww - Dl/2)));
  print("(ii) disc(q1) - 64/Dl = ", simplify(poldisc(q1, ww) - 64/Dl), ";  L(r1) L(r2) - 4/Dl = ", simplify((LA*r1 + MU)*(LA*r2 + MU) - 4/Dl));
  \\ width^2 = disc/lam^4 versus (Dl/2)^2:  width > Dl/2  <=>  64/(Dl lam^4) > Dl^2/4  <=>  256 > Dl^3 lam^4 = 16 e1^4 / Dl
  print("(iii) Dl^3 lam^4 - 16 e1^4/Dl = ", simplify(Dl^3*LA^4 - 16*E1^4/Dl));
}
