read("q22_moments.gp");
{
  my(k = 7, B = bnd(k, e1, e2, lam, mu));
  my(Bs = vector(8, i, subst(subst(B[i], e1, 0), lam, 0)));
  \\ substitute e2 = ep - 1, mu = 2/ep, and p_a = (ep/2)^a q_a  (q_a = power sums of mu*w)
  my(sb(P) = my(Q = subst(subst(P, e2, ep - 1), mu, 2/ep)); for(a = 1, 9, Q = subst(Q, eval(Str("p", a)), (ep/2)^a * eval(Str("q", a)))); Q);
  my(S1 = sb(Bs[1]), S4 = sb(Bs[4]), S5 = sb(Bs[5]), S8 = sb(Bs[8]));
  print("B1+B5: ", numerator(S1 + S5));
  print("B1-B5: ", numerator(S1 - S5));
  print("B4+B8: ", numerator(S4 + S8));
  print("B4-B8: ", numerator(S4 - S8));
}
