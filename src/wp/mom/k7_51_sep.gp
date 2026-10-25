\\ exact separator for k = 7, 5 cubics + 2t+1 + one invisible member (degree >= 7): forced identities as for 6 cubics + 2t+1
\\ (k7_mixed.gp shows the moment form is m-independent), with M < (S - (m-1)/32)^2 + (m-1)/1024 for m = 5.
{
  my(m = 5, Kf = 1355/32 - 62*S - 2*S^2, Mf = 623/32 - 28*S - S^2, B = (S - (m-1)/32)^2 + (m-1)/1024 - Mf);
  print("roots K: ", polrootsreal(Kf), "  roots Mmax-M: ", polrootsreal(B));
  my(r = 3343/5000); print("r = ", r, "  K(r) = ", subst(Kf, S, r), "  (Mmax-M)(r) = ", subst(B, S, r),
    "  K' = ", deriv(Kf, S), "  B' = ", deriv(B, S), " at m/32: ", subst(deriv(B, S), S, m/32));
}
