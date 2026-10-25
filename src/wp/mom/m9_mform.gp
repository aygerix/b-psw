\\ m9 line: explicit monomial-moment form of the four boundary equations (k = 9, 11) and EXACT identity check
\\ against eqs(m, lin) of cubic_moments.gp (residual 0).  Library: m9_lib.gp.
read("m9_lib.gp");
{
  foreach([[9, 0], [8, 1], [11, 0], [10, 1]], c, my(m = c[1], lin = c[2], k = m + lin, r = mform(m, lin), Z = r[1], E = eqs(m, lin), c0 = 2^(1-k));
    \\ E1 = Z0, E2 = -(k-1) Z0 - Z1 (see m9_moments.gp); reflected likewise
    my(Zr = [reflect(Z[1]), reflect(Z[2])]);
    my(res = [E[1] - evalmom(Z[1], m), E[2] - evalmom(-(k-1)*Z[1] - Z[2], m), E[3] - evalmom(Zr[1], m), E[4] - evalmom(-(k-1)*Zr[1] - Zr[2], m)]);
    print("==== k = ", k, ", cubic members m = ", m, ", linear member: ", lin);
    print("identity residuals (moment form vs eqs): ", res);
    print("moment symbols: ", syms(Z[1] + Z[2]));
    print("Lam_n (moment form):"); for(n = 3, k - 1, print("  Lam_", n, " = ", r[2][n]));
    my(Zp = [(Z[1] + Zr[1])/2, (Z[2] + Zr[2])/2], Zm = [(Z[1] - Zr[1])/2, (Z[2] - Zr[2])/2]);
    print("Z0+ = ", Zp[1]); print("Z1+ = ", Zp[2]); print("Z0- = ", Zm[1]); print("Z1- = ", Zm[2]));
}
