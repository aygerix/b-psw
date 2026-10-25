\\ m9 line: the V'(-1) equation (E2 = -( (k-1) Z0 + Z1 )) split into even/odd parts, k = 9, 11 (both cases):
\\ the only nonlinear terms left are low products; prints them (these are the "forced" relations).
read("m9_lib.gp");
{
  foreach([[9, 0], [8, 1], [11, 0], [10, 1]], c, my(m = c[1], lin = c[2], k = m + lin, r = mform(m, lin), Z = r[1], Zr = [reflect(Z[1]), reflect(Z[2])]);
    my(Zp = [(Z[1] + Zr[1])/2, (Z[2] + Zr[2])/2], Zm = [(Z[1] - Zr[1])/2, (Z[2] - Zr[2])/2]);
    print("==== k = ", k, ", m = ", m, ", lin = ", lin);
    print("  even part of (k-1) Z0 + Z1 [= -(V'(-1) - (1-k) 2^(1-k)) even]: ", simplify((k-1)*Zp[1] + Zp[2]));
    print("  odd  part of (k-1) Z0 + Z1:                                   ", simplify((k-1)*Zm[1] + Zm[2])));
}
