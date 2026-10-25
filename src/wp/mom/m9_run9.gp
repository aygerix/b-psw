\\ driver: search exact feasible points for k = 9 (both cases), then verify exactly.
read("m9_exact9.gp");
{
  setrand(20261025);
  foreach([[9, 0], [8, 1]], c, my(m = c[1], lin = c[2], ZZ = Zs(m, lin), z0m = ZZ[3], z1m = ZZ[4], a0, a1, Phi, phi, found = 0, tries = 0);
    print("==== k = 9, m = ", m, ", lin = ", lin);
    \\ S enters the odd equations only through M1_0 * M1_1 (check) -> S-free combination
    a0 = simplify(polcoef(z0m, 1, M1_0)); a1 = simplify(polcoef(z1m, 1, M1_0));
    if(poldegree(z0m, M1_0) != 1 || poldegree(z1m, M1_0) != 1 || variables(a0) != [M1_1] || variables(a1) != [M1_1], error("odd structure"));
    Phi = simplify((a1/M1_1) * z0m - (a0/M1_1) * z1m);
    if(polcoef(Phi, 1, M1_0) != 0, error("S not eliminated"));
    print("odd equations: Z0- = ", z0m); print("               Z1- = ", z1m);
    print("S-free combination Phi = ", Phi);
    foreach(syms(Phi), z, if(Vec(Str(z))[#Vec(Str(z))] == "0", error("Phi involves an even moment")));
    phi = msub(Phi, [[ww, vv]], []);
    print("single-member phi(w, v) = ", factor(phi));
    my(cfg = [[3, 1, 1], [3, 0, 3], [2, 1, 2], [4, 1, 0], [2, 2, 0], [5, 0, 1], [6, 0, 0], [1, 1, 3], [1, 2, 1], [3, 1, 0], [2, 1, 1], [4, 0, 1], [1, 2, 0], [2, 0, 3], [5, 0, 0]]);
    while(found < 3 && tries < 400000,
      tries++;
      my(cf = cfg[1 + random(#cfg)], r);
      if(cf[1] + 2*cf[2] + cf[3] + 3 != m, next);
      r = try9(m, lin, ZZ, phi, cf[1], cf[2], cf[3]);
      if(r != 0, found++;
        print("FOUND after ", tries, " tries, configuration [nA, nP, nBf] = ", cf);
        print("  members [w, v] (rational part): ", r[1]);
        print("  group-B cubic C(x) = ", r[2]);
        report9(r[1], r[2]);
        print("  EXACT check of the four boundary equations (eqs, brute force) in the splitting algebra: ", verify9(m, lin, r[1], r[2]))));
    print("tries: ", tries, ", found: ", found));
}
