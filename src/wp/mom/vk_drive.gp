\\ vk_drive.gp: echo 'SOLFILE="m9_sol_9_0.gp"; read("vk_drive.gp")' | gp -q
read("cubic_moments.gp");
read("vk_indep.gp");
read(SOLFILE);
{
  my(m = CASE[1], lin = CASE[2], n = 2*m, E, ALL, bad = 0, Jn, dets = List(), best, sel1, sel2, X0 = concat(U0, V0));
  print("==== vk independent check: ", SOLFILE, "  m = ", m, ", lin = ", lin, ", k = ", m + lin);
  E = eqs(m, lin);
  ALL = concat(vector(m, i, eval(Str("u", i))), vector(m, i, eval(Str("v", i))));
  \\ (2) exact comparison at random rational points
  setrand(20261025);
  for(tr = 1, 6,
    my(P = vector(n, j, (random(2*10^6) - 10^6) / (1 + random(10^4))));
    if(substvec(E, ALL, P) != allbnd(P[1..m], P[m+1..n], lin), bad++));
  print(" (2) eqs (Laurent truncation) vs exact division quotient at 6 random rational points: ",
        if(bad, Str(bad, " MISMATCHES"), "all EQUAL"));
  \\ numerical Jacobian and minor selection
  Jn = numjac(U0, V0, lin);
  forsubset([n, 4], s, listput(dets, [abs(matdet(vecextract(Jn, "1..4", Vec(s)))), Vec(s)]));
  dets = vecsort(Vec(dets), 1, 4);
  sel1 = dets[1][2];
  sel2 = 0; for(j = 2, #dets, if(#setminus(Set(dets[j][2]), Set(sel1)) >= 2, sel2 = dets[j][2]; break));
  print(" numerical residual of the stored solution (max |E|): ", normlp(allbnd(U0, V0, lin), oo));
  print(" best slice: ", vector(4, j, ALL[sel1[j]]), " |det| = ", 1.0*dets[1][1]);
  print(" alternative slice: ", vector(4, j, ALL[sel2[j]]));
  \\ (2b) symbolic identity on the m9 slice (fixed coordinates rounded with 10^8, unknowns = sel1, as in m9_kraw)
  my(fix = setminus([1..n], Set(sel1)), Fv = vector(#fix, j, ALL[fix[j]]), Fq = vector(#fix, j, rat(X0[fix[j]], 10^8)),
     Fthem = vector(4, r, substvec(substvec(E[r], Fv, Fq), vector(4, j, ALL[sel1[j]]), [y1, y2, y3, y4])),
     Xs = vector(n, c, rat(X0[c], 10^8)));
  for(j = 1, 4, Xs[sel1[j]] = [y1, y2, y3, y4][j]);
  my(Fme = allbnd(Xs[1..m], Xs[m+1..n], lin));
  print(" (2b) m9 slice polynomials, eqs vs division quotient: ", if(Fthem == Fme, "IDENTICAL", "DIFFERENT"));
  \\ (3) own Krawczyk certificates
  print(" (3a) own certificate on best slice, denominator 2^30, rho = 1e-70:");
  my(r1 = vkraw(U0, V0, lin, sel1, 2^30, 1/10^70));
  print("   VERDICT 3a: ", if(r1[1], "EXISTS", "NOT CERTIFIED"));
  print(" (3b) own certificate on alternative slice, denominator 3^19, rho = 1e-70:");
  my(r2 = vkraw(U0, V0, lin, sel2, 3^19, 1/10^70));
  print("   VERDICT 3b: ", if(r2[1], "EXISTS", "NOT CERTIFIED"));
}
