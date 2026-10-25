\\ driver: gp -q m9_kcert.gp with SOLFILE set, e.g.  echo 'SOLFILE="m9_sol_9_0.gp"; read("m9_kcert.gp")' | gp -q
\p 160
read("m9_kraw.gp");
read(SOLFILE);
{
  my(m = CASE[1], lin = CASE[2], r);
  print("==== Krawczyk certificate: k = ", m + lin, ", cubic members ", m, ", linear member ", lin, "  (", SOLFILE, ")");
  r = krawczyk(m, lin, U0, V0, 10^8);
  print("  fixed coordinates (rational, denominator 10^8): ", r[4]);
  print("  box centre of the unknowns (rational, 110 digits; shown to 30): ", precision(1.0 * r[3], 30));
  print("  CERTIFICATE ", if(r[1], "VALID: a real admissible family-candidate satisfying all four boundary equations EXISTS", "NOT VALID"));
}
