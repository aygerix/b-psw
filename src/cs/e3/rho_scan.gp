\\ src/cs/e3/rho_scan.gp -- E3 confirm branch: end-to-end scan of the rho-family on prime specialisations.
\\ For rho = a/b and x = X/d, keep X with f1, f2, f3 all integers > 3, distinct and (pseudo)prime; then
\\   n = p1 p2 p3:  Selfridge D, Jacobi (D/p_j), sprp2(n), slprp(n), Korselt (p_j - 1 | n - 1),
\\   Fermat defect dF = lcm odd(ord_p 2) / gcd(., n - 1), Lucas defect dL = lcm odd(z(p)) / gcd(., n + 1),
\\   per-member Lucas defect dL_j = odd(z_j) / gcd(odd(z_j), n + 1).
\\ Every n with all three factors prime is written to OUTN for an independent check with src/bpsw.py.
\\ Usage: gp -q rho_scan.gp with globals AMAX, BMAX, DLIST, XMAX, MMAX, OUTN set by the caller (see run line).
read("scan_core.gp");
fam(r) = my(B = 2*x^2 - (4*r+1)*x + 2*r - 2); [1 + B, 1 + B*x/(2*r), 1 + B*(x-1)/(2*r-1)];
scanfam(a, b, d, Xmax, Mmax, fh) = scanpolys(fam(a/b), Str(a, "/", b), d, Xmax, Mmax, fh);
{
  my(fh = fileopen(OUTN, "w"), nf = 0, t0 = getabstime());
  for (b = 1, BMAX, for (a = -AMAX, AMAX, if (a == 0 || gcd(a, b) != 1 || 2*a == b, next);
    foreach (DLIST, d, my(c = scanfam(a, b, d, XMAX, MMAX, fh)); nf++;
      if (c < 0, print("  rho=", a, "/", b, " d=", d, ": modulus > MMAX, skipped")))));
  fileclose(fh);
  print("families x denominators scanned: ", nf, "  time(ms) ", getabstime() - t0);
  printstats();
}
