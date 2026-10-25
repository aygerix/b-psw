\\ src/cs/e3/p23_scan.gp -- E3: end-to-end scan of the point P23 (I1 = 1, I2 = -4/3, I3 = 3): members 2 and 3 have
\\ gcd degree 2 (defect 1 each), member 1 is dead; total Lucas defect 4 (the minimum, tied with the rho-family).
\\ Globals from the caller: DLIST, XMAX, MMAX, OUTN.
read("scan_core.gp");
{
  my(f = [x^2 - x - 1/3, 3*x^3 - 3*x^2 - 4*x + 1, -3*x^3 + 6*x^2 + x - 3], fh = fileopen(OUTN, "w"), t0 = getabstime());
  foreach (DLIST, d, my(c = scanpolys(f, "P23", d, XMAX, MMAX, fh)); print("  d=", d, ": tuples ", c));
  fileclose(fh); print("time(ms) ", getabstime() - t0); printstats();
}
