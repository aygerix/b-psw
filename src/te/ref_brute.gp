\\ ref_brute.gp -- independent brute-force reference for the Theorem E engine (see SPEC.md in this directory).
\\
\\ For a Selfridge-type discriminant D (D == 1 mod 4, P = 1, Q = (1-D)/4) and a bound X, list the inert primes
\\ p <= X (kronecker(D,p) = -1, p not dividing 2QD) and test EVERY triple p1 < p2 < p3, n = p1*p2*p3 (S = 1):
\\   fermat : n is a strong probable prime to base 2;
\\   lucas  : n is a strong Lucas probable prime for the FIXED (P,Q) = (1,(1-D)/4) (not Selfridge);
\\   bpsw   : n passes the exact BPSW test of src/bpsw.py (strong Fermat base 2 + Selfridge Method A strong Lucas);
\\   weak   : lcm(ord(p), z(p)) | n/p - 1 for all three p | n.
\\ plus two one-torus versions of weak (cheap extras, non-empty unlike weak at these X; weak = weakord cap weakz):
\\   weakord: ord(p) | n/p - 1 for all three p   (necessary for fermat: the W1 survivors of --mode fermat);
\\   weakz  : z(p) | n/p - 1 for all three p     (necessary for lucas: the W1 survivors of --mode lucas).
\\ Nothing from the engine's theory is used: ord and z are computed by naive iteration, the Lucas test uses the
\\ 2x2 matrix power (a different algorithm from bpsw.py and from the engine), and every predicate is applied
\\ literally to n.
\\
\\ Two run modes:
\\   full = 1 : every predicate is evaluated on every triple (no filter at all).
\\   full = 0 : "filtered" -- the expensive tests are only run on triples passing a trivially necessary condition:
\\              strong Fermat base 2  => 2^(n-1) == 1 mod n  => ord(p) | n-1 for each p | n;
\\              strong Lucas (D/n=-1) => U_(n+1) == 0 mod n  => z(p) | n+1 for each p | n;
\\              exact BPSW            => strong Fermat base 2, so bpsw is only evaluated on the fermat list.
\\ Validation done (2026-10-24): ref_selftest(10^6): bpsw == isprime for all odd n < 1e6, and the composite
\\ strong-Lucas(Selfridge) list (58 numbers) and the 2-spsp list (46) below 1e6 are identical to those of src/bpsw.py;
\\ lucasmat == literal port of bpsw.py's lucas_UV on 20000 random inputs; full and filtered modes give identical
\\ lists at X = 2000, 6000 for D = 5, -7, 13; a third (Python, bpsw.py-based) brute force agrees at X = 2000.
\\ Reference outputs: ref_out/ (X = 2000, 6000 in full mode; X = 20000 in filtered mode). The weakord/weakz files
\\ come from a filtered-mode re-run (those two predicates are evaluated on every triple in both modes); that re-run
\\ reproduced all fermat/lucas/bpsw/weak files exactly. The X = 2000, 6000 lists equal the p3 <= X part of X = 20000.
\\
\\ Usage: gp -q ref_brute.gp   then   ref_run(D, X, outdir, full)   (or ref_selftest())
\\    or: REF_D=5 REF_X=2000 REF_FULL=1 REF_OUT=ref_out gp -q ref_brute.gp     (runs and quits)
\\ Output files: outdir/D<D>_X<X>_<pred>.txt with '#' header lines and data lines "n p1 p2 p3", sorted by n.

\\ ---------- primality-test building blocks (exact, n odd) ----------

\\ strong probable prime base 2 (same definition as strong_fermat2 in src/bpsw.py)
sfermat2(n) =
{
  my(d = n - 1, s = valuation(n - 1, 2), x);
  d = d >> s;
  x = Mod(2, n)^d;
  if (x == 1 || x == -1, return(1));
  for (i = 1, s - 1, x = x^2; if (x == -1, return(1)));
  0;
}

\\ [U_k, V_k, Q^k] mod n for U(P,Q) via the matrix [[P,-Q],[1,0]]^k = [[U_(k+1), -Q U_k],[U_k, -Q U_(k-1)]]
\\ and V_k = 2 U_(k+1) - P U_k.
lucasmat(k, P, Q, n) =
{
  my(A = Mod([P, -Q; 1, 0], n)^k, U = A[2,1], U1 = A[1,1]);
  [lift(U), lift(2*U1 - P*U), lift(Mod(Q, n)^k)];
}

\\ literal port of lucas_UV from src/bpsw.py (used only for cross-checking lucasmat in ref_selftest)
lucasUV_port(k, P, Q, n) =
{
  my(D = P*P - 4*Q, U = 0, V = 2, Qk = 1, inv2 = (n + 1)/2, b = binary(k), t);
  for (i = 1, #b,
    t = U*V % n; V = (V*V - 2*Qk) % n; U = t;
    Qk = Qk*Qk % n;
    if (b[i],
      t = (P*U + V)*inv2 % n; V = (D*U + P*V)*inv2 % n; U = t;
      Qk = Qk*Q % n));
  [U, V, Qk];
}

\\ strong Lucas test with P = 1 and the given D (Q = (1-D)/4), n odd, using n - (D/n)
slucasD(n, D) =
{
  my(Q = (1 - D)/4, eps = kronecker(D, n), e, t, W, U, V, Qk);
  if (eps == 0, return(0));
  e = n - eps; t = valuation(e, 2); e = e >> t;
  W = lucasmat(e, 1, Q, n); U = W[1]; V = W[2]; Qk = W[3];
  if (U == 0 || V == 0, return(1));
  for (i = 1, t - 1,
    V = (V*V - 2*Qk) % n; Qk = Qk*Qk % n;
    if (V == 0, return(1)));
  0;
}

\\ Selfridge Method A exactly as selfridge_D in src/bpsw.py: perfect squares -> none (0 here);
\\ D = 5, -7, 9, -11, 13, ...: first D with Jacobi(D, n) = -1 (Jacobi 0 is skipped, i.e. the search continues)
selfridgeD(n) =
{
  my(D = 5);
  if (issquare(n), return(0));
  while (kronecker(D, n) != -1, D = if (D > 0, -D - 2, -D + 2));
  D;
}

slucas_selfridge(n) =
{
  my(D = selfridgeD(n));
  if (D == 0, return(0));
  slucasD(n, D);     \\ kronecker(D,n) = -1 here, so this uses n + 1 as in bpsw.py
}

bpsw(n) = n > 1 && n % 2 == 1 && sfermat2(n) && slucas_selfridge(n);

\\ ---------- per-prime data by naive iteration ----------

ord2_naive(p) = { my(x = 2 % p, k = 1); while (x != 1, x = 2*x % p; k++); k; }

\\ least k > 0 with p | U_k(1, Q): iterate U_(k+1) = U_k - Q U_(k-1)
zrank_naive(p, Q) =
{
  my(u0 = 0, u1 = 1, k = 1, t);
  while (u1 % p, t = (u1 - Q*u0) % p; u0 = u1; u1 = t; k++);
  k;
}

inert_primes(D, X) =
{
  my(Q = (1 - D)/4, L = List());
  forprime (p = 3, X, if (kronecker(D, p) == -1 && (2*Q*D) % p != 0, listput(L, p)));
  Vec(L);
}

\\ ---------- output ----------

writelist(fname, hdr, L) =
{
  my(f = fileopen(fname, "w"), v = vecsort(Vec(L), 1));
  filewrite(f, hdr);
  filewrite(f, Str("# count=", #v));
  for (i = 1, #v, filewrite(f, Str(v[i][1], " ", v[i][2], " ", v[i][3], " ", v[i][4])));
  fileclose(f);
  #v;
}

\\ ---------- the brute force ----------

ref_run(D, X, outdir, full) =
{
  my(Q = (1 - D)/4, P, N, O, Z, M, LF = List(), LL = List(), LB = List(), LW = List(), LWo = List(), LWz = List(),
     t0 = getwalltime(), p1, p2, p3, m, n, nt = 0, hdr, nf = 0, nl = 0);
  if (D % 4 != 1, error("D must be 1 mod 4"));
  P = inert_primes(D, X); N = #P;
  O = vector(N, i, ord2_naive(P[i]));
  Z = vector(N, i, zrank_naive(P[i], Q));
  for (i = 1, N,
    if (O[i] != znorder(Mod(2, P[i])), error("ord mismatch at ", P[i]));
    if ((P[i] - 1) % O[i], error("ord does not divide p-1 at ", P[i]));
    if ((P[i] + 1) % Z[i], error("z does not divide p+1 at ", P[i]));
    if (lucasmat(Z[i], 1, Q, P[i])[1] != 0, error("U_z != 0 at ", P[i])));
  M = vector(N, i, lcm(O[i], Z[i]));
  for (i = 1, N - 2, p1 = P[i];
    for (j = i + 1, N - 1, p2 = P[j]; m = p1*p2;
      for (k = j + 1, N, p3 = P[k]; n = m*p3; nt++;
        if (full,
          \\ every predicate on every triple
          if (sfermat2(n), listput(LF, [n, p1, p2, p3]));
          if (bpsw(n), listput(LB, [n, p1, p2, p3]));
          if (slucasD(n, D), listput(LL, [n, p1, p2, p3]))
        ,
          \\ trivially necessary prefilters, then the literal tests
          if ((n - 1) % O[k] == 0 && (n - 1) % O[j] == 0 && (n - 1) % O[i] == 0, nf++;
            if (sfermat2(n), listput(LF, [n, p1, p2, p3]);
              if (slucas_selfridge(n), listput(LB, [n, p1, p2, p3]))));
          if ((n + 1) % Z[k] == 0 && (n + 1) % Z[j] == 0 && (n + 1) % Z[i] == 0, nl++;
            if (slucasD(n, D), listput(LL, [n, p1, p2, p3])))
        );
        if ((p2*p3 - 1) % M[i] == 0 && (p1*p3 - 1) % M[j] == 0 && (m - 1) % M[k] == 0,
          listput(LW, [n, p1, p2, p3]));
        if ((p2*p3 - 1) % O[i] == 0 && (p1*p3 - 1) % O[j] == 0 && (m - 1) % O[k] == 0,
          listput(LWo, [n, p1, p2, p3]));
        if ((p2*p3 - 1) % Z[i] == 0 && (p1*p3 - 1) % Z[j] == 0 && (m - 1) % Z[k] == 0,
          listput(LWz, [n, p1, p2, p3])))));
  hdr = Str("# ref_brute.gp D=", D, " P=1 Q=", Q, " X=", X, " S=1 mode=", if (full, "full", "filtered"),
            " inert_primes=", N, " triples=", nt);
  my(cf = writelist(Str(outdir, "/D", D, "_X", X, "_fermat.txt"), Str(hdr, " pred=fermat (strong 2-prp)"), LF),
     cl = writelist(Str(outdir, "/D", D, "_X", X, "_lucas.txt"), Str(hdr, " pred=lucas (strong Lucas, fixed D)"), LL),
     cb = writelist(Str(outdir, "/D", D, "_X", X, "_bpsw.txt"), Str(hdr, " pred=bpsw (exact, Selfridge A)"), LB),
     cw = writelist(Str(outdir, "/D", D, "_X", X, "_weak.txt"), Str(hdr, " pred=weak (lcm(ord,z) | n/p-1 at all p)"), LW),
     cwo = writelist(Str(outdir, "/D", D, "_X", X, "_weakord.txt"), Str(hdr, " pred=weakord (ord | n/p-1 at all p)"), LWo),
     cwz = writelist(Str(outdir, "/D", D, "_X", X, "_weakz.txt"), Str(hdr, " pred=weakz (z | n/p-1 at all p)"), LWz));
  print("D=", D, " X=", X, " mode=", if (full, "full", "filtered"), " inert=", N, " triples=", nt,
        " fermat=", cf, " lucas=", cl, " bpsw=", cb, " weak=", cw, " weakord=", cwo, " weakz=", cwz,
        if (full, "", Str(" prefilter_fermat=", nf, " prefilter_lucas=", nl)),
        " time=", strtime(getwalltime() - t0));
  [cf, cl, cb, cw, cwo, cwz];
}

\\ ---------- self-test ----------
\\ (1) lucasmat == literal bpsw.py port on many (k, D, n);
\\ (2) bpsw(n) == isprime(n) for all odd 3 <= n < B (no BPSW pseudoprime below 2^64 is known/exists);
\\ (3) writes the composite strong-Lucas-Selfridge and 2-spsp lists below B to outdir for comparison with bpsw.py.
ref_selftest(B, outdir) =
{
  my(bad = 0, Ls = List(), Lf = List(), f, n, D, Q, k, a, b);
  for (i = 1, 20000,
    n = 2*random(10^12) + 3; D = [5, -7, 13, -11, 17, -15, 21][random(7) + 1]; Q = (1 - D)/4;
    k = random(10^15) + 1;
    a = lucasmat(k, 1, Q, n); b = lucasUV_port(k, 1, Q, n);
    if (a != b, bad++; print("lucas mismatch k=", k, " Q=", Q, " n=", n)));
  forstep (n = 3, B, 2,
    if (bpsw(n) != isprime(n), bad++; print("bpsw != isprime at ", n));
    if (!isprime(n),
      if (slucas_selfridge(n), listput(Ls, n));
      if (sfermat2(n), listput(Lf, n))));
  f = fileopen(Str(outdir, "/selftest_slucas_selfridge_composites.txt"), "w");
  for (i = 1, #Ls, filewrite(f, Str(Ls[i]))); fileclose(f);
  f = fileopen(Str(outdir, "/selftest_spsp2_composites.txt"), "w");
  for (i = 1, #Lf, filewrite(f, Str(Lf[i]))); fileclose(f);
  print("selftest: B=", B, " mismatches=", bad, " composite strong-Lucas(Selfridge)=", #Ls,
        " (first: ", Vec(Ls)[1..min(8, #Ls)], ") composite 2-spsp=", #Lf, " (first: ", Vec(Lf)[1..min(8, #Lf)], ")");
  bad;
}

{
  if (getenv("REF_SELFTEST"),
    ref_selftest(eval(getenv("REF_SELFTEST")), getenv("REF_OUT")); quit);
  if (getenv("REF_D"),
    ref_run(eval(getenv("REF_D")), eval(getenv("REF_X")), getenv("REF_OUT"), eval(getenv("REF_FULL"))); quit);
}
