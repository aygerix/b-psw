\\ mkplants.gp -- planted-control generator for the Theorem E engine, --mode plant (see SPEC.md here).
\\
\\ A plant is a triple p1 < p2 < p3 <= X of inert primes for D (kronecker(D,p) = -1, p not dividing 2QD,
\\ Q = (1-D)/4), a split part S (1, or a product of one or two primes q with kronecker(D,q) = +1, q not dividing
\\ 2QD) and W2 moduli A | n-1, B | n+1 (A = B = 1 when S = 1), n = p1 p2 p3 S. The overrides are
\\     M'(p_i) := gcd(p_i^2 - 1, n/p_i - 1),
\\ so (W1) M'(p_i) | n/p_i - 1 holds at all three primes by construction, and (W2) n == 1 mod A, n == -1 mod B
\\ holds by the choice of A, B. All plants in one file are pairwise disjoint (no prime is reused).
\\
\\ Classes (thresholds as in the task statement):
\\   a  : M'(p3) large:  idx3 := (p3^2-1)/M'(p3) <= 1000                       (branch A at p3)
\\   b  : M'(p3) small:  idx3 >= min(1e8, (p3^2-1)/48);  p2/lcm(M'(p2),M'(p3)) <= 100   (branch C, short loop)
\\   c  : both small:    idx3 >= min(1e8, (p3^2-1)/48);  p2/lcm(M'(p2),M'(p3)) >= min(1000, p2/48)
\\                       (every third c plant is "extreme": lcm(M'(p2),M'(p3)) <= 12, i.e. a very long loop)
\\   d/a, d/b, d/c : the same three constructions with a split part S > 1 and nontrivial A, B.
\\
\\ File format (SPEC.md): "O p M" override lines and "T p1 p2 p3 S A B" lines; each plant is preceded by one
\\ comment line "# plant <k> class=<..> Sfac=<q1*q2> idx3=.. ..." with its metrics. All comment lines start '#'.
\\
\\ Override rule variants: the default (header token mrule=gcd) is M' = gcd(p^2-1, n/p-1) exactly as specified.
\\ About 1/5 of such M' have v2(M') > max(v2(p-1), v2(p+1)), so they are NOT of the form lcm(a, b) with a | p-1,
\\ b | p+1 (a real lcm(ord, z) never is like that). With MKP_LCMFORM=1 (header token mrule=lcmform) the 2-part
\\ is cut down to max(v2(p-1), v2(p+1)); W1 still holds, and every override is a possible lcm(ord', z').
\\
\\ Usage: gp -q mkplants.gp, then  mkplants(D, X, N, seed, fname)  (N plants per class a, b, c, d)
\\        and  checkfile(fname, D, X)  (re-parses the written file and re-verifies every plant).
\\    or: MKP_D=5 MKP_X=10^6 MKP_N=30 MKP_SEED=1 MKP_OUT=plants/x.txt [MKP_LCMFORM=1] gp -q mkplants.gp
\\                                                                              (generate + checkfile, quit)
\\        MKP_CHECK=plants/x.txt MKP_D=5 MKP_X=10^6 gp -q mkplants.gp          (check only)
\\ Independent pure-Python re-check: python3 checkplants.py D X FILE...
\\ Files made 2026-10-24 (plants/): D5_X1e6 (seed 101), D5_X1e7 (102), D-7_X1e6 (103), the same with _lcm
\\ (seeds 201-203, independent plants), D5_X2e4_small (N = 10, seed 301) and D5_X3e3_small (N = 8, seed 302).
\\ Forcing branch A (--cap infinity) on planted p3 with small M' costs about sum p3^2/lcm(M'(p3),A,B) candidates:
\\ ~8e12 (X = 1e6 files), ~7e14 (X = 1e7), 1.2e9 (X2e4_small), 2.9e7 (X3e3_small): use the small files for that control.

D0 = 0; Q0 = 0; X0 = 0; USED = Map(); LCMFORM = 0;

isinert(p) = p > 2 && (2*Q0*D0) % p != 0 && kronecker(D0, p) == -1 && isprime(p);
issplit(q) = q > 2 && (2*Q0*D0) % q != 0 && kronecker(D0, q) == 1 && isprime(q);
isfree(p) = !mapisdefined(USED, p);

rnd(lo, hi) = lo + random(hi - lo + 1);
rndlog(lo, hi) = max(lo, min(hi, floor(exp(log(lo) + random(1.)*(log(hi) - log(lo))))));

shuffle(v) =
{
  my(w = v, j, t);
  forstep (i = #w, 2, -1, j = random(i) + 1; t = w[i]; w[i] = w[j]; w[j] = t);
  w;
}

\\ random unused inert prime in [lo, hi] (uniform or log-uniform start point)
inert_in(lo, hi, lg) =
{
  my(p);
  if (hi < lo, return(0));
  for (t = 1, 10^5,
    p = nextprime(if (lg, rndlog(max(lo, 2), hi), rnd(lo, hi)));
    if (p <= hi && isinert(p) && isfree(p), return(p)));
  0;
}

split_in(lo, hi) =
{
  my(q);
  for (t = 1, 10^5, q = nextprime(rndlog(lo, hi)); if (q <= hi && issplit(q), return(q)));
  error("no split prime in [", lo, ",", hi, "]");
}

\\ override rule: M'(p) = gcd(p^2-1, n/p-1) (default, as in the task), or with LCMFORM = 1 its largest divisor of
\\ "lcm form" (2-part cut down to max(v2(p-1), v2(p+1))), i.e. a value that lcm(a|p-1, b|p+1) can take.
Mg(p, n) = gcd(p^2 - 1, n/p - 1);
lcmform(p, M) = M >> max(0, valuation(M, 2) - max(valuation(p - 1, 2), valuation(p + 1, 2)));
Mp(p, n) = if (LCMFORM, lcmform(p, Mg(p, n)), Mg(p, n));
idxmin(p3) = min(10^8, (p3^2 - 1) \ 48);

\\ an override M is "lcm form" iff M = lcm(a,b) with a | p-1, b | p+1, iff v2(M) <= max(v2(p-1), v2(p+1))
islcmform(p, M) = valuation(M, 2) <= max(valuation(p - 1, 2), valuation(p + 1, 2));

\\ ---------- constructions; each returns [p1,p2,p3] or 0 ----------

\\ class a: d | p3^2-1 with (p3^2-1)/d = k <= 1000; m == S^-1 mod d, m < p3^2, m = p1 p2 inert, p1 < p2 < p3
build_a(S) =
{
  my(p3, N3, ks, k, d, r, J, js, m, f, p1, p2, lg);
  for (t = 1, 20000,
    lg = random(2);
    p3 = inert_in(if (lg, min(1000, X0 \ 10), X0 \ 10), X0, lg); if (!p3, next);
    N3 = p3^2 - 1;
    ks = select(x -> x >= 2 && x <= 1000, divisors(N3)); if (#ks == 0, next);
    k = if (random(4), ks[max(1, #ks - random(min(5, #ks)))], ks[random(#ks) + 1]);
    d = N3 / k;
    if (gcd(d, S) != 1, next);
    r = lift(Mod(S, d)^-1);
    J = (p3^2 - 1 - r) \ d;
    js = shuffle(vector(J + 1, i, i - 1));
    for (i = 1, #js,
      m = r + js[i]*d; if (m < 15 || m >= p3^2, next);
      f = factor(m);
      if (#f~ != 2 || f[1, 2] != 1 || f[2, 2] != 1, next);
      p1 = f[1, 1]; p2 = f[2, 1];
      if (p2 >= p3 || !isinert(p1) || !isinert(p2) || !isfree(p1) || !isfree(p2), next);
      if ((p3^2 - 1) / Mp(p3, p1*p2*p3*S) > 1000,
        if (LCMFORM, next, error("class a construction failed at ", [p1, p2, p3, S])));
      return([p1, p2, p3])));
  0;
}

\\ class b: p3 large; d2 | p2^2-1 with p2/100 <= d2 <= p2/3; p1 == (p3 S)^-1 mod d2, p1 < p2 inert;
\\          then require M'(p3) small (idx3 >= idxmin) and p2/lcm(M'(p2),M'(p3)) <= 100
build_b(S) =
{
  my(p3, p2, N2, ds, d2, r, J, js, p1, n, M2, M3);
  for (t = 1, 20000,
    p3 = inert_in(X0 \ 2, X0, 0); if (!p3, next);
    p2 = inert_in(max(100, X0 \ 20), p3 - 1, random(2)); if (!p2, next);
    N2 = p2^2 - 1;
    ds = select(x -> 100*x >= p2 && 3*x <= p2 && gcd(x, p3*S) == 1, divisors(N2)); if (#ds == 0, next);
    d2 = ds[random(#ds) + 1];
    r = lift(Mod(p3*S, d2)^-1);
    J = (p2 - 1 - r) \ d2; if (J < 0, next);
    js = shuffle(vector(J + 1, i, i - 1));
    for (i = 1, #js,
      p1 = r + js[i]*d2; if (p1 < 3 || p1 >= p2, next);
      if (!isinert(p1) || !isfree(p1), next);
      n = p1*p2*p3*S; M2 = Mp(p2, n); M3 = Mp(p3, n);
      if ((p3^2 - 1)/M3 < idxmin(p3), break);         \\ M'(p3) not small enough: new p3
      if (p2 > 100*lcm(M2, M3), if (LCMFORM, next, error("class b construction failed")));
      return([p1, p2, p3])));
  0;
}

\\ class c: random triple with p3 in [X/2, X], p2 in [X/4, p3), p1 log-uniform; both M' small
build_c(S, extreme) =
{
  my(p3, p2, p1, n, M2, M3, L23);
  for (t = 1, 10^6,
    p3 = inert_in(X0 \ 2, X0, 0); if (!p3, next);
    p2 = inert_in(max(100, X0 \ 4), p3 - 1, 0); if (!p2, next);
    p1 = inert_in(3, p2 - 1, 1); if (!p1, next);
    n = p1*p2*p3*S; M2 = Mp(p2, n); M3 = Mp(p3, n); L23 = lcm(M2, M3);
    if ((p3^2 - 1)/M3 < idxmin(p3), next);
    if (extreme, if (L23 > 12, next), if (p2 < min(1000, p2 \ 48)*L23, next));
    return([p1, p2, p3]));
  0;
}

\\ nontrivial divisor of x (x = n-1 or n+1) built from its prime factors below 2^20, size about T (it can exceed
\\ T when the smallest odd prime factor does), with odd part > 1; 0 if impossible
pick_div(x, T) =
{
  my(f = factor(x, 2^20), pp = List(), i0 = 0, A, rest);
  for (i = 1, #f~, if (f[i, 1] < 2^20, for (e = 1, f[i, 2], listput(pp, f[i, 1]))));
  pp = Vec(pp);
  \\ start from the smallest odd prime factor (odd part > 1), then add random further factors while <= T
  for (i = 1, #pp, if (pp[i] > 2, i0 = i; break));
  if (!i0, return(0));
  A = pp[i0];
  rest = shuffle(vector(#pp - 1, j, pp[if (j < i0, j, j + 1)]));
  for (i = 1, #rest, if (A*rest[i] <= T, A *= rest[i]));
  if (A < 3 || x % A || A >= x, return(0));
  A;
}

\\ ---------- verification of one plant (used before writing, and by checkfile) ----------
\\ returns 0 if OK, else an error string. cls in {"a","b","c","c!","d/a","d/b","d/c","d/c!"} (! = extreme)
verify(p1, p2, p3, S, A, B, M1, M2, M3, cls) =
{
  my(n = p1*p2*p3*S, base, fS, L23);
  if (!(p1 < p2 && p2 < p3 && p3 <= X0 && p1 >= 3), return("order/range"));
  if (!isinert(p1) || !isinert(p2) || !isinert(p3), return("not inert prime"));
  if (S < 1, return("S < 1"));
  if (S == 1,
    if (A != 1 || B != 1, return("S = 1 needs A = B = 1"))
  ,
    fS = factor(S);
    for (i = 1, #fS~, if (fS[i, 2] != 1 || !issplit(fS[i, 1]), return("S not a squarefree product of split primes")));
    if (#fS~ > 2, return("S has more than two primes"));
    if (A < 3 || B < 3 || A >= n - 1 || B >= n + 1, return("A or B trivial"));
    if ((n - 1) % A || (n + 1) % B, return("W2 fails"));
    if (2 % gcd(A, B), return("gcd(A,B) does not divide 2")));
  if ((p1^2 - 1) % M1 || (p2^2 - 1) % M2 || (p3^2 - 1) % M3, return("M' does not divide p^2-1"));
  if ((n/p1 - 1) % M1 || (n/p2 - 1) % M2 || (n/p3 - 1) % M3, return("W1 fails"));
  if (M1 != Mp(p1, n) || M2 != Mp(p2, n) || M3 != Mp(p3, n), return("M' does not follow the override rule"));
  base = if (#cls >= 3, strsplit(cls, "/")[2], cls);
  if (#cls >= 3 && S == 1, return("class d needs S > 1"));
  if (#cls < 3 && S != 1, return("classes a,b,c need S = 1"));
  L23 = lcm(M2, M3);
  if (base == "a", if ((p3^2 - 1)/M3 > 1000, return("class a: idx3 > 1000")));
  if (base == "b",
    if ((p3^2 - 1)/M3 < idxmin(p3), return("class b: idx3 too small"));
    if (p2 > 100*L23, return("class b: p2/lcm > 100")));
  if (base == "c" || base == "c!",
    if ((p3^2 - 1)/M3 < idxmin(p3), return("class c: idx3 too small"));
    if (p2 < min(1000, p2 \ 48)*L23, return("class c: loop too short"));
    if (base == "c!" && L23 > 12, return("class c!: lcm > 12")));
  if (base != "a" && base != "b" && base != "c" && base != "c!", return("unknown class"));
  0;
}

plantline(k, cls, p1, p2, p3, S, A, B, M1, M2, M3) =
{
  my(n = p1*p2*p3*S, L23 = lcm(M2, M3), Lin = lcm([M2, M3, A, B]), nl = []);
  if (!islcmform(p1, M1), nl = concat(nl, p1));
  if (!islcmform(p2, M2), nl = concat(nl, p2));
  if (!islcmform(p3, M3), nl = concat(nl, p3));
  Str("# plant ", k, " class=", cls, " Sfac=", if (S == 1, "1", strjoin(apply(x -> Str(x), factor(S)[, 1]~), "*")),
      " idx3=", (p3^2 - 1)/M3, " p2/lcm(M2,M3)=", p2 \ L23,
      " inner=p2/lcm(M2,M3,A,B)=", p2 \ Lin, " candA=p3^2/lcm(M3,A,B)=", p3^2 \ lcm([M3, A, B]),
      " M'=", [M1, M2, M3], " nbits=", #binary(n), " nonlcm=", if (#nl, strjoin(apply(x -> Str(x), nl), ","), "none"));
}

\\ ---------- generator ----------
mkplants(D, X, N, seed, fname) =
{
  my(P = List(), tr, S, A, B, n, M, cls, e, sub, nq, q1, q2, f, cnt = Map(), k = 0, t0 = getwalltime(), lo, hi);
  D0 = D; Q0 = (1 - D)/4; X0 = X; USED = Map(); setrand(seed);
  if (D % 4 != 1, error("D must be 1 mod 4"));
  \\ classes a, b, c with S = 1
  foreach (["a", "b", "c"], c0,
    for (i = 1, N,
      tr = if (c0 == "a", build_a(1), if (c0 == "b", build_b(1), build_c(1, i % 3 == 0)));
      if (tr == 0, warning("could not build class ", c0, " plant ", i); next);
      cls = if (c0 == "c" && i % 3 == 0, "c!", c0);
      listput(P, [tr[1], tr[2], tr[3], 1, 1, 1, cls]);
      for (j = 1, 3, mapput(USED, tr[j], 1))));
  \\ class d: split part, rotating over a, b, c
  for (i = 1, N,
    sub = ["a", "b", "c"][(i - 1) % 3 + 1];
    for (t = 1, 1000,
      \\ S: one or two split primes; every 5th plant gets a huge S (two primes in [1e12, 1e14]: n > 2^128 typically)
      if (i % 5 == 0, lo = 10^12; hi = 10^14; nq = 2, lo = 3; hi = 10^5; nq = 1 + random(2));
      q1 = split_in(lo, hi); q2 = if (nq == 2, split_in(lo, hi), 1);
      if (q1 == q2, next);
      S = q1*q2;
      tr = if (sub == "a", build_a(S), if (sub == "b", build_b(S), build_c(S, i % 6 == 3)));
      if (tr == 0, next);
      n = tr[1]*tr[2]*tr[3]*S;
      \\ A, B: the c plants get small W2 moduli (long loops survive), a/b plants larger ones
      if (sub == "c", A = pick_div(n - 1, rndlog(3, 30)); B = pick_div(n + 1, rndlog(3, 30)),
                      A = pick_div(n - 1, rndlog(100, 10^6)); B = pick_div(n + 1, rndlog(100, 10^6)));
      if (A == 0 || B == 0, next);
      cls = Str("d/", if (sub == "c" && i % 6 == 3, "c!", sub));
      listput(P, [tr[1], tr[2], tr[3], S, A, B, cls]);
      for (j = 1, 3, mapput(USED, tr[j], 1));
      break));
  \\ verify every plant, then write
  f = fileopen(fname, "w");
  filewrite(f, Str("# mkplants.gp D=", D, " P=1 Q=", Q0, " X=", X, " seed=", seed, " plants=", #P, " mrule=", if (LCMFORM, "lcmform", "gcd"),
                   " (per class target ", N, "; classes a,b,c have S=A=B=1; d = split part)"));
  filewrite(f, Str("# override rule: M'(p) = ", if (LCMFORM, "lcmform(p, gcd(p^2-1, n/p-1)) (2-part <= max(v2(p-1),v2(p+1)))",
                   "gcd(p^2-1, n/p-1)"), ", n = p1 p2 p3 S; all plants pairwise disjoint"));
  for (i = 1, #P,
    my(v = P[i], M1, M2, M3);
    n = v[1]*v[2]*v[3]*v[4];
    M1 = Mp(v[1], n); M2 = Mp(v[2], n); M3 = Mp(v[3], n);
    e = verify(v[1], v[2], v[3], v[4], v[5], v[6], M1, M2, M3, v[7]);
    if (e, error("plant ", v, " fails: ", e));
    k++;
    filewrite(f, plantline(k, v[7], v[1], v[2], v[3], v[4], v[5], v[6], M1, M2, M3));
    filewrite(f, Str("O ", v[1], " ", M1));
    filewrite(f, Str("O ", v[2], " ", M2));
    filewrite(f, Str("O ", v[3], " ", M3));
    filewrite(f, Str("T ", v[1], " ", v[2], " ", v[3], " ", v[4], " ", v[5], " ", v[6]));
    mapput(cnt, v[7], if (mapisdefined(cnt, v[7]), mapget(cnt, v[7]), 0) + 1));
  fileclose(f);
  print("mkplants D=", D, " X=", X, " wrote ", k, " plants to ", fname, " counts=", Mat(cnt),
        " time=", strtime(getwalltime() - t0));
  k;
}

\\ ---------- independent re-check of a written file ----------
checkfile(fname, D, X) =
{
  my(L = readstr(fname), O = Map(), seen = Map(), cls = "", tok, w, nT = 0, nO = 0, bad = 0, cnt = Map(),
     p, M, e, v, nonlcm = 0, big = 0, passes = Map(), maxinner = 0);
  D0 = D; Q0 = (1 - D)/4; X0 = X; USED = Map(); LCMFORM = 0;
  \\ first pass: header (override rule) and overrides
  for (i = 1, #L,
    w = strsplit(L[i], " ");
    if (#w && w[1] == "#", for (j = 2, #w, if (w[j] == "mrule=lcmform", LCMFORM = 1)));
    if (#w && w[1] == "O",
      p = eval(w[2]); M = eval(w[3]); nO++;
      if (mapisdefined(O, p), bad++; print("duplicate override for ", p));
      if ((p^2 - 1) % M, bad++; print("override M does not divide p^2-1 at ", p));
      if (!islcmform(p, M), nonlcm++);
      mapput(O, p, M)));
  \\ second pass: plants
  for (i = 1, #L,
    if (L[i] == "", next);
    w = strsplit(L[i], " ");
    if (w[1] == "#",
      for (j = 2, #w, tok = strsplit(w[j], "="); if (#tok == 2 && tok[1] == "class", cls = tok[2]));
      next);
    if (w[1] == "T",
      v = vector(6, j, eval(w[j + 1])); nT++;
      my(lbad = 0);
      for (j = 1, 3,
        if (mapisdefined(seen, v[j]), lbad++; print("prime ", v[j], " reused"));
        mapput(seen, v[j], 1);
        if (!mapisdefined(O, v[j]), lbad++; print("no override for ", v[j])));
      if (lbad, bad += lbad; cls = ""; next);
      e = verify(v[1], v[2], v[3], v[4], v[5], v[6], mapget(O, v[1]), mapget(O, v[2]), mapget(O, v[3]), cls);
      if (e, bad++; print("plant ", v, " class ", cls, ": ", e));
      if (v[1]*v[2]*v[3]*v[4] >= 2^128, big++);
      mapput(passes, [v[4], v[5], v[6]], 1);
      maxinner = max(maxinner, v[2] \ lcm([mapget(O, v[2]), mapget(O, v[3]), v[5], v[6]]));
      mapput(cnt, cls, if (mapisdefined(cnt, cls), mapget(cnt, cls), 0) + 1);
      cls = ""
    ,
    if (w[1] != "O", bad++; print("unparsed line: ", L[i]))));
  if (#seen != nO, bad++; print("override count ", nO, " != planted primes ", #seen));
  print("checkfile ", fname, ": D=", D, " X=", X, " mrule=", if (LCMFORM, "lcmform", "gcd"), " plants=", nT, " overrides=", nO, " errors=", bad,
        " per-class=", Mat(cnt), " distinct(S,A,B)=", #passes, " n>=2^128: ", big,
        " non-lcm-form overrides=", nonlcm, " max inner loop=", maxinner);
  bad;
}

{
  if (getenv("MKP_CHECK"),
    checkfile(getenv("MKP_CHECK"), eval(getenv("MKP_D")), eval(getenv("MKP_X"))); quit);
  if (getenv("MKP_OUT"),
    LCMFORM = if (getenv("MKP_LCMFORM"), eval(getenv("MKP_LCMFORM")), 0);
    mkplants(eval(getenv("MKP_D")), eval(getenv("MKP_X")), eval(getenv("MKP_N")), eval(getenv("MKP_SEED")),
             getenv("MKP_OUT"));
    checkfile(getenv("MKP_OUT"), eval(getenv("MKP_D")), eval(getenv("MKP_X")));
    quit);
}
