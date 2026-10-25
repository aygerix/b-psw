\\ brute force with a split part S (validation of the engine's S/W2 code paths in the one-torus control modes)
\\ usage: SPL_D=5 SPL_X=6000 SPL_S=209 SPL_OUT=dir gp -q ref_split.gp
\\ For every inert triple p1 < p2 < p3 <= X (kronecker(D,p) = -1, p not dividing 2QD), n = p1 p2 p3 S:
\\  ford : ord(p) | n/p - 1 at the three p and n == 1 mod A   (A = lcm ord(q), q | S)   [= engine W1ONE fermat, --split S A 1]
\\  fz   : z(p)   | n/p - 1 at the three p and n == -1 mod B  (B = lcm z(q), q | S)     [= engine W1ONE lucas,  --split S 1 B]
\\  fermat: n strong 2-prp (tested on triples with ord(p) | n-1 for p | n and A | n-1)
\\  lucas : n strong Lucas, fixed D (tested on triples with z(p) | n+1 and B | n+1)
read("../ref_brute.gp");
run_split(D, X, S, outdir) =
{
  my(Q = (1 - D)/4, P = inert_primes(D, X), N, O, Z, fs = factor(S)[, 1], A = 1, B = 1, n, m, Lo = List(), Lz = List(), Lf = List(), Ll = List());
  for (i = 1, #fs, my(q = fs[i]); if (kronecker(D, q) != 1 || (2*Q*D) % q == 0, error("not split: ", q));
    A = lcm(A, ord2_naive(q)); B = lcm(B, zrank_naive(q, Q)));
  N = #P; O = vector(N, i, ord2_naive(P[i])); Z = vector(N, i, zrank_naive(P[i], Q));
  for (i = 1, N - 2, for (j = i + 1, N - 1, m = P[i]*P[j];
    for (k = j + 1, N, n = m*P[k]*S;
      if ((n/P[i] - 1) % O[i] == 0 && (n/P[j] - 1) % O[j] == 0 && (n/P[k] - 1) % O[k] == 0 && (n - 1) % A == 0, listput(Lo, [n, P[i], P[j], P[k]]));
      if ((n/P[i] - 1) % Z[i] == 0 && (n/P[j] - 1) % Z[j] == 0 && (n/P[k] - 1) % Z[k] == 0 && (n + 1) % B == 0, listput(Lz, [n, P[i], P[j], P[k]]));
      if ((n - 1) % A == 0 && (n - 1) % O[i] == 0 && (n - 1) % O[j] == 0 && (n - 1) % O[k] == 0 && sfermat2(n), listput(Lf, [n, P[i], P[j], P[k]]));
      if ((n + 1) % B == 0 && (n + 1) % Z[i] == 0 && (n + 1) % Z[j] == 0 && (n + 1) % Z[k] == 0 && slucasD(n, D), listput(Ll, [n, P[i], P[j], P[k]])))));
  my(h = Str("# ref_split.gp D=", D, " X=", X, " S=", S, " A=", A, " B=", B, " inert=", N));
  writelist(Str(outdir, "/D", D, "_X", X, "_S", S, "_ford.txt"), Str(h, " pred=ford"), Lo);
  writelist(Str(outdir, "/D", D, "_X", X, "_S", S, "_fz.txt"), Str(h, " pred=fz"), Lz);
  writelist(Str(outdir, "/D", D, "_X", X, "_S", S, "_fermat.txt"), Str(h, " pred=fermat"), Lf);
  writelist(Str(outdir, "/D", D, "_X", X, "_S", S, "_lucas.txt"), Str(h, " pred=lucas"), Ll);
  print("SPLIT D=", D, " X=", X, " S=", S, " A=", A, " B=", B, " ford=", #Lo, " fz=", #Lz, " fermat=", #Lf, " lucas=", #Ll);
}
{ run_split(eval(getenv("SPL_D")), eval(getenv("SPL_X")), eval(getenv("SPL_S")), getenv("SPL_OUT")); quit; }
