# Theorem E engine: specification (two-torus anchor enumeration)

Goal. For a fixed Selfridge discriminant D (P = 1, Q = (1 - D)/4) and a bound X, decide whether there is a
squarefree n passing the exact strong BPSW test (src/bpsw.py) such that

  n = p1 p2 p3 * S,   p1 < p2 < p3 <= X primes, all three INERT for D ((D/p) = -1), p_i not dividing 2QD,
  S = 1 or a fixed "split part" (product of admissible split primes, supplied by the caller).

Every n found by the engine that passes exact BPSW is a counterexample and must be printed with its factors.
A run with no hits proves: no BPSW psp of this shape (for this D, this S, primes <= X).

## Necessary conditions used for enumeration (all PROVED; register rows F1, F2)

For a prime p not dividing 2QD put ord(p) = multiplicative order of 2 mod p, z(p) = rank of apparition of
U_k(1, Q) mod p (least k > 0 with p | U_k). For inert p: ord(p) | p - 1, z(p) | p + 1.
Put M(p) = lcm(ord(p), z(p)) (note gcd(ord, z) | 2) and I(p) = (p^2 - 1)/M(p) (the "index", an integer).

Strong Fermat base 2 implies 2^(n-1) == 1 mod n, hence ord(p) | n - 1 for p | n.
Strong Lucas (with (D/n) = -1) implies U_(n+1) == 0 mod n, hence z(p) | n + 1 for p | n.

(W1) For inert p | n with cofactor c = n/p:  M(p) | c - 1.
     Proof: p == 1 mod ord(p) gives n - 1 == c - 1 mod ord(p); p == -1 mod z(p) gives n + 1 == 1 - c mod z(p).
(W2) For the split part S: n == 1 mod A_S and n == -1 mod B_S, where A_S = lcm of ord(q), B_S = lcm of z(q) over
     the primes q | S. (Caller supplies S, A_S, B_S; S = 1 means A_S = B_S = 1.)

These are the ONLY conditions the enumeration may rely on (plus p1 < p2 < p3 <= X and the inert/primality of the
p_i). Do NOT add 2-adic class filters or other strong-test refinements to the enumeration: completeness of the
search must follow from (W1), (W2) alone. The final test may be anything (exact BPSW in production).

## Enumeration (per D, per split part S)

Let P3 range over the inert primes p3 <= X. Write m = p1 p2 (so m < p3^2) and n = p3 m S.

Branch A (anchor at p3; used when the candidate count is small).
  (W1) at p3: m S == 1 mod M(p3).  (W2): p3 m S == 1 mod A_S, p3 m S == -1 mod B_S.
  CRT-combine into m == r mod L (L = lcm(M(p3), A_S, B_S)); if incompatible or a needed inverse does not exist,
  there is nothing to do for this p3. Enumerate m = r, r + L, ... < p3^2 (m >= 15). Each m gives one candidate
  n = p3 m S. Candidate count ~ p3^2 / L <= I(p3).
  In bpsw mode a cheap prefilter is allowed: Fermat base 2 modulo m (2^(n-1) == 1 mod m), which is necessary.
  Survivors get the exact BPSW test on n. Do NOT require m to factor as p1 p2 in bpsw mode (that only enlarges
  the tested set); in the control modes, see below.

Branch C (tail: p3 whose branch-A count is too large, i.e. M(p3) small).
  Loop over every inert prime p2 < p3 (from the table). Then p1 is fixed by
    (W1) at p3: p1 == (p2 S)^(-1) mod M(p3);
    (W1) at p2: p1 == (p3 S)^(-1) mod M(p2);
    (W2):       p1 == (p2 p3 S)^(-1) mod A_S,  p1 == -(p2 p3 S)^(-1) mod B_S.
  CRT-combine into p1 == r' mod L'; enumerate p1 = r', r' + L', ... < p2 (p1 >= 3), keep inert primes p1 that
  satisfy (W1) at p1 (M(p1) | p2 p3 S - 1), then the final test on n = p1 p2 p3 S.
  For pairs with p2 / L' large the inner loop is long; that is allowed but must be COUNTED and logged
  (log every pair with more than 1e6 inner iterations).

Branch choice. Per p3, estimate cost_A = (number of m) * c_test and cost_C = (number of inert p2 < p3) * c_pair
(+ inner iterations). Use the cheaper branch. A command-line "--cap K" overrides: K = 0 forces branch C for
every p3, K = infinity forces branch A for every p3. Results must not depend on K (this is a required control).

Completeness argument (must be preserved by the implementation): every n of the stated shape satisfies (W1) at
p1, p2, p3 and (W2). If p3's branch is A, its m = p1 p2 lies in the enumerated progression below p3^2. If it is
C, its p2 is visited and its p1 lies in the enumerated progression below p2 and passes the (W1)-at-p1 check.

## Arithmetic and data

 - Table: for every inert prime p <= X (p not dividing 2QD), ord(p) and z(p) EXACTLY (factor p - 1 and p + 1 by a
   segmented sieve, then strip prime factors while 2^(e/l) == 1 resp. U_(e/l) == 0 mod p). See src/us/usplit.c
   for a 1-sided version of this sieve and for lucasU / ord2 / zrank.
 - Sizes: X up to 1e10 at least. Then p^2 ~ 1e20 > 2^64, M(p) can exceed 2^64 (it is < p^2/2), m < p3^2 exceeds
   2^64, n = p3 m S can exceed 2^128 when S is large. Use unsigned __int128 where it suffices and GMP (installed:
   /usr/include/x86_64-linux-gnu/gmp.h) or 2-limb Montgomery where it does not. Every fast-path arithmetic routine
   must be cross-checked against GMP on random inputs in a self-test mode ("--selftest").
 - Memory: at X = 1e10 there are ~2.3e8 inert primes; keep the table compact (e.g. p as u64 or delta-coded,
   indices (p-1)/ord and (p+1)/z as u32 with an overflow side table). Target < 6 GB at X = 1e10.
 - Parallelism: OpenMP over p3 (branch A) and over tail primes (branch C). 8 cores on the VM, other jobs running;
   default to --threads 4.
 - Exact BPSW: implement with GMP exactly as src/bpsw.py (strong Fermat base 2; Selfridge Method A: D = 5, -7, 9,
   -11, ..., first with Jacobi(D, n) = -1, Jacobi 0 skipped; perfect squares rejected; strong Lucas P = 1,
   Q = (1 - D)/4). Self-test: agree with a reference on all odd n < 1e6, on the known strong Lucas pseudoprimes
   (5459, 5777, 10877, 16109, 18971, ...), on the 2-spsps below 1e6, and on 1e5 random 100-bit odd n.

## Modes (for controls)

 --mode bpsw     production: final test = exact BPSW on n (prefilter allowed). Report every hit with factors.
 --mode weak     final test = (W1) at all three primes + (W2) using the table M (no BPSW). Requires factoring m in
                 branch A: only allowed for small X (<= 1e7); m must be p1 p2 with p1 < p2 < p3 inert primes.
 --mode plant    like weak, but M(p) is overridden for the primes listed in --plants FILE, and the final predicate
                 is "n is one of the planted triples" (hash lookup); report which plants were found / missed.
 --mode fermat   one-torus control: M(p) := ord(p) only, final test = strong Fermat base 2 on n AND m = p1 p2 with
                 p1 < p2 < p3 inert primes. (Comparable with a brute force over inert triples.)
 --mode lucas    one-torus control: M(p) := z(p) only, final test = strong Lucas with FIXED D (not Selfridge) on n
                 AND m = p1 p2 inert as above.

Plant file format (text): lines "O p M" (override M(p) := M; M must divide p^2 - 1) and lines
"T p1 p2 p3 S A B" (a planted solution; S A B = 1 1 1 for no split part). The engine runs one pass per distinct
(S, A, B) appearing in the plant file.

## Output

One summary line per (D, S) pass: X, number of inert primes, branch-A primes and candidates, branch-C primes,
pairs and inner iterations, prefilter survivors, BPSW hits, wall time. Hits: "HIT n = p1 * p2 * p3 * S".
