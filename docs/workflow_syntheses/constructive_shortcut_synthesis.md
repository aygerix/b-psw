### Session 9 constructive-shortcut workflow: synthesis on one metric, mechanism ranking, next experiments [SYNTHESIS; audit status per claim; arithmetic re-checked here]

**Audit scope.** Only the five small-exponent algorithm claims were attacked, and all five survived (refuted = false). Only those are *established* below. Every other claim is tagged NOT ATTACKED: it can rank a direction but cannot close one.

What I checked myself in this pass:
- Theorem U's congruence, the sign lemma and Chernick corollary (a), by hand.
- The k-tree values (55.6 / 68.3 / 37.3), the 445 counting threshold, the bound |x|_1 >= 45, the ratios 40.9x and 97.8x, the deficit 6.51 and the rate 345.2 bits per doubling, all recomputed.
- scipy 1.18.1 with `milp` (HiGHS) is present on dirichlet.

#### Common metric
- G_C(b) = log2( N_req(b; 2^C) / F(b) ): how many doublings of N are missing at group size b bits with a time budget of 2^C.
  - F(b) is the best known number of usable primes with group size at most b.
  - N_req = b + 1 when 0.291 b <= C. This is the existence-limited regime, solved by Becker–Coron–Joux (BCJ) density-1 subset sum.
  - Otherwise N_req = min over ell <= C of ell * 2^(ceil((b+2)/ell) - 1) (k-tree).
- Mechanisms that add primes act on F. Algorithms act on N_req. Below, C = 80 throughout.

The G curve on the best existing frontier (fr_1500_17_os15.out = T1500 to 1e17 plus stragglers; the b = 213.7 point is from part2_200_b100.out):

| b | 139 | 178 | **213.7** | 325 | 686 | 978 | 1279 | 1585 | 1912 | 2201 | 2382 | 2525 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| N | 10 | 13 | 33 | 54 | 141 | 320 | 697 | 1366 | 2512 | 3995 | 4632 | 4451 |
| G_80 | 3.81 | 3.78 | **2.70** | 4.29 | 7.13 | 9.93 | 12.8 | 14.9 | 18.0 | 21.3 | 23.1 | 25.2 |

**Reading (COMPUTED arithmetic on COMPUTED frontiers).**
- The smallest gap is 2.70 doublings (6.5x) at b ~ 214. It sits in the existence regime: BCJ costs 2^62 there, but the expected number of solutions is 2^(33-214).
- G rises steadily with b, to 23.1 doublings at the existence threshold b* ~ 2370 and 25.2 doublings at the full instance.
- Every algorithmic mechanism tested this session (layered small-exponent solving, the Sylow split, lattices) acts only at b >= b*, where the gap is about 9x larger in log terms. Below b*, no algorithm can help: solutions do not exist.
- Two independent measurements agree on the small-b density:
  - Session 3 P1: annealer on the CHM y = 200 twin-smooth list, all sizes, N/B <= 0.155 and flat in y and B.
  - Session 9 sel line: part_opt2 on the side-rule-free census2 records, p < 1e10, true orders including stragglers, 0.1537 at (213.7, 33).
  - So the 6.5x small-b deficit is robust (COMPUTED lower bounds on F; not certified optimal).
- **The structural shortcut in these units:** a mechanism with ΔG <= -2.7 at some b <= 275, or ΔG <= -23 near b ~ 2400. The best lever found by any line is ΔG = -0.9.

#### Mechanism entries (ranked by best achievable G_80; type a = new prime source, b = algorithm, c = identity)

**1. Partition- and class-optimised smooth pools, small b (type a).** [COMPUTED, lower bound on F; NOT ATTACKED; reproduces Session 3 P1]
- Numbers:
  - Class (2,0): 33 primes at 213.7 bits, G_80 = 2.70.
  - Class (1,2): 28 at 222.6 bits.
  - Mod-4 side rule: 20 at 224 bits.
  - Optimising the side partition and 2-adic class gains x1.7-1.9 (ΔG -0.8 to -0.9).
- Correction: part_200.out reports "min (b+1)/N = 6.03 at (325, 54)", but b = 325 is past BCJ reach at C = 80 (0.291 * 325 = 94.6). Its true G_80 is 4.29, not 2.59. An existence ratio alone misleads above b ~ 275.
- Provenance: the (213.7, 33) point exists only on the VM (~/b-psw/src/cs/sel/part2_200_b100.out). Copy it back.
- Missing: x6.5 more mutually consistent primes at about 214 bits. Equivalently, F(b) >= b + 1 for some b <= 275. The census and Theorem S already account for every Frobenian source of such primes.
- Files: src/cs/sel/part_opt2.py, src/cs/sel/census2.c, src/cs/sel/out/part_200.out, src/cs/sel/out/c2_1e10_200.out

**2. Designed row sets (Pierpont q-1 = 2^a 3^b; Hamming q-1 5-smooth) (type a).** [COMPUTED, complete for x <= 1e17 / 1e15; NOT ATTACKED]
- Numbers:
  - Pierpont: N = 56, largest consistent subset 7 against 44-140 bits, so G >= 2.7. 1841 independent sets, 0 hits, E[hits] ~ 2e-6.
  - Hamming: largest consistent subset 21 against 387 bits, G ~ 4.2. 125,252,264 independent sets, 0 hits.
  - Density per consistent set is 0.04-0.05, against 0.15 for intervals. This agrees with the Session 7 result that designed sets do not beat intervals.
- Status: CLOSED in range by exhaustion, not by the metric.
- Missing: nothing to tune. The automatic small exponent exists only for rows with q - 1 | E (PROVED via Dirichlet; NOT ATTACKED).
- Files: src/cs/smallexp/enum_pier.out, src/cs/smallexp/enum_ham.out, src/cs/smallexp/pg_pier.out, src/cs/smallexp/pg_ham5.out, src/cs/smallexp/sunit_stream.c

**3. Pomerance T = 1500 frontier, greedy mod-4 rule (baseline, type a).** [COMPUTED]
- Numbers:
  - G_80 = 3.78 at (178, 13), 23.1 at (2382, 4632), 25.2 at the instance (2525.4, 4451).
  - Growth along Xmax: rho_X = 149.6 bits / log2(4451/3296) = 345 bits per doubling of N (recomputed).
- New data point: the T = 1000 one-sided run finished (~/b-psw/src/cs/sel/os_{1,2}_1000_16.out, VM only).
  - Pool per decade: 138 at 1e10, then 129, 131, 127, 122, down to 99 at 1e15. The T = 1000 pool is past its per-decade peak (COMPUTED).
  - Stragglers: 75 and 51 observed against 66.0 and 60.7 predicted, a ratio of 1.03-1.09.
- Missing: rho < b (b <= 60) on the large-b side.

**4. Frobenian selectors and index stragglers (type a).** [Theorem S, Lemma S1 and the exchange-rate corollary: PROVED by the line, NOT ATTACKED; I hand-checked the q-Sylow step of S1. 1/R law COMPUTED]
- Numbers:
  - Exchange rate exactly 1 bit of density per bit of index removed, for any finite Galois M.
  - Stragglers: 405 observed vs 409.8 predicted (p <= 1e15); 82 vs 87.5 (p <= 1e10).
  - ΔG = -0.10 to -0.15 (N up 7-11%; 5018 -> 5425).
- Missing: a non-Frobenian source whose absorbed rough part r(p) is forced algebraically and grows with p. Theorem S cannot see such a source; only heuristic saturation data bound it.
- Files: src/cs/sel/sel_check.py, src/cs/sel/census.c, src/cs/sel/onesided.c, src/cs/sel/frontier_sel.py, src/cs/sel/out/sel_check_3e6.out, src/cs/sel/out/census_1e10.out, src/cs/sel/out/fr_1500_17_os15.out

**5. Small-exponent layered 0/1 zero-sum algorithms (type b).** [ESTABLISHED: attacked, survived]
- Established results:
  - F_2^r: r + 1 elements (optimal, Davenport constant).
  - F_3^r: (r+1)^2; Davenport constant 2r + 1, brute-forced for r <= 3.
  - (Z/2^k)^r: (r+1)^k suffice. This is tight for this algorithm (explicit adversarial family) and far above the Davenport bound 3r + 1 for (Z/4)^r.
  - QR-relation lemma: for any odd prime l, (r+1)^2 elements give a nonzero relation with square coefficients, so ternary SIS_5 is in P.
  - F_5^r: (r+1)^5, with constant about 0.4 when the full pipeline runs.
- Claim corrections from the audit:
  - "Exactly r+1 pulled" should read "at most r+1". The target version needs the target in the span and has no r+1 bound.
  - test_full.out was generated with seed 7, not the default seed 1.
  - The F_5 exponent 4.5-4.7 is a small-r artifact: lucky one-sided blocks end runs early. The mean ratio at r = 10 is 0.32. True cost is Θ((r+1)^5).
  - "Deterministic polynomial time" in the QR lemma means poly(r, l) in zsolve. poly(r, log l) needs the deterministic-nu repair, implemented in audit_qr/qr_indep.py.
- NOT ATTACKED: the Fermat-prime extension, which is an argument only and implemented only for l <= 5. l = 7 is OPEN, and the obstruction to norm-type constructions is stated.
- Applied to BPSW: ΔG = 0. The best hybrid kills no layer. Killing the rank-238 2-layer costs 7.9 bits of log2 N and saves 238/ell levels: a loss for ell >= 31.
  - Only 22% of the 2526 bits sit in the shallow {2,3,5} layers. The generated group is full there (GF(2) rank 236/237, GF(3) rank 115/115).
  - Families closed under ±1 (exponent 2) are empty: 0 of 203,449 pairs in P_100.
- Missing: a pool whose generated group is shallow, so that N >= prod over layers of f_l(rho) is affordable. For example, a depth-1 3-group of rank 100 needs only about 6800 primes. Alternatively, an algorithm for l = 7 and deeper layers, which runs against SIS hardness for q >= beta n^delta (RECALLED).
- Files: src/cs/smallexp/zsolve.py, src/cs/smallexp/test_zsolve.py, src/cs/smallexp/scale5.py, src/cs/smallexp/layer_profile.py, src/cs/smallexp/prof_T1500.out, src/cs/smallexp/pm1_cliques.py, and the audit directories src/cs/smallexp/audit_e2/, src/cs/audit_l3/, src/cs/smallexp/audit_e4/, src/cs/smallexp/audit_qr/, src/cs/smallexp/audit_x5/

**6. Sylow split + GF(2) on G[2] (type b).** [COMPUTED cost law; NOT ATTACKED]
- Cost law: N_split = (r+4) * N_gen(B - r). The crossover r = b * log2(r+4) was confirmed at b = 8 and b = 10. Random columns give identical numbers.
- ΔG: about -4 at ell = 20 (irrelevant there, since G ~ 114), 0 at ell = 30, positive for ell >= 40.
- At N = 4451 the split cannot run: every block falls below density 1.
- Oracle bound: solving every l <= 101 for free still leaves 611 bits, against at most 418 bits of k-tree capacity.
- Supersedes the notebook's two-torus T7 reading that "splitting off the structured part is free".
- Missing: nothing tunable.
- Files: src/cs/rs/sylow.py, src/cs/rs/split_vs_gen.py, src/cs/rs/extrap.py, src/cs/rs/sylow_T1500.out, src/cs/rs/expA_b8_real.out, src/cs/rs/expA_b10_real.out

**7. Lattice / height structure (type b).** [COMPUTED; height bound PROVED; NOT ATTACKED]
- Real columns are indistinguishable from the uniform, support-preserving and height-preserving controls. BKZ-20 succeeds 39/40 at d ~ 34 and 0/160 at d >= 64, for all four column types.
- Height forces |x|_1 >= 45. Generic counting already forces about 445 (both recomputed). ΔG = 0.
- Missing: any column property that the controls lack. None was found here, nor in W5a-c.
- Per-row diagnostic, the large-b form of the missing ingredient: median row costs 240 bits per doubling; only 8 rows are below 50; per-row support is 1.9% against the 13% needed at b = 50.
- Files: src/cs/rs/lat.py, src/cs/rs/rowrate.py, src/cs/rs/expB_scale_bkz.out, src/cs/rs/expB_enum3040.out, src/cs/rs/expB_low.out, src/cs/rs/expB_dense.out

**8. Algebraic identities: universal forms, mirrors, T_2 members, polynomial Korselt (type c).** [Theorem U, sign lemma, Chernick corollary: PROVED, hand-checked here. Sharing lemma and automation budget: PROVED, NOT ATTACKED. Censuses: COMPUTED]
- Per family:
  - At least 2 dead members.
  - The residual on dead members is 0.83-0.91 of log2 p each, and it is m-independent: z(p_d) | N_d.
  - Hence finitely many BPSW psp per D*, and the dead members must be prime factors of U_{2^k N_d}.
- As a pool: N/B_L = 0.043 / 0.026 / 0.018 on the Lucas half alone, with the Fermat half granted free. That puts B_L at about 4.7e4-6.8e4 bits, so G_80 > 500.
- Automation is capped at log2(p_max) + 1 bits.
- Bounded censuses (quadratic, cubic, mixed): about 2800 all-prime tuples tested end to end, 0 hits.
- Missing: a one-sided family with no dead member, or a full Williams family (R1, Theorems 10-12).
- Files: src/cs/ident/univ_lib.gp, src/cs/ident/runU.gp, src/cs/ident/good1.py, src/cs/ident/two_bad.py, src/cs/ident/poolNB.gp, src/cs/ident/kill_check.gp, src/cs/ident/out/

**9. Trivial generated group: "doubly primitive" primes (type a, extreme case).** [Criterion PROVED, elementary, here; census PENDING]
- Criterion: for D = 5, any odd set of at least 3 inert primes with identical (ord_p 2, z(p)) is a strong BPSW psp. Why each condition holds automatically:
  - (5/n) = -1, so Selfridge picks D = 5.
  - ord | p - 1 and z | p + 1 make the odd congruences automatic.
  - Identical v2 values give one 2-adic class, and v2(n-1) >= j1, v2(n+1) >= j2 follow from p ≡ 1 mod 2^j1 and p ≡ -1 mod 2^j2.
- Existing coverage: Attempt D (both orders <= 1500, all such p <= 469561) lies inside Theorem E's range, whose production run is still going.
- Heuristic: each (ord, z) cell holds about 1/(φ(de) ln(de)) primes, so a 3-prime cell is extremely unlikely.
- Missing: any (ord, z) cell holding 3 primes.

**Ranking (lower G is closer to a counterexample):** 1 (2.70) < 2 (>= 2.7, closed in range) < 3 (3.78) < 4 (ΔG -0.15) << 5, 6, 7 (ΔG 0 against 23-25) << 8 (> 500). Entry 9 is unranked (pending).

#### Next 3 experiments (ranked)

**E1. Certify the small-b frontier exactly** (decides the only regime within 3 doublings).
- Pool: the CHM y = 200 list (data/sn_200.txt, about 99% complete, all sizes) together with the census2 records (p < 1e10, rows <= 200, true orders), plus the E2 primes.
- Model: an ILP with `scipy.optimize.milp` (HiGHS) per D (ten D) and per 2-adic class.
  - Variables: x_p, and y_{q,side,e} for each row q, side and exponent cap e.
  - Constraints: x_p <= y for each requirement of p; y_{q,F} + y_{q,L} <= 1 (a row lives on one side); group bits at most b, for b in {60, 80, ..., 300}.
  - Report both the LP-relaxation bound (the certificate) and the integer optimum.
- Prediction (HEURISTIC, from two independent lower bounds at 0.155 and 0.154): max N/(b+1) between 0.15 and 0.25 for every b <= 275 and every D, i.e. G_80 >= 2.
- Kill: certified N/(b+1) < 1/3 at every b <= 275 (and N < N_req for 275 < b <= 600). That closes smooth-type subset products in the BCJ regime for rows <= 200 (COMPUTED, conditional on list completeness).
- Confirm: any configuration with N >= b + 10 at b <= 275.
  1. Run BCJ or meet-in-the-middle (about 2^62).
  2. Check the result with src/bpsw.py.
  3. Exhaust that class exactly with Schroeppel–Shamir, not the independent-set enumerator that stalled on P_150. This would extend Theorem C.
- Cost: about 1 hour, 1 core, at most 2 GB.

**E2. Doubly-primitive census, D = 5, with d, e <= 1e4.**
- Method: a remainder tree of the primitive parts of Φ_d(2), reduced modulo each primitive Fibonacci part F*_e (about 3e7 bits, at most 200 MB). Factor the gcds, group the primes by (ord, z, j1, j2), and test every cell of size >= 2. Export all primes found as extra columns for E1.
- Prediction (HEURISTIC): the total count grows like about 1.4 ln X, so a few new primes appear above the Attempt D bound; no 2-prime cell; P(a 3-prime cell) < 1e-6.
- Kill: no cell of size >= 2. That closes the trivial-group mechanism for orders <= 1e4 (COMPUTED).
- Confirm: a 3-prime cell is a counterexample outright. Verify it with bpsw.py.
- Cost: at most 30 minutes. Run the other nine D only if a 2-prime cell appears.

**E3. Dead-member decision on the k = 3 shared-A Korselt variety** (the only theorem-level decision left in route c).
- Setup: f_j = 1 + A g_j with A = t^2 + pt + q and cofactors (c, t - r2, λ(t - r3)). Add the 2 Korselt equations, and for each "good" member j a rational a_j with f_j(a_j) = -1 and F(a_j) = -1. Saturate away the degenerate loci (r2 = r3, members coinciding, A with a repeated root).
- Run Singular under `timeout 1800`, memory at most 2 GB. If any stage comes out positive-dimensional, stop and add inequalities instead; this is the q22m_k7 lesson.
- Prediction (HEURISTIC, parameter count): one good member is possible (Lemma 5 has a non-dead Psi14 member); two good members give a finite or empty locus; three good members give the unit ideal.
- Kill: the unit ideal for all three. This PROVES that every k = 3 shared-A Korselt family has a dead member. Next, k = 5, via Theorem 10's topmost-'+'-point argument.
- Confirm: a point with all three members good is the first one-sided identity with no dead member. Scan it with bpsw.gp, compute its Lucas defect, and hand it to R6.

**Deferred (4th): the two-torus per-row yield screen (rs line's next).** Theorem S and the T3 Kummer-independence data already predict failure on one torus for every Frobenian or AGP-type candidate. It also serves the large-b regime, which is 23 doublings away against 2.7 at small b.

**Housekeeping:**
- Copy part2_200_b100.out and os_{1,2}_1000_16.{out,txt} from the VM into the repo.
- Record seed 7 for test_full.out.
- The exponent-5 audit reruns at r = 14 and 16 did not return before the audit was submitted.
- Theorem E production (R5) is still running.