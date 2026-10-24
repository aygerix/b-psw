# Number Theory Research Notebook — Baillie-PSW pseudoprime

## Target (exact statement)
- Find odd composite n = f1*f2 (f1,f2>1) that is a strong probable prime to base 2 AND a strong Lucas
  probable prime with Selfridge Method A (P=1, D = first of 5,-7,9,-11,... with (D/n)=-1, Q=(1-D)/4).
- YES witness: explicit (f1, f2). NO theorem: out of scope (Pomerance's heuristic predicts infinitely many).
- Status in literature (checked 2026-10-17): none known; none < 2^64 (Gilchrist/Feitsma), none < 1e21
  for the "FL" version (Ishmukhametov et al., Izv. VUZ Mat. 2024 no.4). Listed as open by Epoch AI.

## Assumption register (added 2026-10-17 after the Williams-framework error; KEEP CURRENT)
Every lemma must list which rows it uses. A route may be called CLOSED only by a PROVED result whose hypotheses
are all FORCED; a result using a CHOSEN row closes only that sub-route, and the notebook must name the sub-route.
| id  | assumption                                   | status | why                                                        |
| F1  | one 2-adic class (j1,j2) for all p | n       | FORCED | Lemma 1 (strong tests act prime by prime)                  |
| F2  | gcd(odd ord / odd z moduli of n-1 and n+1)=1 | FORCED | n == 1 and n == -1 mod an odd prime is impossible          |
| F3  | n squarefree                                 | FORCED*| *unless p is Wieferich AND Wall-Sun-Sun-type for the D used |
| C1  | restrict to n with (5/n) = -1 (D = 5)        | CHOSEN | n with (5/n) = +1 get D = -7, -11, 13, ...; equally allowed |
| C2  | every prime factor inert for D               | CHOSEN | split primes allowed by Lemma 1 (heuristically rare)       |
| C3F | Fermat-side bounded index: p - 1 | n - 1     | CHOSEN | BPSW needs only odd(ord_p 2) | n - 1 (C3F = Korselt)      |
| C3L | Lucas-side bounded index: p + 1 | n + 1      | CHOSEN | BPSW needs only odd(z(p)) | n + 1 (C3F + C3L = Williams)  |
| C4  | n = prod f_j(t), divisibility as identity in Q[t] | CHOSEN | polynomial families; a_j > 0 is forced WITHIN C4 for prime values f_j(t) -> +inf |
| C5  | smooth p -+ 1 pools (Lemma 2)                | CHOSEN | used by Attempts A, B, D                                    |
| C6  | exponential parameter forcing an order (Y = L_2K, X = 2^(2^a)) | CHOSEN | Attempt C                              |

## Working rules (procedure changes after repeated failures; check before every log entry)
 G1 Assumption gate: before proving or computing anything, write which register rows it uses. If it uses a CHOSEN
    row, first write one line on whether BPSW forces that row. (The Williams error: Lemma 3 already said "bounded
    index" and nobody asked whether BPSW needs it.)
 G2 Obstruction budget: an obstruction gets ONE attempt entry. The next LOG attempt must be constructive (could output
    n), and its derivation must cite what the obstruction says a working object must look like (as Lemma 4 ->
    split torus -> Attempt C).
 G3 Resource cap: a computation that hits a memory or time cap gets at most ONE re-encoding. A second failure is
    recorded as UNFINISHED and cut. No escalation over Q, over GF(p), msolve and Calyx in turn.
 G4 No project-level stop language ("no mechanism left", "STOP", "dead" for the whole project). Every checkpoint lists
    >= 3 live constructive directions, each with first experiment and kill condition; if fewer exist, the next
    action is to derive new ones from the latest obstructions.
 G5 Heuristics rank, never close (see register rule above).
 G6 Source tags on every claim: [PROVED here] [COMPUTED: file] [READ by me: source] [REPORTED by X, unchecked].
    Before a checkpoint, grep it for untagged claims and either verify or tag them.
 G7 New construction => full property profile before moving on: all bases (Korselt?), Selfridge D actually chosen,
    2-adic classes, A_F, B_F, Lucas defect. (Lemma 5 sat for a whole session labelled "base-2 Fermat" when it is
    Carmichael.)

## Verified tools
- `src/bpsw.py`: verbatim implementation of the problem's tests. Reproduces OEIS A001262 (spsp(2) < 1e5)
  and A217255 (strong Lucas psp < 1e5) exactly; all primes < 2e5 pass; no BPSW psp < 2e6.

## Lemma 1 (exact per-prime criterion; PROVED, standard; numerically checked)
Register rows: F3, C1 (stated for D = 5; the same proof works for any Selfridge D with z = rank of apparition; D = -7 checked in Attempt B). It is the source of F1, F2. Closes nothing.
Let n be squarefree, odd, (5/n) = -1 (so Selfridge picks D=5, P=1, Q=-1, U=Fibonacci). Let
z(p) = Fibonacci entry point = order of gamma = alpha/beta mod p; s = v2(n-1), t = v2(n+1). Then
n is a BPSW psp  <=>  there are j1 <= s, j2 <= t with, for every p | n:
   v2(ord_p 2) = j1,  odd(ord_p 2) | n-1,   v2(z(p)) = j2,  odd(z(p)) | n+1.
Proof: U_k = 0 mod p <=> gamma^k = 1; V_k = 0 mod p <=> gamma^k = -1 (sqrt5 invertible, p != 5);
2^{d 2^r} = -1 mod p <=> v2(ord)=r+1 and odd(ord) | d; CRT over squarefree n.
Facts used (proved via alpha^{p+1} = alpha*beta = -1 for inert p):
  inert p = 3 mod 4: v2(z(p)) = v2(p+1);  inert p = 1 mod 4: z(p) odd;  p = 3,5 mod 8: v2(ord_p 2)=v2(p-1).
Check (`src/lemma_check.py`): 0 violations of the 2-adic facts for inert p < 2e5; criterion agrees with
direct test on 20,009 random squarefree n (D=5); half-criteria agree on all 129,580 squarefree n<1e6.
Non-squarefree n: p^2 | n forces z(p^2)=z(p), i.e. p a Wall-Sun-Sun prime (none known) -> n squarefree.

### Consequence: 2-adic signature classes (cannot be mixed inside one n)
inert p = 3 mod 8: (j1,j2)=(1,2) [class F2];  inert p with v2(p+1)=k>=3: (0,k) [F_k];
inert p = 1 mod 4: (v2(ord_p 2), 0) [E_j], e.g. all p = 5 mod 8 lie in E_2.
Split primes ((5/p)=+1) need gcd(odd ord_p2, odd z(p)) = 1 with both | p-1: heuristically rare.

## Lemma 2 (Pomerance-type sufficient condition; PROVED from Lemma 1)
Register rows: F1, F2, F3, C1 (or the D = -7 analogue), C2. Sufficient condition; closes nothing.
Fix D (=5 or -7 with the earlier Selfridge D's non-(-1)). Let S be a set of distinct primes, all inert for D,
all in ONE 2-adic class (j1,j2), with odd(ord_p 2) | L1 and odd(z(p)) | L2, gcd(L1,L2)=1, |S| odd, and
v2-room automatic (true for classes F2, F_k, E_j with j<=2). If prod S = 1 mod L1 and = -1 mod L2
(and the parity condition for earlier D's), then n = prod S is a BPSW psp.  Search space: subset products
in G = (Z/L1)^* x (Z/L2)^* x (Z/2)^c; heuristic #solutions ~ 2^{|P|}/|G| for a pool P.
Note (literature): Chen-Greene (Fib. Quart. 41 (2003) 334-344) did this for the WEAK $620 version
(2^n=2, F_{n+1}=0); |P| ~ log2|G| only at ~1240 bits (|P|=1242), 4838 primes at 1838 bits.

## Attempt A — strong version, D=5, exhaustive smooth-neighbour pools (bounded evidence)
Register rows: C1, C2, C5. Bounded evidence only.
- Pools: p = 2x+1 with x, x+1 z-smooth (CHM closure `src/chm.c`; z=100: 13,340 pairs vs 13,374 published,
  z=41: 868 vs Lehmer's 869 -> closure ~99.7% complete), 5 not | x(x+1), p prime; true orders used.
- z=100: 804 inert primes; classes (2,0):181, (1,2):140, (0,3):98.  z=150: classes (2,0):837, (1,2):556.
- Partition A/B/excluded optimised by annealing along the Pareto frontier (`src/frontier.py`):
    z=100 best (log2|G|, N): (59, 8), (126, 16)       -> expected solutions 2^-110
    z=150 best (log2|G|, N): (212, 26), (272, 30)     -> 2^-242; F2 class (230, 24)
  N/log2|G| ~ 0.11-0.13, no improvement from z=100 to 150.
- Search-cost side (HEURISTIC, standard): Wagner k-tree with 2^a disjoint leaf lists of b primes kills
  b(a+1) bits, needs N >= b*2^a primes, time ~2^(a+b). At 2^38 work: kills <= ~250 bits and needs
  N >= ~4000. Pools at <= 300 bits have N ~ 30.  At the scale where 2^N >= |G| (>= ~1200-2500 bits),
  cost >= 2^80.  => HEURISTIC: at the measured pool densities this route gives no example at feasible cost.
  Deprioritised, NOT closed (it uses C5; denser pools or a structured subset-product algorithm would reopen it).

## Attempt B — D = -7 (Q = 2) variant (bounded evidence)
Register rows: D = -7 in place of C1, C2, C5. Bounded evidence only.
- Requires (5/n)=+1, (-7/n)=-1; inert means (p/7)=-1; prime 5 now usable on either side.
- Criterion (Lemma 1 analogue) checked vs direct test: 7,475 agree / 0 disagree.
- Raw pools ~2x larger (z=100: 1527 inert primes; z=150: 6801), but optimised frontiers no better:
  z=100 (100.5 bits, N=10); z=150 (174.5, 17), (183.8, 14).  Ratio ~0.1.  Choice of D does not help.
- Side observation (PROVED, elementary): for Q=2 and inert p, alpha^{p+1}=2 so ord_p(2) and z(p) are
  both controlled by ord(alpha); but the count of inert p with small ord(alpha)=R is heuristically
  ~2^{omega(R)}/(sqrt(R) log R), i.e. again dominated by smooth p^2-1. No new source of primes.

## Lemma 3 (PROVED): no linear (Chernick-type) polynomial family
Register rows: C4, C3F, C3L, plus same-sign slopes. Closes ONLY linear same-sign Williams families.
There are no k>=2 DISTINCT linear f_j(t)=a_j t+b_j in Q[t] with all a_j of one sign such that
(f_j - 1) | (n - 1) and (f_j + 1) | (n + 1) in Q[t] for all j, where n = prod f_j.
Proof. WLOG a_j > 0. Let r_j be the root of f_j - 1 and r* = max r_j. Each f_j is increasing, so
f_j(r*) >= f_j(r_j) = 1, while n(r*) = 1; hence every f_j(r*) = 1, i.e. all r_j = r*, and
f_j = 1 + a_j u (u = t - r*) with distinct a_j. Then n + 1 vanishes at the k distinct points
u = -2/a_j, so n + 1 = C prod (u + 2 d_j), C = prod a_j, d_j = 1/a_j, while n = C prod (u + d_j).
Comparing u^{k-1}: sum d_j = 2 sum d_j, so sum d_j = 0, impossible.  QED
Relevance: integer conditions odd(f_j(t)-1) | n(t)-1, odd(f_j(t)+1) | n(t)+1 for infinitely many t
(bounded 2-adic/order index) force the Q[t]-divisibility above, so no Chernick-type family of
Williams/Pomerance-type BPSW psp exists. (Same for f_j = g_j(h(t)) with g_j linear.)
Higher-degree families: NOT settled here (argument needs real roots of f_j - 1).

## Heuristic eliminations (HEURISTIC, not proofs)
- Two prime factors: q-1 must absorb lcm(odd ord_p 2, odd z(p)) ~ p^2, then q needs both orders ~sqrt(q);
  expected count sum_p O(1/(p^2 log p)) -> essentially none.
- Primes with both ord_p(2), z(p) << p: expected O(log log) per (L1,L2) beyond smooth-p^2-1 primes.
- Mersenne / Wagstaff composites are automatically spsp(2), but Lucas side then needs every factor rho to
  have z(rho) | (power of 2)*(structured) -- probability ~ prod 1/rho.

## Strategic checkpoint (2026-10-17) — decision: STOP constructive search
- New information: exact criterion; first strong-version (single 2-adic class) pool measurements; D=-7
  variant no better; linear Chernick-type families impossible.
- Theory of success: none credible within budget. Pomerance-type examples should exist only at
  |G| >~ 2^1200 (weak version, Chen-Greene) and later for the strong single-class version; search >= 2^80.
- Caveats: frontier points are best-found by annealing (lower bounds on achievable N), CHM pair lists are
  ~99.7% complete (32 pairs > 2^63 dropped at z=150). Neither caveat can close a ~240-bit deficit.
- What would change the decision: (i) an efficient high-density subset-product algorithm for groups
  with many small odd-order cyclic factors; (ii) a mechanism producing many inert primes with BOTH
  ord_p(2) and z(p) small without p^2-1 being smooth; (iii) a higher-degree polynomial family.
- No candidate n was found; no submission file was written.

## Checkpoint superseded (2026-10-17, user directive): research resumes. New attempts below.

## Lemma 4 (PROVED, computer-assisted): no k=3 polynomial family of any degree
Register rows: C4, C3F, C3L. Closes ONLY k = 3 Williams families.
Setting: f1,f2,f3 in Q[t] nonconstant, n=f1f2f3, with f_j-1 | n-1 and f_j+1 | n+1 in Q[t]
(the Williams/Lucas-Carmichael identity that a Chernick-type family needs; 2-power and
constant denominators are units in Q[t], so this covers the strong-BPSW relaxation by constants).
1. Condition <=> P_j := n/f_j == 1 mod (f_j^2-1). deg(f_j^2-1)=2d_j <= deg P_j = D-d_j => 3d_j <= D,
   so d1=d2=d3=d and P_j - 1 = c_j (f_j^2-1) with constants c_j != 0. (Constant f_j impossible.)
2. So t -> (f1,f2,f3) is a polynomial curve on Q1∩Q2∩Q3, Q1: yz-1=c1(x^2-1), Q2: xz-1=c2(y^2-1),
   Q3: xy-1=c3(z^2-1). Singular (src/alg/k3_elim.sing): eliminating y,z and requiring K to vanish
   identically in x gives the c-locus  c1+c2+c3 = 3, c1c2c3 = 1  (a nodal cubic, rational).
3. Parametrised c=(-1/(m(m+1)), -m^2/(m+1), (m+1)^2/m); primary decomposition over Q(m)
   (src/alg/k3_comp.sing) gives ONE component: x = -m y + (m+1) z and
   m^2y^2 - m(m+1)yz + (m+1)^2 z^2 = m^2+m+1, i.e. a^2-ab+b^2 = N(m-w) with a=my, b=(m+1)z:
   a norm-one torus of Q(sqrt-3) -- an ellipse (definite form), no real point at infinity.
   The node c=(1,1,1) gives only the lines x=y=z (real) and x=wy=w^2z (complex).
4. A polynomial curve has a real point at infinity => impossible. Also positivity: c_j>0 with
   sum 3, product 1 forces c=(1,1,1) by AM-GM.  QED.
Extracted shape for a working family: k>=5 needed; the k=3 obstruction is that the solution
variety is an ANISOTROPIC torus. A split torus would give exponential (Lucas-sequence) families.

## Attempt C — exponential (cyclotomic) half-families (2026-10-17)
Register rows: C4, C6, and C3F (Lucas-exponential members) or C3L (Fermat-exponential members); C1 for the mod-5 argument.
Idea: control one order by an exponential parameter and the other side by a polynomial identity.
 - Fermat-exponential: p_j = Phi_{d_j}(X), X = 2^(2^a) (ord_p 2 | 2^a d_j), need Phi_j + 1 | prod_{i!=j}Phi_i - 1.
   Necessary: Res(Phi_a+1, Phi_b) = +-1 both ways. RESULT [COMPUTED: src/alg/cyc_sym.gp -> cyc_sym.out, rerun
   2026-10-17 with a,b <= 160]: 0 symmetric pairs. The pair condition is necessary for every k, so Fermat-exponential
   families are excluded when ALL d_j <= 160 (any k). d_j > 160: UNCHECKED (bounded computation, not a proof).
 - Lucas-exponential: p_j = Psi_{d_j}(Y), Psi_d = min poly of 2cos(2pi/d), Y = L_{2K}
   (z(p) | 2K d_j since Psi_d(L_{2K}) w^{phi/2} = Phi_d(w), w = alpha^{2K}).
   Need Fermat identity Psi_j - 1 | prod_{i!=j} Psi_i - 1.
   RESULT (src/alg/cyc_sets.gp, d<=100): exactly the triples {10,14,18}, {20,28,36}, {40,56,72}.
   LEMMA 5 (PROVED, direct check): for n = Psi10(Y)Psi14(Y)Psi18(Y),
     Register rows: C4, C3F (Fermat side only; no C3L). A CONSTRUCTION, closes nothing.
     Psi10 = Y^2-Y-1, Psi14 = Y^3-Y^2-2Y+1, Psi18 = Y^3-3Y-1, each Psi_j - 1 divides n - 1 in Z[Y]
     (monic quotients), because Psi10-1=(Y-2)(Y+1), Psi14-1=Y(Y-2)(Y+1), Psi18-1=(Y-2)(Y+1)^2 and all
     three Psi take values +-1 at Y in {2,-1,0} with product 1 (and derivative condition at -1).
   So for every integer Y >= 3 with all three values prime, n is a CARMICHAEL number (Korselt: p_j - 1 | n - 1 and n
   squarefree), i.e. a Fermat pseudoprime to EVERY base coprime to n [PROVED: identity in Z[Y]; COMPUTED:
   src/fam/lemma5_count.gp, 76 Y in [3,20000), Korselt holds for all]. Beeger-Duparc form: p = Psi10,
   h2 = Y^2-2, h3 = (Y-1)^2, Delta = h2 h3 - p^2 = 2Y-3, Psi14 = 1 + (p-1)(p+h3)/Delta, Psi18 = 1 + (p-1)(p+h2)/Delta
   [CHECKED by me: p+h3 = Y(2Y-3), p+h2 = (2Y-3)(Y+1)]. The STRONG base-2 test additionally needs one 2-adic class
   (F1). Negative Y (two values negative) gives |b| + 1 | n - 1 instead: not Carmichael. See the Lemma 5 profile.
   Screen Y = L_M, M <= 600 (src/fam/screen.gp): weak Fermat only for M <= 3 (values not all prime);
   no Lucas pass. Obstruction: Lucas control forces Y = L_{2K} (only units of the norm-1 torus),
   so values grow like phi^{2K}; P(all three prime) ~ K^-3 => finitely many expected. Also for
   K = 2^a, Y == 2 mod 3,5,7 forces n == 1 mod 315 while z(p_j) contains 5,7,9 => mod-315 clash.
 - Structural fact found: Psi_d(Y) +- 1 splits into cyclotomic Psi_e's for (almost) every d
   (prime d: Psi_d(2cos t) = sin(dt/2)/sin(t/2), and Phi_d(x) -+ x^m are cyclotomic quadrinomials).
   => a FULL Williams family (Psi_j - 1 | n-1 AND Psi_j + 1 | n+1) is a finite search over sets of d
   (cyclotomic-unit relations). Running: src/fam/will.gp (k <= 6, d <= 90).
 - RESULT will.gp: graph of pairs with unit resultants for both Psi_a -+ 1 vs Psi_b (d <= 90): 303 edges;
   cliques of size 3..6: 1337/4055/8638/13131; full Williams identity: NONE.

## Lemma 6 (PROVED): real special points of a Williams polynomial family are common
Register rows: C4, C3F, C3L, positive leading coefficients. A tool for Williams families only.
Setting: f_1..f_k in R[t] nonconstant, positive leading coefficients, n = C * prod f_i (C > 0 const),
Williams condition: for each j, P_j := n/f_j == 1 at every root of f_j^2 - 1 (with multiplicity).
Claim (C = 1): if t0 is real and f_j(t0) = +-1 for some j, then f_i(t0) = f_j(t0) =: E for ALL i,
with equal multiplicity; E = -1 only if k odd.
Proof: induction downward from +inf over the finitely many real special points. On (R,inf) all f_i > 1.
At the largest special point R, f_j(R) = 1 and f_i(R) >= 1, so P_j(R) = prod_{i!=j} f_i(R) >= 1 with
equality iff all f_i(R)=1. Local expansions f_i = 1 + b_i s^{m_i} (s = t-R) have b_i of one sign
(all come from the same region), so P_j - 1 has order min_{i!=j} m_i, which must be >= m_j for all j
=> all m_i equal. Then all f_i lie in the same region (1,inf), (-1,1) or (-inf,-1) on the next
interval, and the same argument applies at the next special point (|P_j| <= 1 or >= 1 with equality
iff all |f_i| = 1; sign condition forces all f_i equal to E).  QED
Corollaries:
 (a) With a constant cofactor C >= 2 (fixed primes) there is NO family having any real special point
     (P_j(R) = C*prod >= C > 1). In particular no linear/odd-degree family with constant primes.
 (b) At most one linear member; if f_i - 1 = a_i*G for a common G, reduce to Lemma 3 (simple proof:
     for the largest b_i = 1/a_i, every factor 1 - 2b_i/b_j <= -1, so |prod| > 1).
 (c) A surviving family must have its +-1 special points mostly non-real (e.g. f_i > 1 on R).
 - Faithful BPSW in GP (src/fam/bpsw.gp): reproduces A001262 (<1e5) and A217255 (<1e5) exactly,
   all primes < 2e5 pass, no composite passes < 3e5.
 - Half-family search extended (src/fam/half.gp): d <= 150, cliques of size 3..5 under the
   pairwise unit filter; Fermat-identity families = ONLY {10,14,18}*2^i (i=0..3), i.e. the single
   family {10,14,18} composed with Y -> Y^2-2. All contain d divisible by 5 =>
   LOCAL OBSTRUCTION for Y = L_{2K}: Y = 5F_K^2 + 2(-1)^K == +-2 mod 5; K even gives n == 1 mod 5
   while 5 | z(p_1) needs n == -1 mod 5; K odd gives 5 | Psi_10(Y). Dead for all K.
   => Lucas-exponential Psi half-families: exhausted (k<=5, d<=150).
 - Aurifeuillian extension (src/fam/aur.gp): Lucas-exp shapes = irreducible factors (content removed)
   of Psi_d(5Z^2 -+ 2), Z = F_k (L_{2k} = 5F_k^2 + 2(-1)^k), d <= 60: 125 shapes, 0 symmetric
   unit-resultant edges => no half-family. (Psi_10(L_{2k}) always has content 5 here: same mod-5 obstruction.)
 - Numerical LM over C (src/alg/lm_williams.py): k=5 quadratic "solution" found =
   f_i = zeta5^i t^2 + c_i (b_i = 0, |c_i| equal): a complex Galois-orbit family, linear in u = t^2
   with root-of-unity slopes. Values are conjugate algebraic integers, not rational primes => unusable;
   consistent with Lemma 3 (which needs real same-sign slopes). Real-coefficient run: lm_real.py.
 - Literature (alphaXiv, 2026-10-17): no construction of a BPSW psp or Williams number found.
   Relevant tools: Alford-Grantham-Hayman-Shallue subset-product (1203.6664; regime |G| small, N huge),
   Fibonacci factors -> Fibonacci psp (2105.13513), absolute Lucas psp tabulation (2306.17691),
   Baillie-Fiori-Wagstaff strengthening (2006.14425).

## Lemma 7 (PROVED): no Williams family of monic integer quadratics (any k)
Register rows: C4, C3F, C3L, monic integer quadratics. Closes ONLY that sub-route.
Unit-resultant fact (any degree): if f_i in Z[t] monic and the Williams identities hold, then
prod_{i!=j} Res(f_j -+ 1, f_i) = prod over roots rho of f_j -+ 1 of P_j(rho) = 1, so every
Res(f_j -+ 1, f_i) = +-1 (integers).
Quadratics f = t^2 + b t + c: by Lemma 6 + Lemma 3, case "f_j - 1 has real roots" forces f_j = 1 + a_j G
for all j (Lemma 3 contradiction); f_j - 1 without real roots => f_j + 1 > 0 too. So all f_j > 1 on R.
Res(f_j - 1, f_i) = Q_j(x,y) := x^2 - b_j x y + (c_j - 1) y^2 with x = c_i - c_j + 1, y = b_i - b_j,
a POSITIVE DEFINITE form (disc b_j^2 - 4(c_j - 1) < 0). y = 0 forces f_i = f_j (also using the +1 side).
Q_j = (x - b_j y/2)^2 + (h_j - 1) y^2 with h_j = min f_j = c_j - b_j^2/4 > 1. Q_j = 1 with y != 0 forces
h_j - 1 <= 1 and then |y| = 1 (h_j = 2 for b_j even, h_j = 7/4 for b_j odd). Hence |b_i - b_j| = 1 for
ALL pairs i != j => k <= 2; but k=2 violates 3d <= D. QED
(Non-monic integer quadratics: Res(f_j - 1, f_i) only divides a power of a_j; not covered — LM search.)
 - CORRECTED Aurifeuillian search (src/fam/aur3.gp; filter = Res supported on primes of lc; global
   clique state): Lucas-exp k odd: 66 shapes, 168 edges, 1 hit; k even: 60 shapes, 213 edges, 2 hits —
   all hits are the SAME family {Psi10,Psi14,Psi18} (or its Y->Y^2-2 composite) at Y = 5Z^2 -+ 2,
   hence still mod-5 obstructed. Fermat-exp (Phi_d(X), Phi_d(2X^2)): 128 shapes, 66 edges,
   5 triangles, 0 hits. => no new exponential half-family.
 - Monic integer search (src/alg/monic_will.gp): all irreducible monic polys deg 2 (|coef|<=8),
   deg 3 (<=4), deg 4 (<=3): 2337 polys. Williams unit graph: 42 edges, 4 triangles, 0 full families
   (k<=6). Fermat-half graph: 761 edges; one maximal hit {Psi10, Psi14, Psi18, g}, g = t^4-3t^2-2t+1
   (NOT cyclotomic; g-1 = t(t-2)(t+1)^2: same special points {2,-1,0}). Polynomial half-families
   exist beyond Psi but are useless without Lucas control (which needs Psi shapes at Y = L_{2K}).
 - Groebner tooling: Singular k=3 + constant C job timed out (25 min). User directed: use Calyx
   (~/Desktop/calyx-math) for Groebner; build delegated to an agent (2026-10-17).
 - Direct scan (src/fam/halfscan.gp), all integers |Y| <= 2e5: family {Psi10,Psi14,Psi18}(Y):
   50 values pass STRONG Fermat base 2 (Fermat identity working as proved), 0 pass BPSW;
   4-member family {.., t^4-3t^2-2t+1}: 1 strong-Fermat pass, 0 BPSW. Lucas side is the binding constraint
   (p+1 = Y(Y-1), (Y-1)(Y^2-2), Y(Y^2-3); n+1 == 0 mod (Y^2-2) only; elsewhere needs index ~ Y^-6).
 - Calyx build (agent, 2026-10-17): BLOCKED by repo guard, 19 GB free < 20 GB (verified with df);
   needs user to free disk or waive guard. Singular remains the Groebner engine meanwhile.

## Attempt D — doubly-small-order (DSO) pool (2026-10-17)
Register rows: C1, C2, C5 (pool primes with ord, z <= 1500). Bounded evidence only.
Primes with ord_p(2) <= 1500 and z(p) <= 1500 = prime factors of gcd(2^a-1, F_b), a,b <= 1500
(src/dso/dso.gp, dso_primes.txt): 396 primes, ALL small (max 469561) — consistent with gcd bounds of
Bugeaud-Corvaja-Zannier type (gcd(2^a-1, F_b) is tiny). Frontier (src/dso/dso_frontier.py, use orders):
class (2,0) [100 primes]: (45.2 bits, 6), (243.4, 22), (277.9, 24); class (1,2): (115.0, 11);
classes (0,1),(1,1): nothing usable. N/log2|G| ~ 0.09 — same as CHM pools. Deprioritised as a pool source
(MEASURED for ord, z <= 1500 only; larger order bounds unmeasured). Not closed.

## Lemma 8 (PROVED): a real Williams family needs shared special points on BOTH sides
Register rows: C4, C3F, C3L (real coefficients; (a) also a_j > 0). A tool for Williams families only.
Setting: f_1..f_k in R[t] (k >= 3), degrees d_j, leading coeffs a_j, n = prod f_j, f_j -+ 1 | n -+ 1.
Let e_m = elementary symmetric polys of the f_j. Then
   E_- := n - 1 - prod(f_j - 1) = -1 - sum_{m<k} (-1)^{k-m} e_m,   E_+ := prod(f_j + 1) - n - 1 = sum_{m<k} e_m - 1.
Both have degree exactly D - d_min with top coefficient (prod a) * sum_{j in J_min} 1/a_j  (J_min = members of
minimal degree), unless that sum is 0.  Lambda_-+ := lcm_j(f_j -+ 1) divides E_-+.
 (a) All a_j > 0 (the case of prime values at large t): deg Lambda_-+ <= D - d_min, i.e. the OVERLAP
     omega_-+ := sum_j deg(f_j -+ 1) - deg lcm >= d_min on EACH side.
 (b) Any signs: if the f_j - 1 are pairwise coprime AND the f_j + 1 are pairwise coprime, then
     n - 1 = prod(f_j - 1), n + 1 = prod(f_j + 1), i.e. S_odd = 1, S_even = 0 (S = sums of e_m by parity of k-m);
     top coefficients give sum_{J_min} 1/a_j = 0 and e_2(1/a_j : J_min) = 0, so sum 1/a_j^2 = 0: impossible over R.
     (Over C it is possible: the zeta_5 family found by LM has sum zeta^-i = e_2 = 0.)  QED
Corollary (real LINEAR families, any signs): (b) kills the pairwise-coprime case, Lemma 3 the common-root case;
  numerics (src/alg/lm_lin.py, lm_lin_defl.py; k=4..9, 300-400 starts each) find only the trivial
  all-equal family (k odd). Mixed-sign linear families: no non-degenerate solution found.
Corollary (real QUADRATIC families, a_j > 0): by Lemma 6 + Lemma 3 every special point is non-real; a shared
  non-real root forces its conjugate, so sharing = proportionality f_i -+ 1 = lam (f_j -+ 1); (a) needs one such
  pair on each side, and they must be different pairs (same pair => f_i = f_j).
## Theorem 9 (PROVED, computer-assisted): no real quadratic Williams family with k = 4, a_j > 0
Register rows: C4, C3F, C3L, quadratic, k = 4, a_j > 0. Closes ONLY that sub-route.
Up to symmetry the patterns are A1: f2-1 = lam(f1-1), f3+1 = mu(f1+1); A2: f2-1 = lam(f1-1), f4+1 = mu(f3+1);
normalisation f1 = t^2 + c, remaining members free; nondegeneracy lam, mu, a_j != 0, lam, mu != 1.
src/alg/gen_quadpat.py -> qp_4_A1.sing, qp_4_A2.sing: Singular std = <1> (dim -1) in both, i.e. no solution
even over C.  Independent check: real LM over ALL k=4 quadratic families (lm_real.py, 400 starts) min residual
0.083, no solution.  (The full unreduced k=4 ideal k4q.sing was killed after ~35 min, exit 137.)
 - Calyx (agent build, rust@4a8d9faba26, 6m49s; binary ~/Desktop/calyx-math/target/release/calyx):
   k3_elim reproduces Singular exactly (all intermediate reduced bases identical, cmpgb.py; checked myself).
   FULL unreduced k=4 quadratic ideal k4q (any signs, complex, a_i != 0; f1 = t^2+c1): Calyx GroebnerBasis = [1]
   in 158 s (k4q_calyx.out/.log). STATUS: COMPUTATIONALLY OBSERVED, not proved — Calyx's modular algorithm
   accepts a result by Monte Carlo (one extra prime). Rigorous version for the case that matters (real,
   a_j > 0) is Theorem 9 above (exact Singular). Exact certificate (lift / non-modular run) still to do.
 - Side check: D = -7 "Lucas-exp" shapes in TWO variables (Y = V_M, Z = 2^M; G_d(Y,Z) = Phi_d(alpha^M, beta^M))
   cannot carry a Fermat-half identity: Res_Y(G_a - 1, G_b) = prod_zeta (Z^{phi(a)/2} G_a(c_zeta,1) - 1) has
   positive Z-degree, but the identity forces prod_{i!=j} Res_Y(f_j - 1, f_i) = 1 (monic in Y). PROVED.
   Register rows: C4, C3F, C6, D = -7. Closes ONLY all-forced D = -7 Lucas-exponential families.
 - k=5 real quadratic, pattern A1: Calyx GB = [1] in 30 s (qp_5_A1_calyx.log; Monte Carlo modular).
   Pattern A2 exceeded the 6 GB memwatch cap; rerun with pairwise-distinctness (random linear forms,
   qp_5_A2d) in progress. Singular runs of both were SIGKILLed under memory pressure.
   Encodings validated: the trivial all-equal family satisfies every generated equation; perturbed points fail.
 - Odd-degree families (a_j > 0): every f_j -+ 1 has a real root, so by Lemma 6 there are real special points
   of BOTH signs, shared by all members => k must be ODD (P_j(s) = (-1)^{k-1} = 1), and the real parts
   R_-+ of f_j -+ 1 are common to all j and proper (all-real => proportional => Lemma 3).
   Cubic k=5 system after normalising r=0, s=-1 (f_i = 1 + t g_i, f_i(-1) = -1): src/alg/gen_cubic.py,
   cub_5d (with pairwise distinctness) running in Calyx.

## Strategic review (2026-10-17, after Attempt D, Lemma 8, Theorem 9)
- Advanced: the polynomial-family route is now closed in every case tested — linear (same-sign proof;
  mixed-sign numerics), k=3 any degree, monic integer quadratics, real quadratics k=4 (exact), k=4 any-sign
  quadratics (Calyx, Monte Carlo), k=5 pattern A1 (Calyx). Exponential half-families: closed (cyclotomic,
  Aurifeuillian, D=-7 two-variable). Pools (CHM, DSO, D=-7): N/log2|G| ~ 0.1.
- Theory of success for an explicit n: NONE credible. The only producer at modest cost would be a real Williams
  family; every structural constraint found (Lemmas 6, 8) points against one. The Pomerance route needs
  |G| >~ 2^1000 with generic search cost 2^(B/log2 N) >> 2^100.
- Kill: no further pool enumeration at larger y (cannot change search cost). Continue (bounded): finish the k=5
  quadratic/cubic eliminations; attempt a GENERAL proof "no real Williams family with a_j > 0" (open; the
  obstruction is complex special points, where Lemma 6's ordering argument does not apply).
- What would change the decision: a real family found in any computation; or a subset-product algorithm
  exploiting structure in G = prod (Z/q)^* beyond Wagner's (none known).
 - Cubic k=5 numerics: (i) LM on the pseudo-remainder equations gave 42 "zeros", ALL spurious (prem scales by
   a_i^m, so a_i -> 0 fakes a solution; direct monic-division check fails by 1e6..1e14). (ii) LM with true
   remainders + barriers (lm_cubic2.py): near-zeros only in the degenerate limit |a_i| -> infinity with
   p_i/a_i -> ~1 (roots of g_i, h_i collapsing onto the shared points 0, -1). No genuine family seen yet.
   Exact cubic and k=5 pattern-A2 eliminations exceed the 6 GB memwatch cap (Calyx over Q and over GF(2^31-1)).
   STATUS: open (no evidence of a real cubic family).

## Item (i) of the stop checkpoint re-examined: dense subset-product algorithms (AGHS 1203.6664, read 2026-10-17)
- AGHS Algorithm 2 (Kuperberg / Flaxman-Przydatek pairing over a subgroup chain): time and space ~4^sqrt(log2|G|),
  but it NEEDS N >~ 4^sqrt(log2|G|) * log|G| symmetric pool elements. Its Carmichael success comes from
  ONE-SIDED divisor pools {d+1 prime : d | Lambda} (N ~ tau(Lambda)/log, about 1e10), where each p is trivial
  modulo Q_i with probability >= 1/(h_i+1), so most coordinates can be filtered out (subgroup trick).
- BPSW needs TWO-SIDED pools (p-1 A-smooth and p+1 B-smooth). The divisor trick does not apply: d | L1 with
  d+1 | L2 is an S-unit-type condition, so N collapses. For primes p <= k^c with p^2-1 k-smooth:
  N <~ k^c rho(c)^2 / log  (maximised over c: about e^{2 sqrt(k)/e}), while log2|G| ~ c k. Kuperberg/Wagner
  need N >= 4^sqrt(log2|G|) ~ e^{1.4 sqrt(c k)}. Also each p has only omega(p^2-1) ~ log trivial coordinates,
  so the subgroup filter is useless.  HEURISTIC CONCLUSION: at NO scale k does a Pomerance-type BPSW pool
  admit a solution findable by known dense subset-product algorithms (existence holds only for large k).
  Item (i) is DEPRIORITISED by this heuristic, NOT closed (relabelled 2026-10-17: a heuristic may rank a route,
  never close it).

## Novelty audit (2026-10-17; web + alphaXiv search, primary sources read where accessible)
Context: the polynomial-family lemmas concern WILLIAMS NUMBERS (Carmichael n with (D/n)=-1, p+1 | n+1 for all p),
H.C. Williams, "On numbers analogous to the Carmichael numbers", Canad. Math. Bull. 20 (1977) 133-143 (READ, pp.142-143):
existence open (still open per McIntosh, Integers 14 (2014) #A59, read). Williams PROVES: such n has k odd, k >= 5,
by an integer inequality argument using only p_j -+ 1 | n -+ 1 (no primality beyond distinct odd integers).
 - Lemma 4 (k=3, any degree): the positive-leading-coefficient case FOLLOWS from Williams 1977 (apply his k != 3
   argument to integer values at large t). New content is only the mixed-sign/real case and the explicit c-locus
   (c1+c2+c3=3, c1c2c3=1; anisotropic torus) -- not found in the literature. Cite Williams.
 - "k must be odd" (Lemma 6, cubic remarks): known for Williams numbers (Williams 1977, via eps(n) = (-1)^k).
 - Lemma 6(a) (fixed cofactor C>=2 kills families with a real special point): cf. McIntosh 2014: finitely many
   super-Carmichael numbers (p-+1 | N-1, a DIFFERENT condition) with d-3 factors fixed. Analogous flavour; no
   polynomial statement found. Polynomial Williams-family statements (Lemmas 3, 6, 7, 8, Thm 9): no prior source found.
 - Lemma 5 is a CARMICHAEL family (Korselt: p_j - 1 | n - 1 in Z[Y]), not merely base-2 Fermat. Checked (sympy):
   76 values Y in [3,20000) with all three prime, Korselt holds for all; Y=3 -> 5*13*17 = 1105, Y=9 -> 71*631*701.
   Beeger-Duparc parameters (Jameson notes, Prop 1.3): h1=(Y^2-Y-1)(Y^2+Y-1), h2=Y^2-2, h3=(Y-1)^2, Delta=h2h3-p^2=2Y-3.
   So it is one polynomial solution of the classical Beeger-Duparc parametrisation; nonlinear polynomial C3 families of
   this kind are plausibly folklore. Only linear (Chernick) families found in print. NOVELTY: UNRESOLVED (OEIS lookup
   blocked by 403; check OEIS for 1105, 31405501, 6047866621 manually).
 - Not accessible (403), closest unchecked risk: Zhang, "Estimating the counts of Carmichael and Williams numbers with
   small multiple seeds", Math. Comp. 84 (2015) 309-337; Ballot (Integers 2026) Chernick-type universal forms for
   absolute Lucas pseudoprimes (linear, Lucas side only).
 - msolve 0.10.1 cross-checks (src/alg/*.ms, mod 65521 unless noted): qp_4_A1 over Q: [-1] (no solution, 0.02 s,
   agrees with Singular). Positive control qp_5_A1_ctl (trivial family allowed): [1,13,-1,[]] = infinitely many,
   as expected. qp_5_A1: [-1] (agrees with Calyx).
 - Parameter count (any k, any class pattern, real quadratics): unknowns 3k-2; proportionality gives 2(k-A)+2(k-B)
   real equations (A, B = numbers of -/+ classes), Williams gives one complex equation per class = 2A+2B. Total 4k,
   so the system is overdetermined by k+2 for EVERY sharing pattern (degree d: by k(d-1)+2). The only known
   solutions are the trivial family f_i = f (k odd) and its Galois twists over C (zeta_k-family: n = u^k,
   f_i = zeta^i u). A real family would need an unexplained coincidence.
 - qp_5_A2d (pattern A2 with pairwise distinctness) and cub_5d (cubic k=5 with distinctness): msolve mod 65521
   did not finish in 600 s; Calyx exceeds the 6 GB cap. The Rabinowitsch distinctness variables (8-10 extra)
   are what makes these hard. Open; low priority given the overdetermination count above.

## State after this round (2026-10-17) — SUPERSEDED by "Checkpoint, end of session 1" below
CONSTRUCTIVE: no mechanism left that can output an explicit n.
  - Pomerance/Erdos subset products: existence only at scales where every known dense subset-product
    algorithm (Wagner, Kuperberg/AGHS) needs N >= 4^sqrt(log2|G|) pool primes, while two-sided pools give
    N <~ e^{0.74 sqrt k}; never feasible (HEURISTIC, quantified above).
  - Polynomial (Williams) families: excluded rigorously for linear same-sign, k=3 (any degree),
    monic integer quadratics, real quadratics k=4; Calyx/msolve (randomized) for all k=4 quadratics and k=5
    pattern A1; overdetermined by k(d-1)+2 in general.
  - Exponential half-families, small-order pools, D=-7 coupling, two-variable Lucas shapes: excluded.
OPEN MATHEMATICAL QUESTION worth pursuing: prove that no real polynomial Williams family with positive
leading coefficients exists (all k, all degrees). Known ingredients: Lemma 6 (real special points), Lemma 8
(complex sharing >= d_min), Lemma 3 (single shared G). Missing: a positivity argument at NON-real special points.

## Checkpoint, end of session 1 (2026-10-17) — READ THIS FIRST WHEN RESUMING
### Last results
 - k=5 pattern A2 with the free member eliminated linearly (src/alg/gen_A2lin.py -> A2lin.sing, A2lin.ms):
   Singular eliminate hit the 1200 s timeout with EMPTY output (A2lin.out, 0 bytes); msolve mod 65521 killed
   after the peer review below. UNFINISHED and CUT (no escalation). Same for qp_5_A2d, cub_5d.
   [COMPUTED: src/alg/A2lin.out (0 bytes), qp_5_A2d*.msout/_calyx.log, cub_5d*.msout/_calyx.log]
 - No candidate n. No submission JSON written.

### Peer review accepted (session "Lead agents lemmas/theorems", 2026-10-17)
 - Lemmas 3, 4, 6-9 and Theorem 9 all assume the bounded-index identity f_j -+ 1 | n -+ 1, i.e. they are about
   WILLIAMS NUMBERS (Carmichael AND p+1 | n+1), a stronger problem than BPSW [CHECKED by me: register tags on each
   lemma above, rows C3F + C3L]. Lemma 4, positive-leading case,
   follows from Williams, Canad. Math. Bull. 20 (1977) 133-143, pp.142-143 (k odd, k >= 5): cite it.
   [REVIEWER-REPORTED: the reviewing session read pp.142-143; I have NOT read the paper myself.]
   DECISION: no more obstruction work in the Williams framework (k=5 quadratic/cubic eliminations stopped).
 - Labels: Calyx GB = [1] (k4q, qp_5_A1) and msolve mod-p [-1] are Monte Carlo/modular = COMPUTATIONALLY
   OBSERVED, not theorems. Exact certificates are low priority now that the framework is cut.
   Theorem 9 rests on exact Singular std over Q (qp_4_A1.out, qp_4_A2.out: dim -1).
 - Lemma 5 is a CARMICHAEL family for Y >= 3 (Korselt in Z[Y]; Beeger-Duparc h2 = Y^2-2, h3 = (Y-1)^2, Delta = 2Y-3)
   [PROVED: identity; Beeger-Duparc form CHECKED by me, see Lemma 5 entry]: the WEAK Fermat test holds for every
   base coprime to n. The STRONG base-2 test still needs one 2-adic class (F1), and the Lucas half binds.
 - What BPSW actually needs (Lemma 1): odd(ord_p 2) | n-1 and odd(z(p)) | n+1, with no bound on the index
   (p-1)/ord_p 2 or (p+1)/z(p) [PROVED: Lemma 1]. Relaxed target ("relaxed Williams family"): f_j - 1 | n - 1 by identity,
   and g_j | n + 1 for a factor g_j of f_j + 1 with z(f_j) | g_j FORCED (as Psi_d(L_2K) forces z | 2Kd).

### Lemma 10 (PROVED from Lemma 1, CRT): core + extension decomposition
Register rows: F1, F2, F3, C1 (any Selfridge D works the same way). Equivalence; closes nothing.
n = n_F * s squarefree, (5/n) = -1, all prime factors in one 2-adic class (j1, j2) with v2-room.
A_F = lcm_{p | n_F} odd(ord_p 2), B_F = lcm_{p | n_F} odd(z(p)); A_S, B_S the same for s. Then
   n is a BPSW psp  <=>  s == n_F^{-1} (mod A_F),  s == -n_F^{-1} (mod B_F),  n_F s == 1 (mod A_S),  n_F s == -1 (mod B_S).
Necessary: gcd(A_F A_S, B_F B_S) = 1 and gcd(n, A B) = 1.
Consequence (direct): for a Carmichael core (A_F | n_F - 1) with Lucas defect delta_L = B_F / gcd(B_F, n_F + 1),
the extension needs s == 1 mod lcm(A_F, gcd(B_F, n_F+1)), a modulus of size >= A_F [PROVED]. That A_F is about the
size of the core primes was UNMEASURED; now MEASURED on the Lemma 5 family: log2 A_F / log2 n has median 0.80 (range
0.16-0.97) [COMPUTED: src/fam/lemma5_profile.out], i.e. A_F is even larger, about the whole core. So
a SMALL delta_L does NOT make a hybrid cheap: the cost is set by A_F, not by delta_L. A hybrid is cheap only if (i) delta_L = 1 exactly (pure family member, s = 1), or (ii) the core primes
have small orders (then the core is just part of an Attempt-A pool, no gain), or (iii) every extension prime
satisfies q == 1 mod lcm(A_F, gcd(B_F, n_F+1)) individually and has structured q -+ 1 (Erdos-style with L divisible by A_F).
Since A_F is large (measured above), this corrects the peer's suggestion "if the defect is a few bits, a tiny subset
product closes it".

### Lemma 5 property profile (G7; 2026-10-17) [COMPUTED: src/fam/lemma5_profile.gp -> lemma5_profile.out]
Register rows: C4, C3F. All Y with |Y| <= 200000 and |Psi10|, |Psi14|, |Psi18| distinct primes >= 3: 682 rows
(Y > 0: 337, Y < 0: 345; Y = -2, n = 105 skipped since 3 | Q).
 - Validation: the Lemma 1 prediction from per-prime orders (one 2-adic class + odd orders dividing n -+ 1) agrees
   with the direct strong tests on all 682 rows, across 13 different Selfridge D. So Lemma 1 holds for general D
   in practice, with z = rank of apparition [COMPUTED].
 - Selfridge D actually selected: D = 5 for only 153/682 (Y > 0: 75/337); -7: 204, -11: 125, 17: 80, -19: 51, ...
   So Attempt C's D = 5 analysis covers a minority of the family.
 - Y > 0: Korselt 337/337 (Carmichael, Fermat defect 1 always); strong Fermat(2) = one Fermat 2-adic class = 50/337;
   one Lucas 2-adic class 47/337; Lucas defect 1: 0/337; strong Lucas 0; BPSW 0.
 - Y < 0: Korselt 0/345; instead |Psi14|+1 and |Psi18|+1 divide n-1 for all of them [COMPUTED: lemma11_check.out];
   Fermat defect 1: 0/345. Not a Carmichael branch.
 - Lucas defect delta_L: log2 delta_L / log2 n has median 0.48, range 0.185-0.611. Small defects only at tiny n:
   Y = -3 (2.3 bits), 3 (3.2), -4 (7.1), 17 (7.8); for n > 2^33 every defect is >= 21.6 bits (Y = 66), typically about n^0.48.
   ANSWER to the reviewer's question "is the Lucas defect ever small?": NO, not beyond tiny Y [COMPUTED].
 - Indices: (p - (D/p))/z(p) <= 2 for 1322 of 2046 primes, <= 10 for 1726; (p-1)/ord_p 2 <= 2 for 1160, <= 10 for 1640.
   The orders are generic; nothing in this family forces a large index.

### Lemma 11 (PROVED, elementary; checked numerically): why the Lucas defect cannot vanish in the Lemma 5 family
Register rows: F1, C3F (Korselt), C4; ANY Selfridge D.
 (a) General: if n is a Carmichael number and a strong Lucas psp for D, then every prime p | n that SPLITS for D has
     z(p) a power of 2. Proof: z(p) | p - 1 | n - 1 and odd(z(p)) | n + 1, so odd(z(p)) | 2.
 (b) Lemma 5 family, Y >= 3, n = Psi10 Psi14 Psi18 (all prime) a strong Lucas psp for its D. Then
     z(Psi10(Y)) = z(Psi18(Y)) = 2^j2, i.e. Psi10(Y) * Psi18(Y) divides V_{2^(j2-1)} for that D.
     Proof: split case by (a). Inert case: odd(z(p)) | gcd(p+1, n+1). In Z[Y], n(0) = 1 and n(1) = -3, so
     n+1 == 2 mod Y and == -2 mod (Y-1); hence gcd(Psi10+1, n+1) = gcd(Y(Y-1), n+1) divides 2. Modulo Y^2-3:
     Psi10 = 2-Y, Psi14 = Y-2, Psi18 = -1, so n+1 == 4(2-Y), and gcd(2-Y, Y^2-3) = 1 since Y^2-3 = (Y-2)(Y+2)+1;
     hence gcd(Psi18+1, n+1) = gcd(Y(Y^2-3), n+1) divides 8. So odd(z) = 1 for both; Lemma 1 (one Lucas class)
     makes the two 2-powers equal.  QED
     [COMPUTED check, src/fam/lemma11_check.gp -> .out: both gcd bounds and (Y^2-2) | n+1 hold on all 337 rows;
     no row has odd(z) = 1 for Psi10 or Psi18.]
 Only Psi14 has a polynomial factor of p+1 dividing n+1 (Y^2-2); Psi10 and Psi18 have none.
 WHAT A WORKING OBJECT MUST LOOK LIKE (G2 extraction):
  - In a pure relaxed family (S1), EVERY member f_j needs a polynomial factor g_j of f_j + 1 (inert) with
    deg gcd(g_j, n + 1) large in Q[Y]. Otherwise z(f_j) is forced to be a power of 2 by the argument above. This is
    a polynomial condition to impose directly in the Beeger-Duparc search (R1), not a measurement.
  - Lemma 5 itself survives only through Y with Psi10(Y) * Psi18(Y) | V_{2^m} for the selected D: a finite search
    per (D, m) over prime factors of V_{2^m} of the shapes Y^2-Y-1 and Y^3-3Y-1 with the same Y.

### Sharp sufficient conditions (targets; each would output n)
 (S1, pure) A Carmichael polynomial family f_1..f_k in Z[Y] (Korselt identity) and a parametrisation Y = Y(K)
   (exponential or polynomial) with polynomials g_j such that z(f_j(Y(K))) | g_j(K) is forced and
   g_j | n + 1 holds identically; then any K with all f_j(Y(K)) prime gives n (after Selfridge-D and
   2-adic class checks) [PROVED from Lemma 1]. Attempt C = this with g_j Lucas-exponential; the only Psi families
   found ({10,14,18}*2^i, d <= 150, k <= 5 [COMPUTED: src/fam/half.out]) fail for D = 5 at every K by the mod-5
   argument [PROVED: Attempt C].
 (S2, extension) A Carmichael core and a pool P of primes q == 1 mod M_F = lcm(A_F, gcd(B_F, n_F+1)) with
   odd(ord_q 2) | L1, odd(z(q)) | L2: n exists iff some subset of P has product hitting the target of Lemma 10 in
   G = (Z/delta_L)^* x (Z/L1)^* x (Z/L2)^* [PROVED from Lemma 10]; count 2^|P| / |G| [HEURISTIC];
   DEFICIT = log2|G| - |P| bits. Not yet measured for any extension pool (Attempt-A pools: 241.7 bits at z = 150,
   class (2,0), point (271.7 bits, 30 primes) [COMPUTED: data/frontier_z150.json]).

### Live constructive directions (>= 3), each with first experiment and kill condition
 R1 Pure relaxed families (S1). [First measurement DONE 2026-10-17: profile + Lemma 11 above. NEXT: Beeger-Duparc
    search with the Lemma 11 condition imposed: every member's f_j + 1 must share a high-degree factor with n + 1.]
    Original first step: along Lemma 5 (the 76 Y in [3,20000) with all three values prime, Korselt holding for all 76 [CHECKED by me:
    src/fam/lemma5_count.gp; Y = 3, 9, 17 give 5*13*17, 71*631*701, 271*4591*4861], plus negative Y)
    compute the Selfridge D actually chosen, the 2-adic classes, A_F, B_F, delta_L and whether B_F or delta_L factor
    polynomially in Y. Then enumerate low-degree Beeger-Duparc polynomial families (h2, h3 of degree <= 2,
    Delta(Y) | (p-1)(p+h3) identically) and repeat. Kill: delta_L shows no structural (polynomial or
    Lucas-exponential) dependence on Y in >= 3 families, so no forced-z substitution exists for them.
 R2 Mixed families: some members Lucas-forced (Psi_d(Y), Y = L_2K for D=5), the others free polynomials in Y
    chosen so that the whole product satisfies Korselt in Z[Y]; only the free members contribute to delta_L.
    Other D: for D = -7 (Q = 2) the all-forced version is DEAD (PROVED: two-variable resultant obstruction, Theorem 9
    side check, Z = 2^M). Other Q != +-1: UNCHECKED (the same argument plausibly applies; verify before relying on it).
    The mixed version is NOT covered by that proof for any D. First: extend src/alg/cyc_sets.gp to search for
    free Z[Y] members completing a Korselt identity with 1-2 Psi_d members, d <= 60. Kill: no completion at degree <= 4.
 R3 Carmichael-extension hybrid (S2, case (iii) of Lemma 10). First: for the smallest Lemma-5 cores
    (Y = 9: 71*631*701, etc.) compute A_F, B_F, delta_L, M_F; count primes q = 1 + m*M_F (m smooth, q < 2^64)
    with small odd(z(q)) and odd(ord_q 2), and measure the deficit of S2 against Attempt A at equal |G|.
    Kill: the deficit is no better than CHM pools at the same log2|G|.
 Closed BY PROOF, only for the stated sub-route (each uses a CHOSEN restriction, see Assumption register):
   Williams-identity families that are linear same-sign, k=3 (any degree), monic integer quadratic, real quadratic
   k=4 with a_j > 0 (Lemmas 3, 4, 7, Thm 9); D=-7 all-forced two-variable Lucas shapes; Lucas-exponential triples
   at Y = L_2K for D = 5, every K, for the Psi families found (mod-5 argument, Attempt C).
 Bounded computation only: Fermat-exponential pairs a,b <= 160 [COMPUTED: src/alg/cyc_sym.out]; Psi_d
   Fermat-identity sets d <= 150, k <= 5 [COMPUTED: src/fam/half.out]; k=4 and k=5-A1 quadratic ideals
   (Calyx/msolve, Monte Carlo [COMPUTED: src/alg/k4q_calyx.out, qp_5_A1*.out]).
 Deprioritised by measurement/heuristic, NOT closed: Pomerance pools (CHM, DSO, D=-7; measured N/log2|G| about 0.1
   at z <= 150; infeasibility at all scales is HEURISTIC); further obstruction work in the Williams framework.

### Resume instructions
 - Tools: src/bpsw.py (check_submission(path) = verifier), src/fam/bpsw.gp (sprp2, selD, lucUV, slprp, bpsw),
   GP/PARI, Singular, msolve 0.10.1, Calyx (~/Desktop/calyx-math/target/release/calyx; memwatch caps any calyx
   process at 6144 MB, user config, do not override). Converters: src/alg/sing2ms.py, src/alg/sing2calyx.py.
 - R1's first measurement is DONE (Lemma 5 profile, Lemma 11). Next: R1 Beeger-Duparc search with the Lemma 11
   gcd condition imposed (constructive, per G2), then R3, then R2. Optional cheap side search: Y with
   Psi10(Y) * Psi18(Y) | V_{2^m}.
 - Process state at checkpoint: no task-owned solver or worker running (memwatch is the user's persistent service).

## Session 2 (2026-10-24) — new results
### Lemma 12 (PROVED, elementary): constraints from the selected D
Register rows: F1. Let D* be the Selfridge D of n (the list is every D == 1 mod 4 with |D| >= 5, ordered by |D|).
 (a) Which D can be selected. Write D* = s^2 D0 with D0 squarefree, a product of prime discriminants l* (l* = +-l == 1 mod 4).
     Every l* with 5 <= |l*| < |D*| comes earlier in the list, so (l*/n) is 0 or +1 for those. Hence
     D* is in {l* : l >= 5 prime} or {-3 l* : l >= 5 prime} or {-3 s^2 : s >= 3 odd}.
     Proof: if D0 had two prime-discriminant factors of size >= 5, both are +1 (a 0 gives (D*/n) = 0), so (D*/n) = (-3/n) or 1.
     Then -3 l1* (or l1*) comes earlier with the same symbol, and s > 1 repeats an earlier D0.  QED
 (b) mod-l lemma. If |D*| = l is prime then (D*/n) = (n/l) = -1, so n is not 1 mod l. Hence l does not divide ord_p 2 for
     any p | n. If moreover D* > 0 (l == 1 mod 4), n is not -1 mod l either, so l does not divide z(p).
     For D = 5: NO prime factor of a BPSW psp has 5 | ord_p(2) or 5 | z(p). This kills every Lucas-forced member
     Psi_d(L_2K) with 5 | d, for every K, in one line. That includes the whole {10,14,18}*2^i family for D = 5, and
     subsumes Attempt C's mod-5 argument. [PROVED]
### Lemma 13 (PROVED): one-side-automatic products and an equal-orders criterion (D = 5)
 (a) If all p | n (squarefree, >= 2 primes) have the same Fibonacci entry point z(p) = N and (5/n) = -1, then n is
     a strong Lucas psp (Selfridge D = 5). Proof: p == (5/p) mod N (split: N | p-1, inert: N | p+1), so
     n == (5/n) = -1 mod N, so N | n+1 and every p has the same v2(z) <= v2(n+1).  QED
     (Dual, classical: any product of >= 2 primes with the same ord_p 2 = N is a strong psp base 2.)
 (b) Hence: >= 2 primes with the SAME pair (ord_p 2, z(p)) = (N1, N2), gcd(N1, N2) | 2, an odd number of them inert
     ==> their product is a BPSW psp. (5 cannot divide N1 N2 by Lemma 12(b).) This is the "doubly small order"
     coincidence. Attempt D's 396 DSO primes would have shown such a pair as a tiny subset-product solution.
### F3 loophole (non-squarefree n) checked [COMPUTED: src/wief/lw.gp]
p^2 | n needs p Wieferich (base 2) AND z_D(p^2) = z_D(p) for the selected D ("Lucas-Wieferich"). Known Wieferich primes:
 - 1093: Lucas-Wieferich for D in {1065, 4877, 5485, 8777, -13451, -14103, 15105, 18193} (|D| <= 20000).
 - 3511: Lucas-Wieferich only for D = 7205.
By Lemma 12(a) only prime discriminants can be selected: D* = 4877 or -13451 for 1093, nothing for 3511.
Selecting D* = 4877 needs (l*/n) in {0, +1} for every prime discriminant with |l*| < 4877: about 650 Jacobi conditions.
That adds about 650 bits of parity constraints to any construction. STATUS: the loophole is real, but with the known
Wieferich primes it is unusable at feasible size. A Wieferich prime p that is Lucas-Wieferich for a small prime
discriminant D would reopen it; expected first hit near |D| ~ 2p.
### Korselt/Williams identities as zero-sum problems (PROVED reformulation; COMPUTED classification)
Register rows: C4, C3F (+C3L for Williams). For monic f_j and n = prod f_j: lcm_j(f_j - 1) | n - 1 iff n(theta) = 1 to order
m(theta) at every root theta of every f_j - 1 (m = max multiplicity). The member with f_j(theta) = 1 contributes
nothing, so the condition is uniform in j. For f_j = Psi_d these are cyclotomic units, and the condition becomes
additive: unit exponents (and torsion) of Psi_d(theta_e) sum to 0. At double roots, also sum Psi_d'/Psi_d(theta_e) = 0.
Williams side: n(theta) = -1 at the roots of f_j + 1 (torsion bit 1), and the two special sets must be disjoint.
 - Tools: src/zs/vec2.gp (unit/derivative data), src/zs/zs2.py (DFS with unit-compatibility pruning; exact).
 - Positive control: recovers {10,14,18}, {20,28,36}, {40,56,72}, {80,112,144}.
 - RESULT [COMPUTED: src/zs/vec2_1200_12.out]: for d <= 1200 with all special points in fields of degree <= 12
   (e <= 400), the ONLY Korselt Psi-families (k <= 9) are {10,14,18}*2^i, i <= 3. There is NO Williams Psi-family (k <= 12).
   Bounded computation, not a proof for larger d or degree.
 - Consequence: for D = 5, the all-forced Lucas route (S1 with Psi members at Y = L_2K) is DEAD in this range,
   because the only family contains 10*2^i (Lemma 12(b)).
 - Mixed families (src/zs/zs3.py): Lucas-forced members (type F, 5 not | d) together with Williams members (type W,
   Psi_d + 1 | n + 1), both types present, k <= 9: NONE in the same three data ranges [COMPUTED].

### Strategic review (2026-10-24, after Lemmas 12-13, the F3 check and the zero-sum classification)
 - Distance to result: the polynomial and exponential family line is now closed IN RANGE on every axis.
   Pure Korselt Psi-families: only {10,14,18}*2^i. Williams Psi-families: none. Mixed forced/Williams: none.
   All-forced D = 5: impossible for that family by Lemma 12(b). Q != +-1: no one-variable Psi structure
   (Q = -t^2 or Q = 4 gives homogeneous two-variable forms, and the "-1" breaks the identity).
   Exactness: every closure above is either a proof or a bounded computation, tagged as such.
 - Theory of success: none credible for an explicit n at laptop scale.
   (i) Subset products are capped by the two-sided S-unit structure: for fixed y the pool {p : p^2 - 1 y-smooth}
   is finite, and known high-density algorithms need N >= 2^(c sqrt B). Index relaxation (cofactors of p -+ 1
   absorbed by the order index) gains only about 3x [HEURISTIC estimate, see below].
   (ii) Exponential families give finitely many candidates per family, and only one family exists in range.
   (iii) Non-squarefree n is blocked by Selfridge ordering.
 - Index-relaxation estimate [HEURISTIC]: allowing p - 1 = 2Au, p + 1 = 4Bu' (A, B smooth, u, u' <= U, 2 a u-th power,
   gamma a u'-th power) multiplies the pool by about 2-3x at y = 100, p ~ 1e12, U = 1000 (per side: P(u = 1 | cofactor
   <= U) ~ 0.15 versus E[1/u] ~ 0.24).
   That is not enough against the ~170x gap at 300 bits recorded in Attempt A. Deprioritised; not run.
 - Kill: further Psi/zero-sum scaling (more d, more fields). Exponential half-families. Williams-identity searches.
 - Keep (live): see "Live directions (session 2)" below.

### Literature check (2026-10-24, alphaXiv; looking for tools or constructions only)
 - Helmreich-Webster, 2306.17691 [READ by me]: tabulate ABSOLUTE Lucas pseudoprimes (Williams' criterion p - (5/p) | n - (5/n))
   for d = 5 up to 2^64 (4,247,414 numbers; most have k = 2). Two-factor ones are exactly n = p(p+2) with (5/p) = -1 and
   (5/p+2) = +1, so n + 1 = (p+1)^2 and the Lucas side is automatic. But the base-2 side forces odd(ord_p 2) | gcd(p-1, n-1) | 2,
   i.e. p = 3 [PROVED here, two lines]. Split primes in any absolute Lucas psp with (5/n) = -1 have p - 1 | n + 1, so for
   spsp(2) they must have ord_p 2 a power of 2 (Fermat-number factors) [PROVED here].
 - Hall, 2601.19817 [READ]: 200 Arnault composites (~350 bits, spsp to bases <= 11) all fail strong Lucas, with no near-misses.
 - Draziotis-Martidis-Tiganourias, 2002.07095 [READ]: birthday attack for subset products (2^(n/2)); nothing structural.
 Conclusion: no construction or algorithm beyond what the notebook already uses.

### Live directions (session 2), ranked by plausibility x testability
 P1 Pool-growth discriminating experiment (the only route with guaranteed existence). Extend the CHM closure to y = 200
    (then 250 if feasible). Measure the class-restricted inert pool N(y) and the full-pool group size B(y), and fit
    log N against sqrt(B). Hypothesis H: N(y) eventually exceeds 2^(c sqrt B(y)) for the best subset-product algorithm
    (Wagner / Kuperberg-AGHS). If H holds, compute the crossover scale and its cost; if it fails, the Pomerance route has
    no theory of success with known algorithms. First step: time src/chm at y = 200.
 P2 Close the exponential-family question fully for D = 5. Twisted shapes P_d(Y) = i^-deg Psi_d(iY) (4 | d) at Y = L_M,
    M odd, and Aurifeuillian factors of L_5m, both run through the zero-sum tool. Cheap; low plausibility.
 P3 Algorithmic: handle the 2-primary part of G (and any Jacobi/F_2 constraints) by linear algebra inside a k-tree
    pipeline, so the effective B drops by the 2-rank (~20-30% of B at y <= 150). Test on synthetic instances; matters
    only together with P1.
 P4 (theory, low priority) Uniqueness of Korselt Psi-families: are they exactly the doubling orbit {10,14,18}*2^i?
 P5 (passive) F3: a new Wieferich prime that is Lucas-Wieferich for a small prime discriminant would reopen p^2 | n.

## Session 2, continued (2026-10-24)
### Lemma 12(c) (PROVED; checked): the selected D is a prime discriminant or -15
Register rows: F1, F3. For a BPSW psp n the Selfridge D* is a prime discriminant l* (l >= 5) or D* = -15 (Q* = 4).
Proof. By 12(a) the other cases are D* = -3s^2 (s >= 3 odd) and D* = -3l* (l >= 7).
 - In both, (D*/n) = -1 forces (-3/n) = -1 (for -3l*: l* precedes D*, so (l*/n) = 1).
 - Then for every prime 5 <= l' below s^2 (resp. below l), both l'* and -3l'* precede D* and
   (-3l'*/n) = -(l'*/n). Neither may be -1, so (l'*/n) = 0, i.e. l' | n.
 - D* = -3s^2: Q* = (1+3s^2)/4 is odd, prime to 3 (Q* == 1 mod 3), and 1 < Q* < s^2. Every prime factor of Q*
   divides n, so gcd(n, Q*) > 1. Mod such a prime U_k == V_k == 1 for k >= 1, so the strong Lucas test fails.
 - D* = -3l*, l >= 11: 5 | n and 7 | n, but v2(ord_5 2) = 2 and v2(ord_7 2) = 0, contradicting F1.
   l = 7: D* = 21, Q* = -5 and 5 | n, so the Lucas test fails.  QED
 Side result (same argument): if D* = l* with l >= 23, then (-3/n) is 0 or 1.
 [COMPUTED: src/lem/l12c.gp. 1,282,319 odd composites up to 3e6. When D* is not a prime discriminant it is -15 (26,785 n),
 or -27/21 with gcd(n, Q*) > 1 (604 / 3,514). No violations.]

### Lemma 14 (PROVED, elementary): the criterion with the cofactor m = n/p
For squarefree n and p | n with m = n/p, the odd-part conditions of Lemma 1 at p are equivalent to:
 - p inert: odd(ord_p 2) * odd(z(p)) | m - 1;
 - p split: odd(ord_p 2) | m - 1 and odd(z(p)) | m + 1.
Proof: n - 1 == m - 1 mod (p - 1). Also n + 1 == 1 - m mod (p + 1) (inert) and n + 1 == m + 1 mod (p - 1) (split).
 Consequences:
 - Last-prime lemma. Every prime q | n divides gcd(2^(2^j1 odd(m-1)) - 1, U_(2^j2 odd(m -+ 1))), a number
   computable from m and D*. Since gcd(2^a - 1, F_b) is small (Attempt D; Bugeaud-Corvaja-Zannier type), the
   other primes determine q up to a short list. A search over cofactors m <= M0 covers every n with n/q <= M0,
   with q unbounded. Cost is about M0^2 bit operations. Estimate: for M0 <= 1e6 this mostly re-covers n < 2^64.
   NOT RUN; low priority.
 - D* in {-7, -15} (Q* = 2 or 4, so N(alpha) is a power of 2). At every p (inert or split), the odd-part
   conditions are equivalent to alpha^n == conj(alpha) mod p up to 2-power torsion. For inert p this reads
   odd(ord_p alpha) | m - 1. So BPSW with D* = -7 is the quadratic Frobenius test for x^2 - x + 2 (odd parts) plus the
   2-adic classes. Restatement only: the count is unchanged, since ord_p alpha has size about p^2.

### Exponential families for D = 5: the remaining cases closed in range
 - Forced Psi members plus ONE free Williams member f in Q[Y] (src/fw/fw1.gp, fw2.gp).
   Reformulation [PROVED]: f - 1 = c D1, f + 1 = c D2, where D1, D2 are monic divisors of P - 1 (P = product of the
   forced Psi), D2 - D1 = kappa is a constant, c = 2/kappa, and lcm(Psi_d - 1) | fP - 1. This is a finite check per
   forced set. If all special points of the forced set are roots of D1, the forced set is itself Korselt.
   The divisor-pair step was checked on two hand-built cases (src/fw/ctl.gp).
   RESULT [COMPUTED: src/fw/r3_60.out, fw2_r4_30.out]: 0 families for r = 2 (d <= 60), r = 3 (d <= 60) and
   r = 4 (d <= 30).
   - Free members satisfying the Williams half do exist: 48 distinct non-Psi f for r = 3 (d <= 40), 29 for
     r = 4 (d <= 30). Many are Dickson-type, e.g. f = 2 D_5(Y) - 1.
   - All of them fail the Korselt condition of the forced members.
 - Twisted members P_d(Y) = i^(-deg) Psi_d(iY) with 4 | d, at Y = L_M, M odd. Here z | Md is forced [PROVED:
   u = i phi^M gives iY = u + 1/u]. Korselt twisted families with 5 not | d, d <= 200, k <= 6: NONE
   [COMPUTED: src/fw/twist.gp; the untwisted control recovers {10,14,18} and {20,28,36}].
 - Fermat-forced Phi_d(2^K) families. A split member needs p - 1 | n + 1 (Williams' criterion for split p). With
   d | p - 1 and odd(d) | n - 1 this forces d = 2^m. Then Phi_d - 1 = X^(2^(m-1)) needs n(0) = -1, but n(0) = 1.
   So split members are impossible [PROVED], and the all-inert case is cyc_sym.out (0 pairs, d <= 160).
 - Other Selfridge D have no one-variable forcing. For Q* != +-1 the forced members are binary forms in
   (V_K, Q^(K/2)). A Korselt identity among distinct irreducible forms would force prod_(i!=j) F_i^(deg F_j) =
   F_j^(D - deg F_j) [PROVED], which unique factorisation rules out.
 STATUS: the exponential/polynomial-identity world (S1 and its mixed and twisted variants) is closed for every
 Selfridge D in the computed ranges. D = 5 is the only D with one-variable forcing.
 - Williams families whose special points are all real are impossible: by Lemma 6 every member takes the same
   value at each real special point, so all f_j - 1 = a_j G and f_j + 1 = a_j H, giving f_j equal [PROVED from Lemma 6].
   This covers every member class whose f -+ 1 split into Psi_e (cyclotomic special points), including the
   Dickson-type free members found above.
   Check [COMPUTED: src/rw/ratw.py]: all 17 irreducible members with special points in {-2..2} (complete for that
   class, degree <= 4), Williams families k <= 9: none.

### Literature tools (2026-10-24, alphaXiv; read by me: Grantham 1903.06820, Khashin 1807.07249)
 - Grantham Thm 4.3: a Frobenius psp for f is a Fermat psp to base f(0). For x^2 - x + 2 this gives base 2, which
   matches the D* = -7 reformulation in Lemma 14. Grantham offers $6.20 for a Frobenius psp for x^2 + 5x + 5 that is
   2 or 3 mod 5. He calls it "equally challenging" as PSW.
 - Khashin (Frobenius test, z = 1 + sqrt(c) or 2 + sqrt(c), c = first non-residue):
   - no Frobenius psp below 2^64;
   - every one has a prime factor > 3000;
   - Prop 3.4: no prime factor p with n/p < 2^21.
   Prop 3.4 uses the last-prime gcd (his Cor 2.10, our Lemma 14 consequence). He reports the gcd "decomposes into small
   primes", with 1-3 candidates per cofactor up to several million [READ, not re-checked]. This matches the
   BCZ-type smallness that sinks every small-cofactor construction.

### Strategic review (2026-10-24, session 2b)
 - Distance to result:
   - New rigorous results: Lemma 12(c) and Lemma 14.
   - Sub-routes closed by proof: real-special-point Williams families, split Fermat-forced members, multi-D forcing.
   - Closed by bounded computation: forced + free and twisted exponential families.
   - No new mechanism that could output n.
 - Mechanism census:
   - Pools: two-sided S-unit cap; the AGHS heuristic gap is a log factor in the exponent.
   - Identities: overdetermined; real special points dead; the cyclotomic coincidence gives one family, dead mod 5.
   - Few factors: Sum 1/p^2-type finite.
   - Forcing: one side only (Fermat via 2^K, Lucas via L_2K, and only for D = 5).
   Every mechanism needs, per prime, BOTH ord_p 2 and z(p) small or shared. That event is the DSO / S-unit
   coincidence, rare by Bugeaud-Corvaja-Zannier-type gcd bounds.
 - The missing ingredient, stated precisely:
   - EITHER a source of many primes whose two orders are simultaneously controlled. Every known source controls
     one order algebraically and leaves the other to smoothness or luck.
   - OR a subset-product algorithm for the structured groups (Z/A)^* x (Z/B)^* that works at density
     N ~ 2^(sqrt(B)/log B), below Wagner/Kuperberg's 2^(c sqrt B).
   Neither exists in the literature I have read. P1 (more CHM pool growth) cannot supply either. KILLED as pool
   enlargement; reason recorded in session 1's AGHS item.
 - Kill:
   - P1 (pool growth).
   - P2 (done: twisted none; Aurifeuillian Lucas factors carry 5 | z, so they die by Lemma 12(b)).
   - P3 (only useful with P1).
 - Keep or open (live; first experiment; kill condition):
   L1 Last-prime gcd search for BPSW (Lemma 14).
      Covers every BPSW psp having a prime q with n/q <= M0, q unbounded.
      First step: M0 = 2e4, D* in {5, -7, -11, 13, -15, 17, -19, -23}.
      Kill: gcds always small and no candidate. Then record the theorem and extend M0 only while the cost per
      decade stays under a day.
   L2 Theory: classify Korselt Psi-families (P4).
      Ingredients: the top common special point (Lemma 6 argument, Korselt half) and cyclotomic unit relations at the
      other special points.
      First step: prove that every member has Phi_d(1) = 1 (d not a prime power) and that 2 is a simple common
      special point, then bound the special-point fields.
      Kill: no progress after two derivation attempts. Then record the partial lemma.
   L3 Two-order source.
      Question: does any number-theoretic construction control ord_p 2 and z(p) at once without smoothness?
      Candidates to test:
       (a) primes dividing values of norm forms tying 2 and phi together, e.g. N(phi^k - 2^m) (inert p | it has
           z | 2k and ord | 4m). These are DSO primes, so first measure their number and size.
       (b) Mersenne primes M_q, q == 3 mod 4: inert for D = 5, ord = q, z(M_q) = 2^q exactly
           [PROVED: v2(z) = v2(p+1) for inert p == 3 mod 4]. They are DSO with B = 1. Products fail F1 (different
           2-adic classes 2^q); a single M_q needs class-(0,q) partners, i.e. p == -1 mod 2^q with ord_p 2 odd.
      First step: tabulate inert primes p | N(phi^k - 2^m) for k, m <= 200 and compare with the CHM pools.
      Kill: no prime > 2^40 and no density gain over CHM.

### L1 result: last-prime gcd search (src/lp/lastprime.gp) [COMPUTED]
 - Method:
   - For each Selfridge D, enumerate cofactors m = n/q: odd, squarefree, all primes in one 2-adic class,
     pairwise F2-compatible, prime to DQ.
   - Compute G = gcd(2^E1 - 1, U_E2) exactly, with q split iff (D/m) = -1.
   - Factor G and run bpsw(m q) on every prime q | G.
 - Controls [COMPUTED: src/lp/ctl_lp.gp]: both divisibility rules hold for every prime factor of all 72 squarefree
   spsp(2) and all 71 squarefree strong Lucas psp below 2e6 (0 violations).
 - THEOREM (computer-assisted): let n be a BPSW psp whose Selfridge D is one of 5, -7, -11, 13, -15, 17, -19, -23.
   Then n/q > 1e5 for every prime q | n. In particular n = pq needs min(p,q) > 1e5, and the larger factor is
   unbounded a priori.
   [COMPUTED: src/lp/lp_1e5.out. About 10.5k admissible cofactors per D; every gcd factored completely; 12k-20k
   candidate q per D; 0 hits; 31 s.] A run to 1e6 for all prime discriminants |D| <= 101 and -15 is in progress
   (src/lp/lp_1e6.out).
   (The run also includes D = -3, which Selfridge never selects; harmless.)
 - Gcd composition [COMPUTED: src/lp/lp_big.gp]:
   - The largest gcds (up to 1e127) are products of 20-30 SMALL primes (< 1e6). They occur when m - 1 is smooth,
     e.g. m = 89461 = 137*653 with E = 89460 gives 27 primes.
   - So for smooth E the gcd is just the DSO / CHM pool of primes whose orders divide E. Large-index DSO primes are
     rare, consistent with the heuristic count Sum_(d|E1, e|E2) 4/(de log de) = O(1).
   - The method gives a rigorous infinite-region exclusion, not a new pool.
 - Twisted families, d <= 600 (5 | d allowed), k <= 8: the only Korselt families are {40,56,72}, {80,112,144},
   {160,224,288} and {320,448,576}, i.e. 4*{10,14,18}*2^i. All have 5 | d, so z is divisible by 5 and Lemma 12(b)
   kills them for D = 5. [COMPUTED: src/fw/twist600.out]
 - Unfactored parts from the M0 = 1e6 run are post-processed with src/lp/lp_fix.gp (full factorisation of G(m)).
   D = -15: 4 such m. Each was a composite of 15-16 digits with prime factors just above the 1e7 trial bound.
   All resolved; 0 hits.

### Lemma 15 (PROVED; checked): the index budget of a 3-factor BPSW psp with all primes inert
Register rows: F1, F2 (via Lemma 14). Any Selfridge D*.
Let n = p1 p2 p3 be a BPSW psp, p1 < p2 < p3, all three inert for D*. Put O_i = odd(ord_{p_i} 2) * odd(z(p_i)) and
I_i = (p_i^2 - 1)/O_i. Note I_i = 2^(j1+j2) * [(p_i - 1)/ord] * [(p_i + 1)/z], i.e. the two order indices times the class.
Then
   I1 I2 I3 >= prod(p_i^2 - 1) / N,   N = 1/2 sum_cyc (p_i^2 - 1)(p_j - p_k)^2,
and hence I1 I2 I3 > (2/3)(p1^2 - 1)(1 - p2^-2)(1 - p3^-2). So some prime has index >= about 0.87 p1^(2/3).
For closely spaced primes the bound is much larger, about p^4/Delta^2.
Proof.
 - Lemma 14 (inert case) gives p_j p_k - 1 = c_i O_i with c_i a positive integer.
 - Hence prod(c_i / I_i) = R := prod_cyc(p_j p_k - 1) / prod(p_i^2 - 1).
 - Identity: prod_cyc(p_j p_k - 1) - prod(p_i^2 - 1) = 1/2 sum_cyc (p_i^2 - 1)(p_j - p_k)^2 = N > 0.
   To see it, expand both sides in s1, s2, s3. The difference is (s2^2 - 3 s1 s3) - (s1^2 - 3 s2), and each bracket is a
   sum of squares.
 - So R > 1, and R prod I_i = prod c_i is an integer. Then prod c_i >= prod I_i + 1, so R - 1 >= 1/prod I_i.
 - Finally N < 1/2 p3^2 (p1^2 + 2 p2^2) < (3/2) p2^2 p3^2.  QED
Scope:
 - It applies only to the all-inert 3-factor case.
 - With one or more split primes, Lemma 14 gives O_i | m_i^2 - 1, the degree count no longer balances, and no
   constraint follows.
 - For k = 2 the same count gives I_p I_q >= (p+1)(q+1) (known heuristic, now exact).
 - For k >= 4 there is no constraint.
Checks [COMPUTED: src/lem/l15.gp]:
 - The identity was checked symbolically.
 - N < (3/2) p2^2 p3^2 held on 1e5 random triples.
 - The integrality argument was tested non-vacuously on all 551,300 odd triples <= 301, taking O_i = gcd(x_i^2 - 1,
   x_j x_k - 1): 0 violations, and 35 cases within a factor 2 of the bound, so the bound is nearly sharp.
 - Live test: no triple of inert primes (D = 5) with n <= 2e7 has all three c_i integral, so it is vacuous there.
Use:
 - A 3-factor all-inert BPSW psp needs a large combined index somewhere, i.e. a DSO-type prime of size >= p1^(2/3).
 - This complements L1 (which covers n/q <= 1e6).

### Remark (HEURISTIC + Chebotarev): Kummer-forced indices cannot manufacture DSO primes
 - Forcing I_F >= k and I_L >= k' for every prime of a family means the primes split completely in
   L = Q(zeta_k, 2^(1/k), sqrt D, zeta_k', beta^(1/k')), with [L:Q] about k phi(k) k' phi(k').
 - The least such prime is expected to be >= [L:Q] (Linnik-type). So O_p = (p^2 - 1)/I >= p^2/(k k') >= p^(3/2):
   at best half-way to DSO, and the remaining order is not smooth.
 - The rank-2 group <2, beta> in K* (beta has norm 1 and is not torsion; 2 is not a unit) makes "both orders small at
   p" a codimension-2 condition in every extension field, so no field change removes it.
 - The only known way around it is exponent arithmetic (smooth orders shared across primes), i.e. the S-unit pools.

### Attempt D2 (2026-10-24): DSO census at highly composite exponents [COMPUTED: src/dso2/]
Register rows: none chosen (a census). Question (missing ingredient (ii)): do large-index primes with BOTH orders
small exist when the orders are allowed to be highly composite? Attempt D only measured ord, z <= 1500.
 - Method: G(E) = gcd(2^E - 1, U_E(1,Q)) contains EVERY prime with ord_p 2 | E and z(p) | E.
   Take E = lcm(1..k), k <= 20 (E <= 2.3e8), factor G completely (src/dso2/hcgcd.gp, cof.gp), then test BPSW
   usability: p inert, or p split with gcd(odd ord, odd z) = 1 (F2).
 - D = 5, k = 20: 240 primes, the largest 232792561 = E + 1.
   - Every prime above ~2e4 is SPLIT with ord_p 2 and z(p) sharing almost all odd primes, e.g.
     ord = z = 2^4 3 7 11 13 17 for p = 22870849. These "diagonal" primes have p - 1 smooth and both indices small.
     They are useless for BPSW: a shared odd l would need l | n - 1 and l | n + 1.
   - The gcd cofactors above the 1e7 trial bound (30 and 105 digits) were products of 8-9 digit diagonal primes.
   - BPSW-usable: 12 primes, all inert, max 1597.
 - Six Selfridge D (5, -7, -11, 13, -15, 17), k = 16 and 20 (src/dso2/usableD.gp):
   235-257 primes each, 13-23 usable, largest usable <= 51481.
 - Conclusion [COMPUTED, bounded]: highly composite exponents produce many DSO primes, but only diagonal split ones.
   No large usable DSO prime exists for orders dividing lcm(1..20), for any of these D.
   This is the codimension-2 obstruction seen directly: one-sided smoothness (p - 1) gives diagonal primes;
   usability needs two-sided smoothness or two independent index coincidences.
   Missing ingredient (ii) is not available from this source.

### Strategic review (2026-10-24, session 2c)
 - New since 2b:
   - L1 theorem to 1e5 (all 8 small D) and, in progress, 1e6 for 26 D.
   - Lemma 15 (3-factor budget).
   - Kummer remark.
   - Attempt D2 census.
   - fw2 r2 d <= 150 (running).
 - Conservation principle (HEURISTIC, with rigorous pieces): every way of making one constraint automatic costs the
   same probability elsewhere. Instances checked here:
   - Smooth n - 1 forces z-roughness. With n == +-1 mod lcm(1..y), all D < y have (D/n) = 1 and the other side
     must avoid every l <= y. [PROVED, elementary]
   - Shared large primes R | p - 1 cost 1/R per pool prime.
   - Proth/Riesel shapes control one side only.
   - Q* == 1 mod p gives z(p) = 3, but only for p | D* + 3 (small). [PROVED]
   - Two-base exponential families f(2^t, phi^2t) cannot force both orders: f would have to divide an X-only and a
     Y-only polynomial. [PROVED]
   - Kummer forcing gives O_p >= p^(3/2). [heuristic + Chebotarev]
   - Highly composite exponents give only diagonal split primes. [COMPUTED]
 - Distance: unchanged. No mechanism outputs n.
   - Few-factor regions beyond L1 have heuristic expected yield <= 1e-6.
   - The pool route needs >= 2^100 work at every scale [HEURISTIC; session 1 AGHS item].
 - Live CONSTRUCTIVE directions (G4), each with first step and kill:
   C1 L1 extension: finish 1e6 x 26 D, then the next decade only for D in {5, -7, -11, 13, -15}.
      Kill: a decade costs > 1 day.
   C2 Structured spsp(2) lottery with exact Lucas side: composite Mersenne M_q (q <= 1e4), Wagstaff (2^q+1)/3,
      Fermat F_k (k <= 16), composite divisors of L_(2^k).
      - Each family is automatically strong in one test (spsp(2), or strong Lucas for L_(2^k)).
      - For q == 3 mod 4, BPSW on M_q is a Lucas-Lehmer run with seed 3: it passes iff all prime factors share
        z = 2^(r+1).
      - Kill: done once (minutes). It is a lottery, recorded as such.
   C3 Exact exhaustion of small smooth pools.
      - Rigorously decide whether any BPSW psp is a product of primes with p^2 - 1 y-smooth, for y = 41 then 100,
        n unbounded.
      - First step: compatibility cliques per class and per side-assignment, then exact subset enumeration (Gray code).
      - Kill: cliques > 45 primes at y = 41.
   C4 (theory feeding constructions) Exceptional cyclotomic units.
      - The Fermat-exponential pair condition Res(Phi_a + 1, Phi_b) = +-1 says Phi_a(zeta_b) + 1 is a unit, i.e.
        an exceptional-unit pair (u, u - 1) in Q(zeta_b).
      - Classifying these for all a, b would close or open the cyclotomic half-family world for every d, not just
        d <= 160.
      - First step: literature check (exceptional units in cyclotomic fields), then derive the classification of
        Phi_a(zeta_b) + 1 being a unit.

### C2 result: structured spsp(2) lottery [COMPUTED: src/lot/lot.gp -> lot.out]
Exact BPSW on the following families, all 0 passes:
 - all composite M_q, q prime <= 1e4 (1207 numbers);
 - all composite Wagstaff (2^q+1)/3, q <= 1e4 (1204);
 - F_5..F_16;
 - all products of >= 2 known prime factors of L_(2^k), k <= 9.
This is a lottery (expected yield about 0), recorded as such; C2 is closed.
fw2 r = 2, d <= 150 (forced Psi pair + one free Williams member):
 - The process died (no final line, last write 04:30) after the pairs starting at d = 118, with 0 hits.
 - Recorded as UNFINISHED and partial (G3). Not re-run: the sub-route is already closed in range by r <= 4 runs.

### Structured subset products (2026-10-24): re-examining missing ingredient (i)
Register rows: C1, C2, C5 (pool setting).
 - Observation [PROVED, elementary]:
   - In G = (Z/A)^* x (Z/B)^*, A, B odd, the target (1, -1) has trivial image in every odd l-primary component, since
     -1 has order 2.
   - So the odd part of the subset-product problem is a ZERO-SUM problem: find nonempty S with sum_{p in S} v_p = 0 in
     G_odd = sum_l G_l, together with an affine F_2-type condition on the 2-part (Jacobi symbols, the C_2 parts, the
     classes).
   - Also (Z/r)^* = C_(r-1) splits into l-parts for l | r - 1. So G_odd is a product of many small elementary-abelian
     l-groups (l small, multiplicity about pi(y)/(l-1)) and a "cyclic" remainder (large l, multiplicity ~1).
 - Measured structure of the best annealed pools [COMPUTED: src/grpstruct.py on data/sn_150.txt]:
   z = 150, class (2,0), N = 13 primes, log2|G| = 91.4 bits:
   - 2-part 28.0 bits (17 components, up to C_16);
   - odd l <= 13: 52.1 bits (3: 10 components incl. Z/9; 5, 7, 13: 3 each; 11: 1);
   - odd l > 13: only 10.3 bits (l = 23, 53).
   The annealer already selects vocabularies whose r - 1 is smooth, so the group is mostly structured.
 - Literature [READ by me: Kothari-Metger-O'Donnell-Shutty-Wu, arXiv 2609.40321; it summarises Kothari-O'Donnell-Wu,
   STOC 2026]:
   - For random vectors in F_q^n, a nonempty zero-sum subset is found in poly time with m >= q^(q^2) n^(1 + floor((q+3)/4))
     vectors (average case, rigorous). Worst case needs m >= q^(q^2) (n+1)^(q-1).
   - For q = 3: m about n^2/3 (rigorous), about 0.1022 n^2 (informal).
   - Wagner-type methods cost 2^O(n / log(m/n)). The information threshold is m about n log2 q.
   So per field the threshold is POLYNOMIAL, not 2^(c sqrt B).
   This is the first algorithmic family that could in principle beat the session-1 AGHS/Kuperberg bound for BPSW pools.
 - Why it does not rescue the route at laptop scale [HEURISTIC, with the numbers above]:
   (a) Algorithms matter only at density N/B > 1. The measured best density is 0.11-0.14 for every B <= 272 at
       z = 100 and z = 150, and N_best(B) is so far LINEAR in B:
       (59, 8), (91, 13), (126, 16), (212, 26), (272, 30).
   (b) For a product of different fields, the only general method I see is staging, one field at a time. Its element
       cost MULTIPLIES across fields: for the z = 150 group, about 1e11 elements versus the information threshold of 91.
       A joint algorithm for products of small fields, or over Z/Q with Q = prod l, is not in the papers read. That is
       the precise open algorithmic question.
   (c) Even with a joint algorithm needing N about B^2, the S-unit extrapolation (raw pair count
       ~ exp(1.89 sqrt(pi(y))), fitted at y = 41, 100, 150 to within 1%) puts the crossover at y about 700-1000.
       There are about 1e10 smooth pairs there, beyond laptop enumeration.
 - Decision:
   - Record (b) as the sharpened form of missing ingredient (i): a joint zero-sum algorithm for random vectors in
     prod_l F_l^(n_l), small l, with polynomial sample threshold.
   - Measure the frontier trend at z = 200 (P1 revived for this reason only; CHM running) to see whether N_best(B)
     bends upward before B ~ 500.

### Lemma 16 (PROVED; checked): composite L_(2^k) are strong Lucas pseudoprimes for Selfridge's parameters
For k >= 2, n = L_(2^k) = V_(2^k)(1,-1) passes the strong Lucas test with Selfridge D = 5.
Proof.
 - Lucas numbers mod 5 have period 4 (2, 1, 3, 4). So n == L_0 = 2 mod 5 and (5/n) = (n/5) = -1, i.e. D* = 5.
 - L_(2^(k+1)) + 1 = L_(2^k)^2 - 1 = (L_(2^k) - 1)(L_(2^k) + 1). By induction from L_2 + 1 = 4, 2^(k+1) | n + 1.
   So n + 1 = d 2^s with s >= k+1 and d odd.
 - V_(2^k) | V_(2^k d) for odd d, hence V_(d 2^k) == 0 mod n, with k < s.  QED
Every composite L_(2^k) (k >= 5) is therefore BPSW iff it is a strong pseudoprime to base 2. This is the exact mirror
of composite Fermat numbers, which are spsp(2) and BPSW iff strong-Lucas.
Check [COMPUTED: src/lot/lucas2k.gp]: k = 5..15 all give D = 5 and strong Lucas = 1, as proved; all fail spsp(2).
k = 16..20 is running (lottery).
Both mirror families are one-side-automatic. Their products mix incompatible 2-adic classes (Fermat factors
(m+1, *) against Lucas factors (0, k+1)).

## Session 3 (2026-10-24, continued)
### P1 result: the frontier at y = 200 and its growth curve [COMPUTED: src/frontier_np.py, frontier_np2.py, frontier3.py]
 - CHM closure at y = 200 finished: 345,275 smooth pairs (21 min; data/sn_200.txt). Inert rows: 13,525.
   Largest classes: (2,0) 3233, (1,2) 2119, (0,3) 1439, (0,4) 995, (3,0) 749, (0,5) 701.
 - Tools. frontier_np.py is a vectorised copy of frontier.py (same annealer, same objective). It reproduces the
   z = 100 frontier exactly, in 2.7 s instead of minutes.
   frontier3.py adds per-prime exponent caps (side A/B/X plus a cap e_q), an optional split-prime pool, and
   restriction of one pair file to any y' < y.
 - Frontier at z = 200, class (2,0) [the best class], 4 seeds x 200k iterations:
   (B, N) = (397.6, 58), (367.5, 57), (335.1, 52), (43.0, 6). Density N/B <= 0.155 everywhere.
   - 50x more annealing than the z = 150 run moved N from 55 to 58 at B ~ 390, so the annealer is not the limit.
   - Exponent caps change nothing measurable (y = 100: N = 14 at B = 107, ratio 0.131).
 - Growth curve (frontier restricted to y' from the one z = 200 file; best ratio over lambda):

   | y' | 60 | 80 | 100 | 120 | 140 | 160 | 180 | 200 |
   |---|---|---|---|---|---|---|---|---|
   | class (2,0) size | 49 | 112 | 185 | 419 | 722 | 1081 | 1814 | 3233 |
   | best N/B | 0.149 | 0.142 | 0.143 | 0.141 | 0.144 | 0.143 | 0.151 | 0.155 |
   | N at B ~ 110-150 | 8 (B 65) | 14 (B 120) | 15 (B 114) | 20 (B 152) | 18 (B 126) | 16 (B 112) | 14 (B 100) | - |

   The class grows 66x from y' = 60 to 200; the best density grows from 0.149 to 0.155.
   At FIXED B ~ 120, N_best stays 14-20 for every y' from 80 to 200.
 - Split primes [COMPUTED]: of 39,875 split primes p = 2x+1 from the y = 200 pairs, exactly ONE has coprime odd
   parts of ord_p 2 and z(p) (and 5 does not divide either): p = 974849, a factor of F_11 (ord = 2^12,
   odd(z) = 7*17). All the others are diagonal. Split primes contribute nothing to pools.
 - Conclusion [COMPUTED, bounded; the extrapolation is HEURISTIC]: within S-unit (smooth-neighbour) pools the
   achievable density is N/B ~ 0.15, flat in y and in B up to B ~ 400.
   - At small B (where generic subset-sum algorithms would be affordable, e.g. BCJ ~ 2^(0.29 B) for B <= 150) the
     expected number of solutions is 2^(N - B) ~ 2^(-100).
   - At large B, even N/B >= 1 would leave a density-1 subset-sum problem, about 2^(0.29 B) work.
   So the pool route fails at EVERY B reachable here: small B gives too few primes, large B gives an intractable search.
   That closes P1 for this project's scale. The missing ingredient stays (i) a structured-group algorithm or
   (ii) a non-S-unit source of two-sided primes.

### Literature (2026-10-24, alphaXiv; read by me) - constructions and tools only
 - Shallue-Webster, arXiv 1806.08697 ("Fast tabulation of challenge pseudoprimes"):
   - No PSW challenge pseudoprimes (psp(2), (5/n) = -1, n | F_(n+1)) with 2 or 3 prime factors below 2^80.
     This covers the D* = 5, t <= 3 part of BPSW below 2^80 (BPSW with D* = 5 implies PSW).
   - Algorithm: pre-product k; the last prime p divides gcd(2^(k-1) - 1, U_(k-eps(k))) when k is small (our Lemma
     14 / L1), else a sieve in the progression p == k^-1 mod L, -k^-1 mod W.
   - 7 "admissible" split primes below 2^40 (gcd(ord, z) <= 2), e.g. 4278255361 = Phi_80(2) with ord 80.
   - 2047 = M_11 is a (2, 23, 131)-challenge pseudoprime. D = 23^2 - 4*131 = 5, and (P, Q) is chosen so that
     gamma has 2-power order mod 23 and mod 89. With per-n parameters the Lucas side is trivial; Selfridge's
     rigid choice is what makes BPSW hard.
   - Rewards: $20 for a PSW challenge psp with an even number of prime factors, $20 with exactly three, $6 divisible by 3.
 - Baillie-Fiori-Wagstaff, arXiv 2006.14425:
   - Theorem 1: every spsp(2) n == 3 mod 4 is a strong lpsp(2^k, 2^(2k-1)). The reason: D = -4^k, so gamma = i is
     torsion. Selfridge never selects such a degenerate D: for P = 1, gamma is torsion only for D = -3.
   - Pomerance's heuristic construction (Section 8): p in [T, T^k], p == 3 mod 8, (5/p) = -1. (p-1)/2 is
     squarefree over primes q < T with q == 1 mod 4; (p+1)/4 is squarefree over primes q < T with q == 3 mod 4.
     Then n == 1 mod Q1 and == -1 mod Q3 with an odd number of factors.
     - The partition is FIXED by q mod 4, and the class is automatically (1,2): v2(z) = v2(p+1) = 2 for inert
       p == 3 mod 4.
     - The authors note a search with this construction "would be very slow".

### DSO-closed families: exhaustive subset check [COMPUTED: src/dso2/subsets.gp -> subsets.out]
Every BPSW psp with Selfridge D* = D whose prime factors ALL have ord_p 2 | E and z_D(p) | E (E = lcm(1..20))
is a product of the usable primes of G(E) = gcd(2^E - 1, U_E). All subsets were enumerated after the Lemma 1
filters (one 2-adic class; odd(ord) coprime to odd(z) over the set), and each survivor got the exact bpsw().
 - D = 5, -7, -11, 13, -15, 17, -19, -23: 13-23 usable primes each, only 5-23 class-consistent subsets, 0 hits.
 - THEOREM [COMPUTED]: no BPSW psp with D* in this list has all prime factors with both orders dividing lcm(1..20).
   n is unbounded within the family (every subset was tested).

### Lemma 17 (PROVED, elementary): split primes of a BPSW psp are thin; a 2-factor BPSW psp has exactly one
Register rows: F1 (classes), F2 (coprime sides). Let n be a BPSW psp with Selfridge D*, and p | n split ((D*/p) = 1).
 (a) odd(ord_p 2) * odd(z(p)) | odd(p - 1), and min(v2 ord_p 2, v2 z(p)) <= 1.
     Proof: both orders divide p - 1, their odd parts are coprime (one divides n - 1, the other n + 1), and
     j1 <= v2(n-1), j2 <= v2(n+1), where one of v2(n -+ 1) equals 1.
     Consequence: I_F(p) * I_L(p) >= (p - 1) / 2^(j1+j2) * 2^(v2(p-1)) >= (p-1)/4, so one index is >= sqrt(p)/2.
     Every split prime of a BPSW psp is a one-sided DSO prime. Call such p ADMISSIBLE.
 (b) n has an odd number of inert primes (since (D*/n) = -1, and no p | D*Q*). If n = pq, exactly one of p, q
     is split.
 (c) Last-prime form. If q is the ONLY inert prime of n, put m = n/q (the split part). Then
     odd(ord_q 2) * odd(z_q) | O := odd(m - 1). So q | gcd(2^(2^j1 A) - 1, U_(2^j2 B)) for a unitary divisor A of O
     (B = O/A); take A generated by the primes of odd(ord_q 2). This is Lemma 14 at q.
     Each of the 2^omega(O) gcds can be computed modulo the smaller number, of about 2^j sqrt(O) bits.
 Census [COMPUTED: src/us/usplit.c, a segmented sieve with a largest-prime-factor filter; 7.6 s to 1e8 for 6 D]:
 admissible split primes below 1e8:
   - D = 5: 1 (61681: ord 40, z 1542). This matches Shallue-Webster, whose next is 363101449.
   - D = -7: 5. D = -11: 3. D = 13: 8. D = -15: 4. D = 17: 4.
   Examples: 6700417 | F_5 (D = 13, ord 64); 2796203 = (2^23 + 1)/3 (D = 17, ord 46); 15467 (D = -11, z = 19).
 Search [COMPUTED: src/us/gq.gp (lastq), drive.gp; control ctl.gp]:
   - Every class-consistent, cross-coprime product m of admissible split primes below 1e8 (m <= 1e16; 1-11 cofactors
     per D) has 0-5 candidate q. 0 BPSW hits.
   - Control: for m = 5581 (D = -11) the solver returns q = 61 with O_61 = 15*31 | 5580, as it must.
 THEOREM B [COMPUTED: src/us/us_1e10.txt (scan to 1e10, 455M primes, 22 min), run10.gp, pairs10.gp, big29.gp]:
 for D* in {5, -7, -11, 13, -15, 17, -19, -23, 29, -31} there is no BPSW psp with Selfridge D* that has exactly one
 inert prime factor and all its split prime factors < 1e10. The size of n and of the inert prime is unbounded.
   - Admissible split primes < 1e10 per D: 5, 5, 4, 11, 6, 4, 2, 4, 8, 6.
   - Every class-consistent, cross-coprime product of them was tested: 62 split parts, 20 candidate inert q,
     0 hits. Only one split part exceeds 1e16: D = 29, m = 32665301 * 1244942269 = 4.07e16, with
     O = 12791 * 54419 * 912853 (0 candidates, 1 s).
   - Independent check: the D = 5 list is EXACTLY Shallue-Webster's admissible primes below 1e10 (61681,
     363101449, 4278255361, 4562284561, 4582537681), found by an independent program.
 COROLLARY: every 2-factor BPSW psp n = pq with D* in this list has its split prime > 1e10.
   - D = 5 beyond 1e10: SW's two further primes below 2^40, 26509131221 and 422013019339 (orders re-verified),
     give 0 candidates (sw2.gp).
   - So the corollary holds with 2^40 for D* = 5 IF SW's list is complete; my own D = 5 scan to 2^40 is running
     (usplit3, Montgomery arithmetic, about 4.5 h).

### Missing ingredient (i) made precise: it contains modular subset sum [PROVED reduction; consequence HEURISTIC]
 - For distinct primes l_1..l_r, the sum of the cyclic groups Z/l_i is cyclic of order M = prod l_i (CRT).
   So any algorithm finding zero-sum (or target-sum) subsets of random elements of
   G = (Z/A)^* x (Z/B)^* in poly(log|G|) time would, restricted to the large-prime part of G_odd (multiplicity 1:
   q == 1 mod l for one q only), solve random modular subset sum mod squarefree smooth M at density ~1.
   No such algorithm is known (lattice methods need density < 0.94 and a planted solution; BCJ-type methods
   cost 2^(0.29 n)).
 - The tractable structure (F_2 linear algebra; F_3^n via Kothari-O'Donnell-Wu) needs SMALL EXPONENT.
   Forcing small exponent (moduli q with q - 1 | E0 fixed, e.g. E0 = 5040) shrinks the S-unit pool. Staging the
   three fields of E0 = 5040 needs about n_3^2 n_5^2 n_7^2 ~ 1e7 samples against pools of ~50.
 - Status: (i) is a subset-sum breakthrough, not a missing trick. Redirected; no further algorithm work planned
   unless a pool with N/B >> 1 at B <= 150 appears. That needs a new prime source, i.e. (ii).

### Lemma 18 (PROVED; checked): the 2-adic classes are governed by (2/p) and (Q*/p)
Register rows: F1 (classes). Any D*, Q* = (1 - D*)/4, p not dividing 2D*Q*, eps = (D*/p).
 (a) v2(ord_p 2) = v2(p - 1) iff (2/p) = -1; otherwise it is smaller.  [standard]
 (b) v2(z(p)) = v2(p - eps) iff (Q*/p) = -1; otherwise it is smaller.
Proof of (b) for inert p: let g generate F_(p^2)^*, h = g^(p-1) (it generates T_p, order p + 1), and write
alpha = g^a. Then alpha^(p+1) = N(alpha) = Q* = g^((p+1)c), so a == c mod (p - 1), say a = c + (p-1)t.
Then gamma = alpha/beta = alpha^2/Q* = g^(2a - (p+1)c) = h^(2t - c). So gamma is a square in T_p iff c is even,
i.e. iff Q* is a square in F_p^*. Split p: gamma = alpha^2/Q* in F_p^*, the same conclusion.  QED
 [COMPUTED: src/lem/l18.gp, all primes p < 2e5 for 10 Selfridge D, about 18,000 inert and 18,000 split each:
  0 violations of (a) or (b).]
Consequences.
 - D* = 5 (Q* = -1) has three families of inert primes:
   - p == 3 mod 8: class (1,2);
   - p == 7 mod 8: class (0, v2(p+1)) with v2(p+1) >= 3;
   - p == 1 mod 4: class (j1, 0), i.e. z(p) is ODD.
   Since all primes of n share one class:
   - class (1,2): every inert prime is == 3 mod 8 (Pomerance's family);
   - class (0,j), j >= 3: every inert prime is == 2^j - 1 mod 2^(j+1);
   - class (j1,0): every inert prime is == 1 mod 4.
   This is the class list the frontier found: (2,0), (1,2), (0,3), (0,4), (3,0), (0,5).
 - D* = -15 (Q* = 4, a square): z(p) never has the full 2-part. D* = -7 (Q* = 2): the Fermat and Lucas sides are
   full or reduced together, so inert p == 3, 5 mod 8 give classes (1,2), (2,1).

### Strategic review (2026-10-24, session 3)
 - New since 2c:
   - P1 closed at this scale: pool density about 0.15, flat in y and B.
   - Lemma 17 (split primes are admissible and thin) with Theorem B (exactly one inert prime, split primes
     < 1e10, ten D).
   - Lemma 18 (2-adic classes from (2/p) and (Q*/p)).
   - DSO-closed theorem.
   - Ingredient (i) reduces to modular subset sum.
 - Where a BPSW psp can still be (for |D*| <= 31, combining L1, Lemma 17, Theorem B and the DSO-closed theorem):
   - n has >= 3 inert prime factors, or a split prime > 1e10;
   - n/q > 1e5 for every prime q | n (1e6 once lp_1e6 finishes);
   - not every prime factor has both orders dividing lcm(1..20).
   Every split prime is one of a thin list (5-11 per D below 1e10).
 - Distance to an explicit n: unchanged. Every remaining region has heuristic expected yield <= 1e-5:
   - 3-inert with q1 q2 > 1e6: about sum 2/Lambda(m);
   - pools: 2^(N - B) with N/B ~ 0.15.
 - Theory of success: none at laptop scale. The two missing ingredients are pinned down:
   (i) subset sum (a known hard problem);
   (ii) a non-S-unit source of inert primes with both orders controlled. For inert p, 2 lies in F_p^* and gamma in
       T_p (they meet only in +-1), so control needs two coincidences. Only exponent sharing (smoothness) or
       Kummer splitting provides them, and both are measured or estimated.
 - Live directions after this review (>= 3):
   E1 Exclusion theorems (running; cheap):
      - Theorem B to 1e11 for nine D (running) and 2^40 for D = 5 (running);
      - L1 at 1e6 (running).
      Kill: a decade costs > 1 day.
   E2 Theory: real polynomial Williams families. A family would give infinitely many BPSW psp (all values == 3
      mod 8, inert). The open step is positivity at non-real special points. First step: express the
      special-point conditions as |n| = 1 on the zero set of n^2 - 1 (all special points satisfy n(theta) = +-1,
      and lcm(f_j^2 - 1) | n^2 - 1). Then look for a Jensen / argument-principle contradiction.
      Kill: two derivation attempts without a new inequality.
   E3 Three-inert exact search beyond L1 with the cheap unitary-divisor gcd (cofactors q1 q2 in (1e6, 1e9]).
      Only if E1 finishes and a C implementation of lastq is written (about 2e7 pairs per D).
      Kill: cost > 1 day per D.

## Session 4 (2026-10-24, continued)
### C3: exact exhaustion of two-sided-smooth pools [COMPUTED: src/c3/exact.py, export.py, enum.c]
Register rows: none chosen. The statement is about every n whose prime factors all lie in the pool.
 - Pool P_y = { p prime : p^2 - 1 is y-smooth } (p = 2x + 1 with x, x+1 y-smooth). This is exactly the pool of the
   Pomerance / Chen-Greene construction. For a Selfridge D, every BPSW psp with D* = D and all prime factors in
   P_y is a product of
   - inert p, or admissible split p (Lemma 17a),
   - in ONE 2-adic class (F1),
   - with no odd prime on both sides (F2): an independent set of the side-conflict graph,
   - an odd number of inert primes, none dividing D*Q*, with Lemma 12(b) applied.
   The enumerator lists EVERY such independent set (n unbounded) and tests the exact Lemma 1 congruences
   (prime-power residues; the 2-adic part mod 2^(max j + 1)). Hits are re-verified with the verbatim bpsw.
 - Completeness of the pair lists. Our CHM closure (src/chm.c) is NOT complete
   [READ by me: Conrey-Holmstrom-McLaughlin 1212.5161; Buzek et al. 2211.04315].
   - y = 41: 868 pairs. It misses x = 63927525375 = 3^3 5^3 7^7 23 (x + 1 = 2^13 11^4 13 41), Lehmer's largest.
     CHM report the same omission. With it the list has 869 pairs = Lehmer's complete count (1964).
     File: data/sn_41_complete.txt.
   - y = 100: sn_100 has 13,340 pairs. The union with the 100-smooth members of sn_200 plus
     9591468737351909375 (79-smooth, above our 2^63 cap; CHM's only miss) has 13,374 pairs. That is exactly
     Luca-Najman's complete count: 13,325 plus 49 in their erratum.
     Every element was re-verified smooth, so our set equals theirs; completeness is CITED.
     File: data/sn_100_complete.txt.
   - Consequence for earlier work: the P1 frontier used CHM lists that miss about 0.25% of pairs. The
     qualitative conclusions do not change.
 - Controls.
   - Fermat-only mode (Lucas data dropped, base-2 strong test at the end): y = 41, at most 3 primes. 31 hits, all
     spsp(2): 2047, 3277, 4033, 4681, 15841, 29341, ... Brute force over all spsp(2) < 3e6 with at most 3 prime
     factors in P_41 \ {5} gives 16 numbers, and the enumerator finds exactly these 16. C and Python agree.
   - Lucas-only mode: y = 41 finds 5777, 10877, 5459, 24569, ... (strong Lucas psp). Every congruence hit is a
     strong Lucas psp for the tested D; the ones that fail the full test have Selfridge D = 5, an earlier D.
   - C and Python agree on every subset count at y = 41.
 - y = 41 (complete list), D* in {5, -7, -11, 13, -15, 17, -19, -23, 29, -31}:
   - pools 82-162 primes, at most 1 admissible split prime per D, classes of at most 37;
   - 55-876 consistent odd subsets per D, 0 congruence hits.
   The side-conflict graph is very dense.

### Theorem C (PROVED, computer-assisted; completeness CITED): no BPSW psp has all prime factors in P_100
Register rows: none chosen. The statement covers every n, with no size bound.
 - Statement. Let P_100 = { p prime : p^2 - 1 is 100-smooth }. Then no BPSW pseudoprime with Selfridge
   D* in {5, -7, -11, 13, -15, 17, -19, -23, 29, -31} has all of its prime factors in P_100.
 - Ingredients:
   - Luca-Najman's complete list of 13,374 pairs (x, x+1), both 100-smooth (CITED; our file
     data/sn_100_complete.txt reproduces the count exactly). This gives 3043 primes p = 2x+1.
   - Lemma 1 (F1, F2), Lemma 12(b), Lemma 17 (admissible split primes, odd number of inert primes).
   - Squarefree: p^2 | n with n spsp(2) forces 2^(p-1) == 1 mod p^2. No prime of P_100 satisfies this
     (checked directly).
   - Exhaustive enumeration with src/c3/enum.c over the exports src/c3/ex100.txt.
 - Result (src/c3/enum100.out): 384 classes over the 10 D; 3,336,445,924 independent sets; 1,668,209,907
   parity-odd congruence checks; 0 hits; complete (no cap).
   - Pools per D: 804 (D=5) to 1571 (D=-19) primes, in 29-48 classes. Largest class: 399 (D=-15, class (2,0)).
   - Independent sets per D: 2.9e6 (D=29) to 1.98e9 (D=-15).
   - Largest independent set: 24. The side-conflict graph keeps every consistent subset to at most 24 primes,
     out of pools of about 1500.
 - What it does NOT say: nothing about n with even one prime factor outside P_100. y = 150 and y = 200 (CHM
   lists, not proven complete) are bounded evidence only (running).

### Pomerance pool by direct generation: saturation and density [COMPUTED: src/c3/pomgen.c -> pomgen.out]
 - Pool: p = 4B - 1, B a product of primes == 3 mod 4 (<= T) and A = 2B - 1 a product of primes == 1 mod 4
   (<= T, not 5). For fixed T this is an S-unit equation 2B - A = 1, so the pool is FINITE. Measured
   saturation: no new primes from 1e12 to 1e18.2 at T = 300, or from 1e15 to 1e17 at T = 500.
 - N / log2|G|, where |G| is the group generated by the moduli actually used:

   | T   | 100   | 200   | 300   | 500   |
   |-----|-------|-------|-------|-------|
   | N/B | 0.064 | 0.079 | 0.102 | 0.161 |

   At T = 500: N = 156 primes, log2|G| = 969.
 - Growth is accelerating (a local exponent of about 0.9 between T = 300 and 500). Read naively, the crossover
   N = log2|G| would come at T of a few thousand, where the pool primes are far beyond 2^64.

### Lemma 19 (PROVED): complete classification of 3-member Williams polynomial families over C
Register rows: C4, C3F, C3L. Supersedes the computer-assisted Lemma 4 with a conceptual proof valid over C.
 - Statement. Let f, g, h in C[t] be nonconstant, n = fgh, with f -+ 1 | n -+ 1 (and the same for g, h).
   Then {f, g, h} = {f, zeta f, zeta^2 f} with zeta^3 = 1. Over R (hence over Q), f = g = h.
 - Proof.
   1. f -+ 1 | n -+ 1 means gh == 1 mod (f^2 - 1), so f^2 - 1 | gh - 1 and gh - 1 != 0 (as k >= 2 and the
      polynomials are nonconstant). So 2 deg f <= deg g + deg h, and cyclically.
   2. Summing, all three inequalities are equalities, so deg f = deg g = deg h = d. Hence gh - 1 = a(f^2 - 1)
      for a constant a. Multiplying by f: n = a f^3 + (1 - a) f, i.e. n = P_a(f). Likewise n = P_b(g) = P_c(h).
   3. Uniqueness of the right component of a given degree (approximate cube root: the polynomial part of
      n^(1/3) determines f up to an affine map) gives g = u f + v and h = u' f + v'.
   4. Compare P_b(u f + v) = P_a(f) as polynomials in f. The f^2 coefficient gives v = 0. Then
      gh = u u' f^2 = a f^2 + 1 - a gives a = 1 and u u' = 1. Similarly fh = b g^2 + 1 - b gives b = 1,
      u' = u^2; and fg = c h^2 + 1 - c gives c = 1, u = u'^2. So u^3 = 1.  QED
 - Coprime normal form for any odd k (PROVED, same idea). If n -+ 1 = prod (f_j -+ 1) (the pairwise-coprime case of
   Lemma 8(b)), then Pi(X) := prod_j (X - f_j) satisfies Pi(1) - Pi(0) = 1 and Pi(-1) - Pi(0) = -1. Hence
     prod_j (X - f_j) = X (X^2 - 1) R(X) + X - n,   with R in C[t][X] monic of degree k - 3.
   - k = 3: R = 1, so prod (X - f_j) = X^3 - n and the f_j are zeta^j n^(1/3).
   - k = 5: R = X^2 + r1 X + r0, so the members are the five roots of
     X^5 + r1 X^4 + (r0 - 1) X^3 - r1 X^2 + (1 - r0) X - n.
     The zeta-family is r0 = 1, r1 = 0. The Dickson case r0 = -4, r1 = 0 (D_5(X,1) = n) needs
     sqrt(f^2 - 4) to be a polynomial, so it gives no polynomial family. It does give the Lucas-type
     algebraic family zeta^j w + zeta^-j w^-1.
 - Conjecture W' (open; evidence: Lemmas 3, 4, 7, 19, Theorem 9 and all LM/Groebner runs). Every complex Williams
   family is "monomial": f_j = c_j u^(e_j) with E/e_j odd (E = sum e_j) and c_j^(E/e_j) = prod c_i.
   Real examples exist, e.g. {v, v, v, v^3, v^3, v^3, v^3} with n = v^15, but their members are powers of
   one polynomial. So W' would imply: no family with distinct irreducible members over Q, for any k and any degree.
   Real case still open: shared NON-real special points on both sides (Lemma 8(a)).

### Negative check: no class-group forcing of Lucas indices [COMPUTED: scratch forced3.gp]
 - Hypothesis. For D = -31 (h = 3, (alpha) = p2^3), alpha = u pi^3 in the Hilbert class field H. Then gamma would
   be a cube in T_p and 3 would never divide z(p) for inert p == 2 mod 3. D = -23 is similar.
 - Test: inert p < 2e5 with 3 | p + 1. The fraction with 3 | z(p) is 0.743-0.758 for all 13 D tested
   (h = 1, 2, 3, 5). This is the random-element value 3/4.
 - Conclusion: the unit u^(1-p) is not a cube in T_p, so there is no forcing. Class field structure gives no free
   Lucas index for any Selfridge D.

### Calibration: the Pomerance pool crosses density 0.6 at T = 1000 [COMPUTED: src/c3/pomgen.c, src/pom/profile.py]
This corrects the session-3 reading "pool density about 0.15, flat in y". That measurement used CHM lists at
y <= 200 only. Direct generation shows the density rising steeply with T.
 - Full pool, D = 5, class (1,2), Pomerance side rule, all primes p <= Xmax. N/log2|G| uses p-1, p+1 exponents.

   | T    | Xmax    | N    | log2|G| | N/log2|G| | saturated?                                          |
   |------|---------|------|---------|-----------|-----------------------------------------------------|
   | 300  | 1e18.2  | 55   | 540.7   | 0.102     | yes (no new primes after 1e12)                      |
   | 500  | 1e18.2  | 156  | 969.3   | 0.161     | yes (no new primes after 1e14)                      |
   | 1000 | 1e15    | 1070 | 2346.7  | 0.456     | no                                                  |
   | 1000 | 1e17    | 1247 | 2507.8  | 0.497     | no (per-decade counts still 99, 78 at 1e15, 1e16)   |

 - True orders (profile.py). Index effects are negligible: 539.7 vs 540.7 bits at T = 300, 2506.8 vs 2507.8 at
   T = 1000. Every pool prime lies in class (1,2) with the expected sides.
 - Greedy frontier: drop the primes attaining the top exponent of one modulus whenever that improves N - bits.
   - T = 300 and 500: the full pool is best (0.100, 0.161).
   - T = 1000 (1e17): best (bits, N) = (1559, 980), ratio 0.629. Further points (1843, 1079) and (1689, 1036);
     the deficit at the best point is 579 bits.
 - Growth: N grew 8x from T = 500 to 1000 while log2|G| grew 2.6x. Extrapolating, the frontier ratio reaches 1 near
   T = 1300-1500, at about 2000-3500 bits with N = 2000-3500 (HEURISTIC extrapolation).
   This is consistent with Chen-Greene's weak-version numbers (density 1 at 1240 bits, 2.6 at 1838 bits), allowing
   for the single-class loss.
 - What it means. Pomerance's heuristic predicts BPSW psp built from pool primes of about 1e10 already in the
   T ~ 1500 pool, roughly 2000-3500 primes and about 10^20000 in size. The pool is not the bottleneck; the
   subset-product search is.
 - Search cost at the crossover (HEURISTIC; standard generalised-birthday analysis).
   - With 2^m leaf lists, Wagner's k-tree needs N >= 2^m B/(m+1) primes and time ~2^(m + B/(m+1)).
   - At B ~ 2500 the minimum over m is about 2^99 time with N >= 2^49 primes. With the actual N ~ 3000, the best
     m is about 3, giving time 2^600.
   - Better random-knapsack algorithms (BCJ-type 2^(0.29 N) at density 1; Flaxman-Przydatek polynomial time only when
     B <= c log^2 N) do not reach B ~ 2500 either.
   - Raising T does not help: N/B grows only polynomially in T while B grows linearly, so B/log(N/B) -> infinity.
 - Sharpened missing ingredient (i): a subset-product algorithm for groups prod (Z/q^e)^* (q <= ~1500) that beats
   generic random-knapsack bounds by hundreds of bits. Equivalently, exploitable non-randomness in the discrete logs
   of two-sided-smooth primes. Known structure: every pool prime is == -1 mod 3, never == -1 mod an F-side q, never
   == +1 mod an L-side q, and is == +-1 on its own support. That frees only the q = 3 coordinate (about 1.6 bits).

### Strategic review (2026-10-24, session 4b, after Theorem C, Lemma 19 and the T = 1000 calibration)
 - Machine rebooted at about 09:00. Lost jobs were restarted with lower memory caps (gp parisizemax <= 3 GB):
   - D = 5 split scan, resumed at 6.31e11 (usplit3 now accepts X0:X);
   - L1 for the remaining positive D;
   - big11 (Theorem B at 1e11, three large cofactors);
   - pomgen calibration at T = 1000/1500/2000;
   - the one-side-automatic family sweep.
 - New since session 3:
   - Theorem C: rigorous exclusion of the whole P_100 pool.
   - Lemma 19: k = 3 Williams families classified over C, plus the coprime normal form.
   - Negative check: no class-group forcing.
   - Calibration: the Pomerance pool density is 0.63 at T = 1000 and rising steeply. This corrects "flat ~0.15".
 - Distance to an explicit n.
   - The pool side is no longer the obstruction: heuristically, BPSW psp exist inside the T ~ 1500 Pomerance pool.
   - The obstruction is now a single, explicit computational problem: a random-like subset-product in a group of
     about 2000-3500 bits.
   - For it, every known algorithm costs >= 2^99:
     - Wagner/Kuperberg/AGHS need N ~ 2^(c sqrt B);
     - KOW-type polynomial-threshold zero-sum algorithms exist only per small field, with constants q^(q^2),
       and need N >= poly(B) >> two-sided pool sizes;
     - "free modulus" variants (all pool primes == 1 mod R) cost a factor R in pool density, which a two-sided pool
       cannot afford.
   - Every non-subset-product mechanism examined is closed or heuristically empty:
     - polynomial families: Lemmas 3, 4, 7, 19, Theorem 9; W' open;
     - exponential families: Attempt C;
     - DSO pools;
     - class-field forcing;
     - Fermat/Mersenne/Lucas automatic families: lottery sweeps running.
 - Theory of success.
   - Constructive: needs either an algebraic family (W' says none over Q) or a structured subset-product algorithm.
     No candidate for either is in hand.
   - Exclusion side: rigorous theorems keep shrinking where small BPSW psp can live:
     - Theorem B: one inert prime, split primes < 1e11;
     - L1: n/q > 1e6;
     - Theorem C: P_100;
     - DSO-closed families.
 - Kill / continue:
   - KILL the C3 y = 150/200 enumeration (CHM lists incomplete, explosive cost, yield ~0). Done 08:40.
   - CONTINUE the exclusion pipeline.
   - CONTINUE the calibration to T = 1500/2000 (gives the crossover directly).
   - CONTINUE E2/W' as theory.
 - Live directions (>= 3):
   E1 Exclusions: Theorem B at 1e11 (9 D) and 2^40 (D = 5); L1 positive D; record.
   E2 Theory: prove W' (complex Williams families are monomial), or at least the real non-coprime case.
   E4 Calibration: measure the Pomerance crossover T*. Estimate the size and factor count of the smallest predicted
      BPSW psp, and publish the explicit crossover pool as a well-defined subset-product instance.
   E5 Algorithms (long shot): a joint zero-sum method for prod_l F_l^(n_l) whose sample threshold is polynomial and
      does not multiply across l. A positive result here is the only thing that would convert E4 into an explicit n.

## Session 5 (2026-10-24, continued after the reboot)
### Theorem B extended to 1e11 for nine D [COMPUTED: src/us/us3_1e11_9D.txt, rows_1e11.gp, run11.gp -> run11.out, biggcd.c -> big7.out, big11b.gp]
 - THEOREM B (1e11). For D* in {-7, -11, 13, -15, 17, -19, -23, 29, -31} there is no BPSW psp with Selfridge D*
   that has exactly one inert prime factor and all of its split prime factors below 1e11. The size of n and of the
   inert prime is unbounded. (D* = 5 stays at 1e10 + SW's two primes until the 2^40 scan finishes.)
 - Inputs. The 1e11 scan finds 61 admissible split primes for the nine D:
   13:12, 29:9, -11:5, -15:8, 17:7, -23:5, -31:6, -7:7, -19:2.
   It reproduces the 1e10 list exactly (cross-validation of the segmented scan). New primes above 1e10:
   -23: 67974381749; -15: 19767865951, 29929785679; -11: 28769448247; -7: 11954925473, 41291876971;
   13: 61620938419; 17: 48695287049, 71309972101, 77059587233; 29: 51354510013.
 - Search. 73 class-consistent, cross-coprime split parts m, 23 candidate inert q, 0 hits.
   Per D (split parts / candidate q): -7: 8/0, -11: 6/1, 13: 15/9, -15: 10/7, 17: 8/3, -19: 2/0, -23: 6/2, 29: 12/1,
   -31: 6/0.
 - Three split parts needed big gcds (the cost filter in run11 skipped them):
   - D = -7, m = 62020897 * 41291876971, O = 3 * 160447363 * 2660227837, class (1,1).
     The four nontrivial unitary splits need gcd(2^W1 - 1, U_W2(1,2)) with W1, W2 up to 1.6e10. GP overflowed a
     3 GB stack, so I wrote src/us/biggcd.c (GMP).
     - Mode N reduces mod 2^W1 - 1 by folding; mode U computes U_W2 exactly.
     - Validated against GP on (-7, 22022, 30030), (13, 1200, 2800), (5, 6000, 2400) and small cases.
     - All four gcds are 1 (49-399 s, at most 3.0 GB RSS).
     - The four trivial splits (min(W1, W2) <= 6) give g = 1 directly.
   - D = -15, m = 4375578271 * 19767865951, class (0,0): 16 splits, 0 candidates.
   - D = -23, m = 1021622741 * 67974381749, class (2,1): 4 splits, 1 candidate, not a psp.
 - COROLLARY (1e11). Every 2-factor BPSW psp n = pq with D* in the nine values above has its split prime > 1e11.

### Compute moved to the GCP VM "dirichlet" (2026-10-24, user directive) — RESUME NOTES
 - All experiments now run on `dirichlet` (the user's GCP VM; 8 cores AMD EPYC 7B12,
   31 GB). This replaces the earlier laptop-only rule. The laptop holds the repo, notebook and analysis.
   - Access: `gcloud compute ssh dirichlet --project <project> --zone <zone>`. Project and zone are kept in the local agent memory, not in this public repo.
   - Work tree: `~/b-psw/src` (copy of src/ without the CHM lists ex100/ex150 and macOS binaries).
   - Rebuilt there: us/usplit3, us/biggcd, c3/pomgen2. Validation: the VM usplit3 scan to 1e8 reproduces all
     145 rows of the laptop 1e11 list below 1e8 (D = -7, -11, 13, -15, 17), byte for byte after sorting.
   - Do not touch `~/f4-gpu-code` (another project). Jobs run under nice 10, at most 7 at once.
 - Queue 1 (src/vm/jobs1.txt, started by vm/setup1.sh; log vm/runq1.log):
   - pomgen2 T = 1000 to 1e18.2, T = 1500 to 1e15, T = 2000 to 1e15 (c3/pomgen_*.out, c3/pom_*.txt);
   - D = 5 split scan [646730219521, 2^40) in four chunks (us/us3_d5_c1..c4.txt). The laptop runs covered
     p < 646730219521 (us/us3_2p40_d5.txt and us/us3_2p40_d5_b.txt; the laptop resume job was stopped);
   - L1 at M0 = 1e6 for D = 41, 53, 61, 73, 89, 97, 101 (lp/lp_1e6_D*.out). D = 37 finished on the laptop (lp_1e6c.out):
     89121 cofactors, 125636 candidate q, 0 hits, unfactored m = 281581, 883411, 890221;
   - autofam2 (cyc/af1..3.out): Mersenne, Wagstaff, Fibonacci sweeps to p = 20000 plus F_5..F_17. Unlike autofam.gp
     there is no ispseudoprime filter. Known prime/PRP indices act as positive controls.
 - W' discriminating patterns (alg/lm_mixed.py). Done on the laptop before the move:
   (1,1,1,1,2), (1,1,1,2,2), (1,1,2,2,2): 150 trials each, 0 solutions (min residuals 3.7e-4, 3.5e-5, 1.1e-4).
   Remaining: (1,2,2,2,2), (1,1,1,3,3), (2,2,2,3,3), (1,1,1,1,1,1,3), (2,2,2,2,2,2,2). They need numpy/scipy on the VM.

### Lemma 20 (PROVED; one step checked for odd k <= 301): W' holds in the coprime normal form when R is constant
Register rows: C4/C3F/C3L (Williams polynomial families). Partial progress on E2.
 - Setting. Take k odd and f_1..f_k in C[t] nonconstant with n = prod f_j, f_j - 1 | n - 1 and f_j + 1 | n + 1, and the
   f_j^2 - 1 pairwise coprime (coprime normal form, session 4). Then
     prod_j (X - f_j) = S(X) - n,   where S(X) = X (X^2 - 1) R(X) + X and R in C[t][X] is monic of degree k - 3.
   Equivalently, the f_j are k distinct polynomials with S(f_1) = ... = S(f_k) = n. Since S(0) = 0 and S(+-1) = +-1,
   this says sum_{odd i < k} e_i(f) = 0 and sum_{even 0 < i < k} e_i(f) = 0.
 - Statement. If R has constant coefficients (R in C[X]), then S = X^k and {f_j} = {zeta^j u} with zeta^k = 1, n = u^k.
 - Proof.
   1. S(f_j) = n with S in C[X] and n nonconstant. So S(X) - n is irreducible over C(n) (it is linear in n), and all
      k of its roots lie in C(t). Its splitting field E satisfies C(n) c E c C(t), and Gal(E/C(n)) = Mon(S).
   2. n is a polynomial in t, so the only place of C(t) over n = infinity is t = infinity. Hence E has a unique
      place over n = infinity. So the decomposition group there is all of Mon(S). Over C it equals the inertia
      group, which is cyclic (tame ramification). So Mon(S) is cyclic.
   3. Mon(S) is transitive on k points, so it is regular of order k, and E = C(f_1).
      - The automorphisms of C(f_1) over C(n) fix the unique place over infinity, so they are affine:
        f_j = l_j f_1 + m_j, with S(l_j X + m_j) = S(X).
      - A degree-k polynomial with an affine symmetry group of order k is a(X - c)^k + b.
   4. S is monic with S(0) = 0, S(1) = 1, S(-1) = -1. Write S = (X + mu)^k - mu^k. Then
      h(mu) = h(-mu) = 0, where h(z) = (1 + z)^k - 1 - z^k.
      - gcd(h(z), h(-z)) = z for every odd k in [3, 301] [COMPUTED: src/wp/mu_check.gp].
      - So mu = 0, S = X^k and f_j = zeta^j f_1.  QED (for k <= 301; step 4 for general k is open but routine).
 - Scope.
   - This is the special case C(e_1(f), ..., e_k(f)) = C(n).
   - The same argument works whenever the subfield K0 = C(e_1(f), ..., e_k(f)) is proper. Then K0 = C(s) with s a
     polynomial, the splitting field is a cyclic Kummer cover C(w), s = w^m + c, and the f_j become polynomials in w
     permuted by w -> zeta_m w.
   - The open case of W' is exactly K0 = C(t): the symmetric functions generate everything, and the Galois argument
     gives nothing.
 - Reformulation for k = 5. W' in the coprime case is equivalent to: distinct f_1..f_5 in C[t] with e3 = -e1 and
   e4 = -e2 force e1 = e2 = 0.
   - The four LM patterns tested so far, (1,1,1,1,2), (1,1,1,2,2), (1,1,2,2,2), (1,2,2,2,2), are exactly the
     Newton-polygon-admissible unequal patterns. The hull vertices lie in {0, 5} union one of {1,3} and one of {2,4},
     because the X^1/X^3 coefficients and the X^2/X^4 coefficients have equal degrees.
