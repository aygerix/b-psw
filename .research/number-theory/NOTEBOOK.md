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
