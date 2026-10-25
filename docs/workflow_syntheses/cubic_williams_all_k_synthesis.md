### Real cubics for odd k >= 9: results of workflow cubic-williams-all-k (2026-10-25) -- READ FIRST [lines m9, xn, ns, ak, sf; verifier scripts m9v_*, m9vmf_*, vk_*, vrf_m9*; all files in src/wp/mom/ unless stated]
Register rows: C4/C3F/C3L, real coefficients, pairwise distinct nonconstant members, degree <= 3. Normal form as in Theorem 11:
f_i = 2t + 1 + t(t+1)(a_i t + c_i), a_i > 0, plus at most one linear member 2t + 1. Notation: v = 1/a, u = c/a, w = u - 1/2,
M(a,b) = sum_i v_i^a w_i^b over the cubic members, S = M(1,0), m = number of cubic members.
Two name clashes in this section. m9's S-free odd relation is written Phi_odd. xn's second-order quantity is written Phi. xn's Re(t) is a
function, not a real part.
 - Tags used in this section.
   - PROVED: an exact argument or exact certificate that survived an independent verifier in this workflow (verdict: not refuted).
   - UNAUDITED: the line reports an exact proof or an exact-rational certificate, but nobody else has checked it. Not citable as a theorem.
   - ELEMENTARY: a short identity or argument re-derived by hand while writing this synthesis. It was not in the verifier set.
   - COMPUTED, HEURISTIC and REFUTED as usual. NEGATIVE means a limitation of a method, not an obstruction.
 - Status of degree <= 3 for odd k >= 9.
   - PROVED state: still OPEN for every odd k >= 9. Theorems 10 and 11 close every other k.
   - PROVED: the method of Theorem 11 cannot close k = 9 or k = 11. In all four cases (9,0), (8,1), (11,0), (10,1) the four boundary
     equations have exact real solutions in the open domain with distinct members.
   - UNAUDITED candidate: no family exists for k = 9 (line xn). The proof combines an a priori bound S < 7.3106 / 6.2990 with 109 exact
     Farkas certificates for all smaller S. If it survives an audit, degree <= 3 is open only for odd k >= 11.
     NUMBERING: the line called this "Theorem 12", but Theorem 12 is already the (2,2) quartic pencil result for k = 5. Call it
     "Theorem 13 (candidate)" until it is audited.
   - UNAUDITED, all k (line ak): no family whose coefficient points (a_i, c_i) are collinear (C2); no ray family, i.e. all u_i equal, with or
     without 2t + 1 (C3); and an explicit half-plane criterion (C4).
   - COMPUTED: recall-calibrated searches found no real family at k = 9 or 11, with or without 2t + 1, and no endpoint approaches zero.
   - Constructive relevance: every result here is structural (route II). No candidate family emerged.
 - Outcome of the sign-free audit: NOT REFUTED. The extension of Lemma 3, Theorem 10 and Theorem 11 to arbitrary leading signs was
   re-derived step by step. There is no gap and no repair is needed. Positive leading coefficients are dropped from their hypotheses
   (PROVED). Details are in the sign-free section below.

### Line m9: moment form and feasibility of the boundary hierarchy, k = 9, 11 [m9_*.gp/.py; verifiers m9vmf_indep.gp, m9vmf_neg.gp, m9v_indep.gp, vk_indep.gp, vk_drive.gp, vrf_m9pt.gp, vrf_m9dom.py, vrf_m9all.gp, m9v_dimcheck.py, m9v_t2check.gp, m9v_r2num.py]
 - Moment form of the four boundary equations [PROVED; m9_check.gp, m9_mform.gp -> m9_mform.out; verifier m9vmf_indep.gp, m9vmf_neg.gp,
   m9vmf_rerun_check.out].
   - Definitions. T(s) = prod_i (1 + x_i(s)) mod s^k, with the factor 1 + s/2 for 2t + 1, and V(t) = t^(k-1) T(1/t). For odd k,
     V(-1) = T(-1) and V'(-1) = -(k-1) T(-1) - T'(-1). So the boundary equations are T(-1) = 2^(1-k) and T'(-1) = 0, for the family and
     for its reflection. log T = sum_n Lam_n s^n with Lam_n = sum_i [s^n] log(1 + x_i).
   - Even/odd split. Reflection acts as M(a,b) -> (-1)^b M(a,b), so the system splits into even parts Z0+, Z1+ and odd parts Z0-, Z1-.
     Moment types: 12 for k = 9 (a <= 3, b <= 5); 20 for k = 11 (a >= 1, 2a + b <= 9, top terms M17, M25, M33, M41).
   - Products for k = 9: S^2, S M12, S M20, M11^2 (even); S M11 (odd).
   - Products for k = 11. The line's list was incomplete; corrected:
     - even: S^3, S^2, S{M12, M20, M14, M22, M30}, M11^2, M11{M13, M21}, M12^2, M12 M20, M20^2;
     - odd: S^2 M11, S{M11, M13, M21}, M11{M12, M20}.
   - Verification. The residual is 0 as a polynomial in all 2m variables for k = 5, 7, 9, 11, both cases. The verifier used independent
     code in two ways:
     - an exact Euclidean-quotient V at 40 random rational points per case;
     - an exact symbolic identity with r active members. Since x_i = O(s^3), each monomial involves at most 2 members (k = 9) or 3
       (k = 11), so checking r = R and R + 1 is a proof.
     All 5 deliberate perturbations (negative controls) were caught.
   - Scope: the four boundary equations only, not the 4k-equation family condition.
 - What the V'(-1) equation forces at k = 9 [the claim as first written was REFUTED in two sentences; the corrected items below are
   PROVED by the verifier's independent script m9v_indep.gp -> m9v_indep.out; sources m9_forced.gp -> m9_forced.out, m9_run9.gp,
   m9v_run9_short.gp].
   - E2 = -((k-1) Z0 + Z1).
   - Odd part: (S - 75/4) M11 = 5 M13 - 21 M21 (9 cubics); (S - 111/8) M11 = (9/2) M13 - 19 M21 (8 cubics + 2t+1).
   - Even part, whose only nonlinearity is S^2:
     - 9 cubics: 5/2 S^2 - 351/16 S - 12 M12 + 27 M20 - M14 + 6 M22 - 4 M30 + 255/32 = 0;
     - 8 + lin: 9/4 S^2 - 15 S - 39/4 M12 + 89/4 M20 - M14 + 6 M22 - 4 M30 + 143/32 = 0.
   - Corrected derivation of Phi_odd. Write Z0- = c0 S M11 + ... - B and Z1- = c1 S M11 + ... + 8B, with B = M15 - 8 M23 + 12 M31 and
     (c0, c1) = (9/2, -35) [8 + lin: (4, -31)].
     - Eliminating B gives the odd V' relation 8 Z0- + Z1-, which still contains S.
     - Eliminating S M11 instead gives the S-free relation Phi_odd = -35 Z0- - (9/2) Z1- = 1149/16 M11 + 13 M13 - 51 M21 - B = 0
       [8 + lin: -31 Z0- - 4 Z1- = 759/16 M11 + 21/2 M13 - 41 M21 - B = 0].
     - Phi_odd keeps B and contains no M(a,0). REFUTED: "Phi is obtained by eliminating the shared top block".
   - Unlike k = 7, M11 = 0 is NOT forced: at the certified (8,1) point M11 = 11.16 and S = 8.497.
     - This is PROVED by the feasibility certificates below, not by m9_run9.gp. That script's exact-point search found 0 points in 400000
       tries, so citing it as evidence is REFUTED.
     - "The higher moments absorb the constraints" is interpretation, not a theorem.
   - k = 11 analogues: as printed in m9_forced.out (PROVED; m9v_indep.gp reads and checks them verbatim).
 - Feasibility of the four boundary equations, k = 9 and 11, all four cases [PROVED, computer-assisted, exact-rational Krawczyk;
   m9_kraw.gp, m9_kcert.gp, m9_sol_{9_0,9_1,11_0,11_1}.gp -> m9_kcert_*.out; verifier vk_indep.gp, vk_drive.gp, VM ~/b-psw/src/wp/mom/vk_out/;
   also m9v_dimcheck.py].
   - Statement. Each case has a real point solving the four equations exactly, with:
     - pairwise distinct members;
     - u_i^2 < 8 v_i and (1 - u_i)^2 < 8 v_i, so v_i > 1/32 and no g_i or h_i has a real root.
     Domain margins min_i [8v - max(u^2, (1-u)^2)]: 0.093 (9,0), 0.140 (8,1), 0.083 (11,0), 0.0041 (10,1).
   - Independent checks.
     - The reruns are byte-identical.
     - vk_drive.gp built new certificates, 8/8 PASS. It used a different formula for V (Euclidean quotient over Z[u,v]), different
       denominators (2^30, 3^19), rho = 1e-70, a Taylor-shift Jacobian bound, and two different slices per case.
     - m9v_dimcheck.py (Python Fractions, no shared code) reproduces all four certificates.
   - Consequence: Theorem 11's method (boundary equations + v > 1/32 + K >= 0 + the M bound) cannot exclude k = 9 or k = 11.
   - ELEMENTARY (verifier): put w_i' = v_i - 1/32. Then Mmax - M = 2 sum_{i<j} w_i' w_j' > 0, so the M bound of Theorem 11 follows from
     v_i > 1/32, and K >= 0 is automatic. Theorem 11 stays valid; its k = 7 separation in fact uses only the forced values and v_i > 1/32.
   - Control (COMPUTED): m9_feas.py converges in 0/24 starts for k = 7 (both cases) and in 1/12 for k = 9.
 - Explicit exact point, 9 cubics [PROVED; m9_symexact.gp -> m9_symexact.out; verifier vrf_m9pt.gp, vrf_m9dom.py, vrf_m9all.gp].
   - The point lies in the reflection-symmetric slice: 4 pairs (1/2 +- sqrt(y_j), v_j) plus one self-dual member (1/2, v_c), with
     - v = (15129/25000, 671/12500, 6201/100000, 1193/20000), v_c = 80537/50000;
     - y_2 = (2863/10^5)^2, y_4 = (1143/10^4)^2;
     - y_1 ~ 2.65243 in Q(sqrt D), where D (printed in m9_symexact.out) is not a square;
     - y_3 = 3640457795177541991529/138017997877500000000 - (3585828512/360831367) y_1 ~ 0.017704.
     Then S = 3.1717.
   - Verification. The residuals are [0,0,0,0] exactly under two different V computations. Certified domain margins: >= 0.0953 for the
     pairs, 12.64 for the centre.
   - Mechanism: 70 Z0+ + 9 Z1+ contains no M14, so it is linear in the y's, while Z0+ is quadratic in them.
   - Corrections.
     - Fields: the (y, v) coordinates lie in Q(sqrt D); the u_i lie in Q(sqrt D, sqrt y_1, sqrt y_3).
     - Count: m9_symexact.out has 6 pair choices. 4 of them give feasible points: [1,3] (above), [1,2], [1,4], [2,3]. [2,4] and [3,4]
       solve the equations but have y < 0. So the line's "four other choices" should read "three other choices".
   - Context: the full remainder (F - 1) mod t prod g_i at this point is about 7e12. It is a boundary point only, far from a family.
   - ELEMENTARY (hand evaluation from the exact data, 4 digits). At this point M12 ~ 3.214, M20 ~ 3.347, and
     T2+ = S^2 - 209/4 S - 11 M12 + 23 M20 + 1785/32 ~ -58.2. At level 2, T2+ = (e2(s) + e2(alpha))/2^10 > 0, so the point fails order 2
     by a wide margin. This agrees with the xn cover, which needs the second-order cut for S in [1.33, 7.5].
 - Local dimension [PROVED as a consequence of the certificates; checked by m9v_dimcheck.py].
   - At each certified point the 4 x 2m Jacobian has rank 4, because ||I - YJ(B)|| < 1 on the box. So the boundary-feasible set is locally
     a real-analytic manifold of dimension 2m - 4: 14 (9,0), 12 (8,1), 18 (11,0), 16 (10,1).
   - The symmetric slice meets the feasible set.
     - 9 cubics: PROVED by m9_symexact.
     - 11 cubics: the cited evidence was insufficient. The certified (11,0) point is provably OFF the slice: at box radius 1e-60,
       u11 - 1/2 ~ 1.7e-9 and v1 - v6 ~ 9.9e-10.
     - It is now PROVED by the verifier's 2 x 2 slice certificate (m9v_dimcheck.py, mode "sym m9_sol_11_0.gp"). There the slice-feasible
       part has local dimension 9.
 - Order 2 at t = -1 [PROVED; m9_t2.gp -> m9_t2.out, m9_t2k11.gp -> m9_t2k11.out, m9_hier.py; verifier m9v_t2check.gp, m9v_r2num.py].
   - Setup. H := (F - 1)/(t prod g) = t^(k-1)[prod(1 + x) - prod x], counting g = 2 for 2t + 1. Its polynomial part is V.
     - The "-" half holds iff the Taylor coefficients r_0..r_(2m-1) of H - V at t = -1 vanish (deg R < 2m and prod g(-1) != 0).
     - The "+" half is the same statement for the reflection.
     - The boundary equations are r_0 and r_1 of each half.
   - Order-2 condition. r_2 = 0 iff e_2(s) = 2^k V''(-1) - 2(k-1)(k-2), where s_i = f_i'(-1) = h_i(-1) = 2 + (1 - u_i)/v_i > 0
     (s = 2 for 2t + 1). For the reflection, alpha_i = g_i(0) = 2 + u_i/v_i. The f_j'' terms cancel, and the formula holds for every set
     of members.
   - Exact identity (verified): T2+ = (V''(-1) + V~''(-1))/2 - (k-1)(k-2) 2^(1-k) - (k-1)(k-2) Z0+ - (2k-4) Z1+. At level 2,
     T2+ = (e_2(s) + e_2(alpha))/2^(k+1) > 0.
   - Reduced forms:
     - 9 cubics: T2+ = S^2 - 209/4 S - 11 M12 + 23 M20 + 1785/32;
     - 8 + lin: T2+ = S^2 - 40 S - 10 M12 + 21 M20 + 1113/32;
     - k = 11: in m9_t2k11.out.
     The top moments cancel because they occur only in T_(k-1) and T_(k-2), where the factor (k-1-j)(k-2-j) is 0.
   - Cosmetic fix: the m9_t2.gp comment and the m9_ord2.py docstring have (2k-2)T'(-1); the correct coefficient is (2k-4). No conclusion
     changes, because T'(-1) = 0 on the boundary set.
 - Level 2 (orders 0, 1, 2 for both halves; 6 equations) is FEASIBLE in all four cases [UNAUDITED: exact-rational Krawczyk; m9_kraw2.gp,
   m9_l2cert.gp, m9_l2cert_11_0.gp, m9_l2cert_9.gp -> *.out; starts from m9_lev2.py, m9_t2pos.py, m9_levrand.py].
   - k = 9: certified points exist only at large S, S = 61.9 (9 cubics, margin 0.91) and S = 50.3 (8 + lin, margin 2.19).
     COMPUTED: every level-2 solution found (8 of 40 large starts) has S in [50, 72] and max v in [13, 21], i.e. a in [0.05, 0.08].
     Small-S starts collapse.
   - k = 11: S ~ 3-5 (11 cubics margin 0.067; 10 + lin margin 0.0041).
   - If the certificates are valid, the level-2 set has local dimension 2m - 6.
 - REFUTED, pending the level-2 audit: "sup T2+ < 0 on the k = 9 boundary set".
   - Its only support was about 20 projected-gradient ascents, all converging to a degenerate corner: members coalescing at
     (u, v) = (1/2, 1/32), where g and h acquire double real roots. There 2^9 T2+ = -20166 (9 cubics) and -17709 (8 + lin).
   - The certified large-S level-2 points have T2+ > 0.
   - The even-only relaxation with 0 <= M11^2 <= S M12 is not obstructed either: Psi = 1.4e5 at S = 202 (m9_even.py).
 - Level 3 (8 equations) NOT reached [COMPUTED, negative, optimizer-limited; m9_climb.py, m9_homot.py, m9_levrand.py and their outputs].
   - The homotopy stalls at lam = 0.18 for k = 9 (five members coalescing near (0.437, 0.0399)) and at lam = 0.58 for 11 + lin (members
     near the corner).
   - No conclusion about feasibility.
 - Self-dual rational ansatz for k = 9: 0 points in 400000 tries per case [COMPUTED, uninformative; m9_exact9.gp, m9_run9.gp]. The
   tangent-line rationalisation over-constrains the ansatz.

### Line xn: necessary conditions beyond the boundary equations; candidate closure of k = 9 [xn_*.py/.gp]
 - Reformulation [PROVED; the same statement as m9's r_j equivalence, verified there; xn_core.py exact_check].
   - V is a polynomial iff E(tau) = Vpoly(tau) - (1 - tau)^(k-1) G(tau) vanishes to order 2k at tau = t + 1 = 0.
   - Truncating at order N gives 2N equations: N = 2 is the boundary system, N = 3 is m9's level 2, N = 2k is the full family condition.
 - Second-order identity [PROVED for H: it is m9's r_2 formula, verified by m9v_t2check.gp].
   - V''(-1) = (k-1)(k-2) 2^(1-k) + 2^(-k) e_2(beta), with beta_i = h_i(-1) = 2 + a_i - c_i. This beta is m9's s_i, and beta = 2 for 2t + 1.
   - The reflection uses gamma_i = g_i(0) = 2 + c_i, which is m9's alpha_i.
 - Third-order identity and V(0) [UNAUDITED; xn_ident.gp symbolic for k = 3, 5, 7; xn_ident9.py exact for k = 9].
   - V'''(-1)/6 = -C(k-1,3) 2^(1-k) - (k-1) 2^(-1-k) e_2(beta) + 2^(1-k)[p(e1) p(e2) - p(e1 e2)], with e1_i = -beta_i/2 and
     e2_i = beta_i + gamma_i/2 - 3 - beta_i^2/4.
   - V(0) = sum_i prod_{j != i} 1/(c_j + 2).
 - Positivity identity [ELEMENTARY: re-derived in this synthesis; also checked by xn_ident].
   - Let A_i = f_i - 1 = t g_i and B_i = f_i + 1 = (t+1) h_i, so B - A = 2. Then F - 1 = 2^(-k)[prod(A + B) - prod(B - A)], and for
     odd k
       H = (F - 1)/(t prod g) = 2^(1-k) sum_j t^(k-1-2j) (t+1)^(2j) e_(2j)(h_i/g_i).
   - Every term is >= 0 on R (h, g > 0 and k - 1 - 2j is even). Hence H(t) >= 2^(1-k) t^(k-1) on R, with equality only at t = -1, and
     H(t) - 2^(1-k) t^(k-1) >= 2^(1-k) t^(k-3) (t+1)^2 e_2(h/g).
   - These inequalities hold for the RATIONAL function H and any members. The full family condition enters only through H = V.
 - Bounds on beta and gamma [ELEMENTARY].
   - h(-1) > 0, g(0) > 0, and |a - c|, |c| < sqrt(8a) < 16, so 0 < beta_i, gamma_i < 18.
   - Define Re(t) = (Ve(t) - 2^(-8) t^8)/(t+1)^2 with Ve = (V + V~)/2. For a family Re(t) >= 0.
   - Phi := 2 Re(-1) = (e_2(beta) + e_2(gamma))/2^10, which equals m9's T2+. It lies in (0, 729/32) for 9 cubics and (0, 585/32) for 8 + lin.
 - A priori bound on S [UNAUDITED; xn_Smax.gp; the identity was checked numerically at boundary solutions by xn_Smax_check.py].
   - On the reflection-even boundary set, with B0 = M20:
     - 9 cubics: Re(-7/5) - (59/110) Phi = -(2/55) S^2 - (788/125) S - (2/55) B0 + 39391023/1100000;
     - 8 + lin: Re(-7/5) - (27/50) Phi = -(1/25) S^2 - (678/125) S - (1/25) B0 + 1294023/50000.
   - With Re >= 0, B0 >= 0 and the Phi bound, S < 7.3106 (9 cubics) and S < 6.2990 (8 + lin). The arithmetic was re-checked in this
     synthesis.
   - The bound uses V = H at t = -7/5, i.e. the FULL family condition, not a finite Taylor level.
 - Farkas cover of the remaining range of S [UNAUDITED. Two implementations, both by the line itself: xn_farkas9.py (Bernstein
   subdivision) and xn_certcheck.gp (Sturm counts), output xn_certcheck_prod.out "COVERED, 0 failed"; certificates xn_cert9_l0.gp,
   xn_cert9_l1.gp].
   - The relaxation is a positive measure nu on the domain with:
     - mass m, a cap v <= S - (m-1)/32, and integral of v equal to S;
     - the reflection-even boundary equations, which are linear in nu once S is fixed, with P = M11^2/2 in [0, S M12/2] (Cauchy-Schwarz);
     - the cuts Phi > 0, tangent lower and upper third-order (N3) cuts, and Ve(t0) >= 2^(-8) t0^8.
   - Coverage: 57 certificates cover S in (9/32, 7.5176] (9 cubics) and 52 cover (1/4, 6.34375] (8 + lin). Both ranges contain the range
     left open by the S-bound.
   - Which constraints each part needs:
     - small S: the boundary equations alone (a level-1 statement);
     - S >= 1.33: also Phi > 0;
     - large S: also the positivity and N3-upper cuts.
   - Corrupted certificates are rejected (negative control).
 - CANDIDATE THEOREM 13 [UNAUDITED; the line's name "Theorem 12" collides]: no real Williams family of degree <= 3 exists for k = 9. The
   proof has three parts:
   - the reduction to the cases (9,0) and (8,1) (PROVED: the case split was confirmed by the feasibility verifier);
   - the a priori S-bound;
   - the Farkas cover.
 - REFUTED (a COMPUTED claim of this line, contradicted by m9's level-2 certificates, which are themselves UNAUDITED).
   - The claim: "the max over the boundary-feasible set of min(2^8 D2, 2^8 D2r) is about -10083 (9 cubics) and -8855 (8 + lin)".
   - On the boundary set T2+ = (D2 + D2r)/2, so this is m9's corner value: 2 x (-10083) = -20166 and 2 x (-8855) ~ -17709. Two independent
     optimizers found the same degenerate corner.
   - Likewise, "N = 3 is not reachable by continuation" is an optimizer limitation: m9 has certified N = 3 points at S ~ 50-62.
   - The candidate proof uses neither claim.
 - Requalify: "the measure LP without the N3-upper cut is feasible with Phi > 0 for S >= 11, which is an artifact". Level 2 is genuinely
   feasible at large S (m9), so fractional splitting is not the whole story. This does not affect the candidate theorem, which handles
   large S through the S-bound.
 - Sanity checks.
   - COMPUTED: the four boundary equations alone are numerically infeasible at k = 5 and 7, consistent with Theorem 11.
   - PROVED (by m9): N = 2 is feasible at k = 9.
 - Failed approaches [NEGATIVE; do not retry; xn_dual9.py, xn_lp9.py, xn_bigS.py, xn_Sbound.py]:
   - a single-atom Lagrangian dual (positive only for S >= 2.9);
   - the G-side identity S = V_5 + 56 (bound about 1e10).

### Line ns: direct numerical search and the census of real "-" half solutions, k = 9, 11 [ns_*.py, ns_out/]
 - Encodings validated [COMPUTED; ns_check.py, ns_val_*.py].
   - The per-member remainders agree with exact Fraction division to <= 1.8e-14, and the reflection checks to <= 8e-15.
   - The notebook's k = 5 half solutions refine to 1e-48.
 - The formulation decides recall [COMPUTED; t5cxc.log, ctl_k5_cx_half*.log].
   - Working form: the polynomial identities F - 1 = t G W- and F + 1 = (t+1) H W+, with W eliminated by variable projection and
     constrained to its Laurent shape (k - 5 free coefficients).
   - On the complex k = 5 "-" half this converges from 465/500 starts. The free-W and remainder formulations converge from about none.
 - Degenerate attractors.
   - ELEMENTARY (hand-checked in synthesis): an exact NON-DISTINCT solution of the full k = 7 system. Since
     (2t+1)^3 = 2t + 1 + t(t+1)(8t + 4), the member (a, c) = (8, 4) together with six copies of 2t + 1 gives F = (2t+1)^9, and
     x^3 -+ 1 | x^9 -+ 1. LM converged to exactly this cluster (r7l0.log).
   - ELEMENTARY: (4t^3 + 4t^2 + 2t + 1) together with four copies of 2t + 1 solves the k = 5 "-" half, because
     (2t+1)^2 = 2(2t^2 + 2t + 1) - 1 gives (2t+1)^4 == 1 mod g.
   - COMPUTED:
     - complex solutions at infinity, with all c_i -> infinity;
     - the reality boundary is a plateau where 96.7% of k = 9 endpoints end; boundary deflation raises planted recall at k = 9 from 0.2%
       to 6.8%;
     - at k = 11, an attractor with a near-double large member.
 - Recall calibration [COMPUTED; one planted x* per start].
   - Full real system, with boundary deflation: k = 7 up to 17.4%, k = 9 6.8%, k = 11 2.8%.
   - Full complex system: k = 7 0.4%, k = 9 0/500.
 - Main result [COMPUTED; non-existence is only HEURISTIC]: no real family at k = 9 or 11, with or without 2t + 1.
   - Starts: 12000 for (9,0), 4000 for (8,1), 4500 for (11,0), 1500 for (10,1).
   - The minimum normalized residual is >= 0.80 everywhere, and >= 0.95 over non-degenerate endpoints. Planted zeros converge to
     1e-12..1e-15.
   - A planted-like family would have given about 400 hits at k = 9 and about 126 at k = 11. Observed: 0.
 - Census of real "-" half solutions in D [COMPUTED].
   - Reduction (PROVED by Lemmas 6 and 3): every family is a "-" half solution in D whose reflection is also one. The half system is
     square, so its solutions are isolated.
   - Distinct solutions up to permutation:
     - k = 7: 1; 6 + lin: 1;
     - k = 9: 6 (hit counts 410, 179, 141, 13, 7, 1); 8 + lin: 7;
     - k = 11: 4; 10 + lin: 4. The k = 11 lists are clearly incomplete.
   - Refinement: all 23 were refined at 50 digits ("-" residual <= 4e-44, nonsingular Jacobian, margins 0.008..0.152).
   - None satisfies "+":
     - the relative "+" residual is 288 .. 2.2e7;
     - not even one member's "+" condition holds (best 7.9e-5);
     - no reflection of a listed solution appears in the list.
   - Structure: 1-2 large members (a ~ 4-18, u in (0.1, 0.7)); the rest have a < 2.7 and c in (-1.7, -0.6).
 - Complex families [HEURISTIC for k = 7, INCONCLUSIVE for k = 9]: none found, but recall is too low to conclude. Whether the real
   obstruction at k = 7 is only an inequality effect stays OPEN.

### Line ak: all-k structure in the member plane [ak_*.gp/.py, ak_out/]
Chart P = (V, ut) with V = 1/a and ut = c/a - 1/2. w_i and v_i are the upper roots of g_i and h_i, X = -c/(2a), R = -2V - u(1-u)/2.
 - C1, exact relative-evaluation identities [UNAUDITED; ak_ident.gp residuals 0; ak_pw.py].
   - Re f_l(w_i) - 1 = a_l <n_i^w, P_l - P_i> and Re(-f_l(v_i)) - 1 = -a_l <n_i^v, P_l - P_i>, with n^w = (2X - 1, R), n^v = (2X + 2, R).
   - Linear member: |2w+1|^2 = 8V + 2ut and |2v+1|^2 = 8V - 2ut.
   - Second-order product coefficient: K = 8V ut^2 - 8V^2 - (1/4 - ut^2)^2.
 - C2 (all k, cubic members only): if the coefficient points are collinear, there is no family [UNAUDITED; the logic was re-read in
   synthesis and no gap found; ak_collin_check.py: 26988 families, 0 failures]. The member with maximal V has a special point at which
   every other factor has real part > 1.
 - C3 (all k, with or without 2t + 1): ray families (all u_i equal) are impossible. This includes the families with every member odd in
   t + 1/2 [UNAUDITED].
   - Synthesis re-read: the argument is fine for k >= 3.
   - For k = 2 with {cubic, 2t+1}, the w-test can give |factor| = 1, but even k is excluded by parity anyway.
 - C4, half-plane criterion [UNAUDITED; immediate from C1]: there is no family if some member has all the others in {<n^w, P - P_i> > 0},
   or all of them in {<n^v, P - P_i> < 0}.
 - C5 NEGATIVE: the pairwise-product conjecture (|f_l(w_i) f_l(v_i)| > 1 whenever a_l > a_i) is FALSE, with 4498 violations in 4e5 pairs.
   No extremal rule based on a alone can work. Mechanism: K < 0, e.g. K = -8V^2 - 1/16 on the line ut = 0.
 - C6 NEGATIVE: fully covered configurations exist, in which every non-real special point lies in {|f_j| < 1} for some other member
   [COMPUTED at 60 digits on exact rational data, min |f| <= 0.99949; not interval-certified; ak_cover_cert.py].
   - Example: a = (2834289/100000, 2836507/100000, 1348957/50000), c = (685267/50000, 1449513/100000, 335351/25000).
   - Hence no Theorem 10-type single-special-point argument works for any k >= 3.
   - Such configurations were found only with some a > 25; none with all a <= 25.
 - C7, reduction plus CONJECTURE WS' [the reduction is ELEMENTARY; the conjecture has only COMPUTED evidence, weak for n >= 5].
   - WS': every finite set of >= 2 distinct admissible cubics (optionally with 2t + 1) has a cubic member P whose vectors
     (log|f_Q(w_P)|, log|f_Q(v_P)|), Q != P, all lie in an open half-plane.
   - Reduction: WS' implies that no real cubic family exists, for any k.
   - With nonnegative weights only, the statement fails once 2t + 1 is present (exact: a member with 8V < 1 - 2|ut| is doubly covered).
   - Note: WS' throws away the argument equations, so it is a claim about the SQUARE modulus-only system (2k equations, 2k unknowns). That
     makes it a strong conjecture.
 - C8 [UNAUDITED; finite differences agree to 4e-8]: second-order safety criterion kappa_Gamma < kappa_min(P, N). kappa_min is odd in N, so
   at every point exactly half of all directions are bad (4680-point scan).
 - C9 [COMPUTED]: the modulus-only square system has no non-degenerate real solution for k = 2, 3 (control recall 18% and 23%). For
   k >= 4 the null result is meaningless (control recall 0/100).
 - C10 NEGATIVE: the single-point height H_w fails (105 of 7256 covering pairs). L_w was derived by hand and is not exact-verified.
 - C11 HEURISTIC: cluster rigidity, i.e. near-diagonal families are impossible for every k. The constants are not explicit.
 - C12 COMPUTED: the Laurent Newton solver recovers the k = 5 universal "-" cofactor V = t^4 + 15/4 t + 45/16. Universality at k = 7 is
   UNTESTED (0 solutions found; recall too low).

### Sign-free audit: Lemma 3, Theorem 10, Theorem 11 without positive leading coefficients -- NOT REFUTED [sf_exact.gp -> sf_exact.out, sf_geom.py -> sf_geom_1.out, sf_lm.py, sf_run.sh -> sf_A.out, sf_B.out]
 - Outcome: every step was re-derived and holds. No repair is needed. Recorded as PROVED:
   - Lemma 6, any leading signs. Real special points are shared, with equal values and equal multiplicities, and the value -1 occurs only
     for odd k. Exact check: "prod_{i != j} e_i = 1 for all j" has only constant solutions, k = 2..9.
   - Lemma 3, any slopes, for members affine in one common nonconstant h of any degree.
     - Pullback: A(h) | B(h) in R[t] implies A | B in R[u].
     - The line through (r-, 1) and (r+, -1) is unique.
     - With repetitions allowed, all members are equal.
   - Theorem 10, any leading signs and any k >= 2 (with repetitions allowed, all members are equal).
     - Step 1 must go through the sign-free Lemma 3, not Lemma 23's "slopes of one sign".
     - Step 3 can be replaced by one line: once all tops coincide, every member is affine in g, and the sign-free Lemma 3 applies.
   - Theorem 11, any leading signs.
     - The parity lemma is sign-free, and the even case is Theorem 10.
     - In the odd case, a_i < 0 would make both discriminants (a - c)^2 - 8a and c^2 - 8a positive, giving extra real special points,
       which is forbidden. So a_i > 0, and the k = 5, 7 inputs hold unchanged.
 - Corollary S, sign rigidity [ELEMENTARY: produced by the audit and re-derived in this synthesis; not separately verified].
   - Statement: if a real Williams family (k >= 2, any signs) has a real special point, then every leading sign equals the common value E
     at the largest real special point R, and E = +1 when k is even.
   - Proof: on (R, inf) every |f_i| > 1 and f_i -> E as t -> R+; and E^(k-1) = 1.
   - Consequences.
     - (a) Mixed leading signs force NO real special point, so every member has even degree and |f_i| > 1 on R.
     - (b) All-negative with a real special point forces k odd, and then -f is a family with positive leading coefficients.
     - So every positive-leading theorem extends to families with a real special point. The only genuinely sign-free regime is mixed signs
       with no real special point.
 - Numerics [COMPUTED, weak; used for no conclusion].
   - Mixed-sign real quadratic families, k = 2..7, every pattern tested: 0 zeros.
   - Complex controls were found in 3/8 starts (k = 3), 1/20 (k = 5) and 0/12 (k = 7).
   - The k = 5 control is misfiled as nonzero by sf_lm.py's absolute cutoff 1e-16. Fix: use a cutoff relative to ||F||^2.
 - Still open, with any signs: even degree >= 4 (including mixed quadratic + quartic) and cubics for odd k >= 9. By Corollary S, a
   mixed-sign family in these cases has no real special point.

### Obstruction of each method (precise)
 - Boundary/moment method (Theorem 11; line m9) [PROVED unless marked].
   - At k = 9 and 11, the four boundary equations plus every inequality implied by the open domain are feasible, and the solution set is
     locally a manifold of dimension 2m - 4.
   - The domain inequalities Theorem 11 uses reduce to v_i > 1/32: the M bound is redundant and K >= 0 is automatic.
   - Why k = 5, 7 worked: the few moment types are forced to values that conflict with v > 1/32 (S and T at k = 5; sum vu = S/2, M and K at
     k = 7).
   - Why k = 9 does not: the new types M14, M22, M30, M15, M23, M31 absorb the constraints, and M11 = 0 is not forced.
   - Level 2 adds 2 equations and remains feasible [UNAUDITED].
   - If the level-2 certificates hold, a k = 9 proof must use information beyond Taylor level 2 at S ~ 50-62, and a k = 11 proof must do so
     at S ~ 3-5.
   - At level 3 the numerics are optimizer-limited: members coalesce.
 - Global positivity + Farkas (line xn).
   - It works at k = 9 [UNAUDITED] because Re(t) there involves only S, B0 and A2 = M12, and a single multiplier cancels A2.
   - For k >= 11, Re(t) contains the even parts T_n^e for n <= k - 3, which bring in A4, B2, C0, P and higher moments. Cancelling A2 alone
     then gives no S-bound, and the 10-function column basis must be extended.
   - Whether positive combinations of several Re(t0) and Phi give an S-bound for k >= 11 is untested.
   - The measure relaxation lets the largest member split fractionally. At k = 9 this is compensated only by the N3-upper and positivity
     cuts together with the a priori S-bound.
 - Direct search (line ns).
   - It cannot prove non-existence: recall is measured only on planted solutions (2.8%-17% per start).
   - A family with a much smaller basin, or with a member close to 2t + 1 (a < 0.01), could be missed.
   - The k = 11 census is incomplete (near-double attractor), and the complex search is underpowered.
 - Member-plane geometry (line ak).
   - Single-special-point arguments are impossible for k >= 3: full coverage occurs near the reality boundary for a > 25 (C6).
   - Product and pairwise rules are impossible (C5, K < 0).
   - Nonnegative weights fail once 2t + 1 is present, so weights of either sign are needed (C7).
   - WS' needs a rule that selects a safe hull vertex beyond first order. Every hull vertex is safe at first order (C8). No such rule is
     proved, and the non-perturbative case is untouched.
 - Sign-free audit: no obstruction. The extension is valid, and by Corollary S it adds content only in the mixed-sign case with no real
   special point.

### Cross-line consistency (checked in the synthesis)
 - Same artifact found twice. xn's "max min(2^8 D2, 2^8 D2r) ~ -10083 / -8855" and m9's "2^9 T2+ = -20166 / -17709" are the same number,
   since T2+ = (D2 + D2r)/2 on the boundary set. Both come from one degenerate corner reached by two independent optimizers, and m9's
   level-2 certificates (UNAUDITED) refute both as suprema.
 - The candidate k = 9 proof is consistent with every m9 point:
   - the exact level-1 point has S = 3.17 and T2+ ~ -58, so the Phi > 0 cut removes it, as the xn cover requires for S in [1.33, 7.5];
   - the certified (8,1) level-1 point has S = 8.497 > 6.34;
   - the level-2 points have S = 50-62 > 7.31. By the (UNAUDITED) S-identity, Re(-7/5) < 0 there, i.e. these points violate global
     positivity and must fail some order >= 3 (hand arithmetic).
 - m9 predicted that "every valid S-bound must use orders >= 3". This is exactly what xn does, through V = H at t0 = -7/5.

### Errata to apply elsewhere in the notebook
 - Theorem 11 section, "Open: k >= 9 cubics": add "the boundary equations alone are FEASIBLE at k = 9 and 11 (PROVED; see the m9
   section)".
 - Theorem 11, Ingredient 3: add "the M bound is implied by v_i > 1/32, since Mmax - M = 2 sum_{i<j} (v_i - 1/32)(v_j - 1/32)".
 - Theorem 11 statement and register line: drop the positive-leading-coefficient hypothesis (sign-free audit).
 - Research avenues list, item R1: update it for these results.
 - Numbering: xn's "Theorem 12" becomes "Theorem 13 (candidate)", because Theorem 12 is the (2,2) quartic pencil theorem.
 - m9 claim texts, corrections listed above:
   - the derivation of Phi_odd;
   - the evidence for "M11 not forced";
   - "four other choices" becomes "three";
   - "real quadratic field" becomes "(y, v) in Q(sqrt D)";
   - the evidence for the 11-cubic slice;
   - the k = 11 product list;
   - (2k-2) becomes (2k-4);
   - "prod g(-1) = 2^k" holds only when g = 2 is counted for 2t + 1.
 - Methodological rule, from the refuted "sup T2+ < 0" and xn's -10083. In boundary-moment problems the feasible set reaches very large v
   (a -> 0, members approaching 2t + 1). Before claiming that a supremum is negative:
   - sample S >> m/32;
   - check that the maximiser is not a degenerate corner (coalescing members, v = 1/32).

### Next steps (ranked by plausibility x testability x value; re-rank as results arrive)
 1. AUDIT candidate Theorem 13 (k = 9). Use an independent auditor and new code. Cost: hours. Payoff: k = 9 closed. Kill criterion: any
    failed certificate or invalid cut.
    - (i) Re-derive the S-identity at t0 = -7/5 from the verified k = 9 moment forms in m9_mform.out (PROVED), not from xn_k9mom.gp.
    - (ii) Re-derive the third-order identity and the direction of every tangent N3 cut.
    - (iii) Re-check every certificate with a third positivity test, e.g. exact interval subdivision in (u, v) instead of (r, theta).
    - (iv) Confirm that the measure relaxation contains every real family: the cap, the mass, and open vs closed domain at the S endpoints
      9/32 and 1/4.
 2. Cheap discriminating cross-checks (minutes each, on the VM).
    - (a) Evaluate every xn relaxation constraint at every m9 point (the exact level-1 points and the certified level-1 and level-2
      points). A point that satisfies all of them with S inside a certified interval would refute the cover.
    - (b) At the 23 ns half solutions, evaluate the reflected boundary equations, T2+ and the reflected positivity. The first condition
      that fails identifies the mechanism that kills "+".
    - (c) At m9's k = 11 level-2 points (S ~ 3-5), compute Re(t0) on a grid of t0 < -1. A k = 11 proof of xn's shape needs positivity to
      fail at these points for some t0.
 3. Port the xn pipeline to k = 11; this is the main frontier if step 1 passes.
    - Inputs: the PROVED k = 11 moment forms (m9_mform.out) and the PROVED r_2 formula.
    - S-bound: an LP over nonnegative combinations of Re(t0) at several t0 < -1, Phi and the bounded V'''(-1).
    - Then a Farkas cover with columns v^p W^j of higher weight, checked by two independent checkers.
    - If more is needed, use the strengthened positivity H - 2^(1-k) t^(k-1) >= 2^(1-k) t^(k-3) (t+1)^2 e_2(h/g).
    - Test set: m9's k = 11 level-2 points. Every certificate on their S range must use a positivity or N3 cut.
 4. Look for a uniform-in-k mechanism, the missing ingredient for all odd k. Two concrete leads:
    - (a) The exact positivity identity together with 0 < beta, gamma < 18 bounds every Taylor coefficient of V at -1 by polynomials in
      beta and gamma. Look for an S-bound sum_t0 lambda(t0) Re(t0) whose structure does not depend on k.
    - (b) m9 observed that the reduced T2+ at k has the product skeleton S(M12 - 2 M20) + M11^2/2 of Z0+ at k - 2. Test exactly in gp, at
      k = 9, 11, 13, whether Taylor order j at k equals the boundary system at k - 2j up to linear terms. Such a recursion would give an
      induction on k.
 5. WS' (line ak): try to falsify it first.
    - Solve the square modulus-only system (C9) at k = 5..9 with planted-recall calibration. A real modulus-only solution in D refutes WS'.
    - Run DE at n = 5-8 on multi-cluster and large convex configurations.
    - If WS' survives: prove the small-diameter case (oddness of kappa_min plus C11 rigidity), then the spread case with C1.
    - In parallel: audit C2-C4 (cheap, and they are all-k theorems), and extend C2 to collinear cubics plus 2t + 1.
 6. Certified census (line ns), a proof route independent of the moments.
    - Pilot at k = 7 (1 known half solution): certified enumeration of all real "-" half solutions in the closure of D, by interval
      branch-and-bound or certified homotopy, followed by the "+" check. Estimate the cost at k = 7 before attempting k = 9.
    - Complete the k = 11 census (pair deflation, 20000+ starts) to feed step 2(b).
 7. m9 level 3 with a repulsion term that penalises member distance < delta: map S_min(J) for k = 9, J = 2, 3, and try Krawczyk at level 3
    for k = 11. Mainly useful as test data for step 3.
 8. Parked: complex families (monodromy enumeration of complex "-" half solutions at k = 5, 7) and k = 7 cofactor universality.
 - Housekeeping.
   - Fix the sf_lm.py cutoff and the (2k-4) comments.
   - Collect the VM jobs still running: ns x9l1s (complex k = 9 + lin, log ns_out/x9l1s.log), and the ak_cover4 loops (timeout 700 s each)
     writing ak_cover4_{full,wt,wta}.out.
   - Nothing from this workflow is committed, and the lines did not edit NOTEBOOK.md or LOG.jsonl.

[Not for the notebook] Every relative path above resolves under src/wp/mom/ (ns logs in its ns_out/, ak outputs in its ak_out/).