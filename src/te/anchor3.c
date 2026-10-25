/* anchor3.c -- Theorem E engine: two-torus anchor enumeration.  Specification: src/te/SPEC.md.  Change log: FIXES block at end of file.
 *
 * Build (Makefile line):   anchor3: anchor3.c ; gcc -O2 -fopenmp -o anchor3 anchor3.c -lgmp -lm
 *
 * Usage:  anchor3 --D 5 --X 1e8 [--mode bpsw|weak|plant|fermat|lucas] [--cap K] [--split S A B]...
 *                 [--plants FILE] [--threads T] [--nosieve] [--sieveB B] [--ctest c] [--cpair c]
 *                 [--chunkA n] [--chunkC n] [--progress] [--quiet] [--selftest] [--bench]
 *         anchor3 --D 5 --sample LO:HI   (cost sampling: table and p3 restricted to [LO, HI], branch A only)
 *         anchor3 --D 5 --X 1e10 --planonly   (exact branch plan + CPU estimate from the cost model, no search)
 *   --cap K   : K = 0 forces branch C for every p3, K = -1 forces branch A for every p3 (K = infinity of the
 *               spec), K > 0 puts p3 in branch A iff its branch-A candidate count is < K.  Default: cost
 *               balance, branch A iff count_A * (ctest or ctestW for m >= 2^63) <= (#inert p2 < p3) * cpair.
 *   --split S A B : one pass with split part S, A_S = A, B_S = B (repeatable; default one pass S=A=B=1).
 *   --nosieve : disable the small-prime factorisation of the Fermat-mod-m prefilter (control: prefilter
 *               survivor counts must be identical with and without it).
 *   --quiet   : silences control-mode hit lines (FERMAT/LUCAS/WEAK/PLANTFOUND) only; bpsw-mode HIT lines always print (FIX F1).
 * ======================================================================================================
 * COMPLETENESS ARGUMENT (restated from SPEC.md, with pointers to the code; tags [..] are greppable)
 * ======================================================================================================
 * Target shape: n = p1 p2 p3 S, p1 < p2 < p3 <= X inert primes for D not dividing 2QD, S the split part,
 * n passing exact strong BPSW with Selfridge discriminant D.  Proved necessary conditions (SPEC, F1, F2):
 *   (W1) for each p_i:  M(p_i) | n/p_i - 1,  M(p) = lcm(ord_p 2, z_D(p));
 *   (W2) n == 1 mod A_S and n == -1 mod B_S.
 * The table (build_table/seg_process) holds EVERY prime 3 <= p <= X with (D/p) = -1 and p not dividing Q
 * (tag [T-INERT]); ord and z are exact (prime factors of p -/+ 1 stripped while 2^(e/l)=1 / U_(e/l)=0,
 * tag [T-ORD], [T-Z]).  Every inert p3 of the table is planned exactly once (run_pass block loop, tag
 * [P-ALL]) and goes to branch A or branch C (decide(), tag [P-DECIDE]); results must not depend on K.
 *
 * Branch A (anchor p3).  setupA() solves  m S == 1 (mod M(p3))  [W1 at p3]  and  m p3 S == nu (mod W),
 * W = lcm(A_S,B_S), nu = CRT(1 mod A_S, -1 mod B_S)  [W2], i.e. m == r (mod L).  Every m = p1 p2 of a
 * solution lies in this class, is odd, satisfies 15 <= m (p1 >= 3, p2 >= 5) and m < p3^2 (p1 < p2 < p3).
 *   - progression bound m < p3^2 ............ tag [A-BOUND]  (count = (bound-1-m0)/step + 1, bound = p3^2)
 *   - lower bound m >= 15 .................... tag [A-M15]
 *   - odd m only (step 2L when L is odd) ..... tag [A-ODD]   (even m gives even n: no mode accepts it)
 *   - chunks [k0,k1) of one progression cover 0..count-1 exactly once (tag [A-CHUNK]).
 * Branch C (tail p3).  For EVERY inert p2 < p3 of the table (tag [C-ALLP2]; items cover j = 0..i3-1),
 *   p1 == (p2 S)^-1 (mod M(p3)), p1 == (p3 S)^-1 (mod M(p2)), p1 p2 p3 S == nu (mod W).
 *   The congruences are CRT-combined in this order until the modulus L reaches p2 (then at most one p1 < p2
 *   remains); congruences not combined are verified directly on every enumerated p1 (tags [C-DEF2],
 *   [C-DEFW]).  So the enumerated set contains every solution p1.
 *   - progression bound p1 < p2 .............. tag [C-P1LTP2]
 *   - lower bound p1 >= 3 .................... tag [C-P1GE3]
 *   - p1 must be an inert table prime ........ tag [C-P1TAB]; (W1) at p1 checked .......... tag [C-W1P1]
 * Every 'continue' / early return that skips work, and its justification:
 *   [JA1] gcd(S, M(p3)) > 1: m S == 1 (mod M(p3)) has no solution (missing inverse).
 *   [JA2] gcd(p3 S, W) > 1: m p3 S == nu (mod W) has no solution since gcd(nu, W) = 1 (missing inverse).
 *   [JA3] the two congruences for m are CRT-incompatible: (W1 at p3) and (W2) cannot both hold.
 *   [JA4] every m of the class is even: n even, rejected by every mode's final predicate.
 *   [JA5] sieve: l | m and ord_l(2) does not divide n - 1, so 2^(n-1) != 1 mod l, hence != 1 mod m: the
 *         allowed Fermat-mod-m prefilter would reject this m (factorised evaluation of the same prefilter;
 *         control: --nosieve must give identical prefilter-survivor counts).  sieveA() itself never skips a
 *         candidate except by marking it under exactly this condition (its 'continue's skip marks only).
 *         With S = 1 the prefilter runs in 4 independent lanes (queue + flush); every queued m is tested.
 *   [JA6] prefilter 2^(n-1) != 1 (mod m) (bpsw, fermat modes; implied by the final predicate).
 *   [JA7] control modes: m is not p1 p2 with p1 < p2 < p3 inert table primes (shape requirement of the
 *         control's final predicate; never used in bpsw mode).
 *   [JP0] decide(): setupA() proved the system unsolvable (JA1-JA4) -> p3 skipped (not with --cap 0).
 *   [JC1] gcd(S, M(p3)) > 1: p1 p2 S == 1 (mod M(p3)) unsolvable for every p2 (missing inverse).
 *   [JC2] gcd(p3 S, W) > 1: p1 p2 p3 S == nu (mod W) unsolvable (missing inverse).
 *   [JC3] gcd(p2, M(p3)) > 1: p1 p2 S == 1 (mod M(p3)) unsolvable (missing inverse).
 *   [JC4] gcd(p3 S, M(p2)) > 1: (W1) at p2, p1 p3 S == 1 (mod M(p2)), unsolvable (missing inverse).
 *   [JC5] (W1 at p3) and (W1 at p2) CRT-incompatible.
 *   [JC6] gcd(p2, W) > 1: (W2) unsolvable (missing inverse).      [JC7] (W2) CRT-incompatible.
 *   [C-DEF2] deferred (W1) at p2 fails.   [C-DEFW] deferred (W2) fails.
 *   [C-P1TAB] p1 not an inert table prime (shape).   [C-W1P1] (W1) at p1 fails.
 * Pointers to code lines (line numbers in this file):
 *   [T-INERT] 502  [T-ORD] 524  [T-Z] 533  [P-ALL] 864  [P-DECIDE] 742  [A-BOUND] 735  [A-M15] 730
 *   [A-ODD] 732  [A-CHUNK] 874  [C-ALLP2] 880  [C-P1LTP2] 829  [C-P1GE3] 827  [C-P1TAB] 834  [C-W1P1] 836
 *   [C-DEF2] 831  [C-DEFW] 832  [JA1] 703  [JA2] 709,716,725  [JA3] 711,718,725  [JA4] 733  [JA5] 761,779
 *   [JA6] 770,782  [JA7] 773,793  [JP0] 743,871  [JC1] 802  [JC2] 804  [JC3] 809  [JC4] 814  [JC5] 816
 *   [JC6] 821  [JC7] 823
 * ======================================================================================================
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <gmp.h>
#include <omp.h>
#include <sys/resource.h>

typedef uint64_t u64; typedef uint32_t u32; typedef uint16_t u16; typedef uint8_t u8;
typedef unsigned __int128 u128; typedef __int128 i128;

static double wall(void){ struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + 1e-9 * ts.tv_nsec; }
static double tcpu(void){ struct timespec ts; clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ts); return ts.tv_sec + 1e-9 * ts.tv_nsec; }
static double pcpu(void){ struct rusage ru; getrusage(RUSAGE_SELF, &ru); return ru.ru_utime.tv_sec + ru.ru_stime.tv_sec + 1e-6 * (ru.ru_utime.tv_usec + ru.ru_stime.tv_usec); }
static void die(const char *msg){ fprintf(stderr, "FATAL: %s\n", msg); fflush(stdout); exit(2); }
static char *u128str(u128 x, char *buf){ char t[48]; int k = 0; if(!x) t[k++] = '0'; while(x){ t[k++] = (char)('0' + (int)(x % 10)); x /= 10; }
  for(int i = 0; i < k; i++) buf[i] = t[k - 1 - i]; buf[k] = 0; return buf; }
static void mpz_set_u128(mpz_t z, u128 x){ mpz_set_ui(z, (unsigned long)(u64)(x >> 64)); mpz_mul_2exp(z, z, 64); mpz_add_ui(z, z, (unsigned long)(u64)x); }
static int mpz_get_u128(u128 *x, const mpz_t z){
  if(mpz_sgn(z) < 0 || mpz_sizeinbase(z, 2) > 128) return 0;
  u64 w[2] = {0, 0}; size_t cnt = 0; mpz_export(w, &cnt, -1, sizeof(u64), 0, 0, z);
  *x = ((u128)w[1] << 64) | w[0]; return 1; }

/* ================================ basic integer arithmetic ================================ */
static inline u64 gcd64(u64 a, u64 b){ while(b){ u64 t = a % b; a = b; b = t; } return a; }
static u128 gcd128(u128 a, u128 b){
  while(b){ if(!(a >> 64) && !(b >> 64)) return gcd64((u64)a, (u64)b); u128 t = a % b; a = b; b = t; }
  return a; }
/* (a*b) mod N, 1 <= N < 2^127 */
static u128 mulmod128(u128 a, u128 b, u128 N){
  if(N == 1) return 0;
  if(a >= N) a %= N; if(b >= N) b %= N;
  if(!(N >> 64)) return ((u128)(u64)a * (u64)b) % N;
  if(!(N >> 95)){ u128 r = 0; for(int s = 64; s >= 0; s -= 32){ u128 c = (b >> s) & 0xffffffffu; r = ((r << 32) + a * c) % N; } return r; }
  u128 r = 0; for(int i = 126; i >= 0; i--){ r <<= 1; if(r >= N) r -= N; if((b >> i) & 1){ r += a; if(r >= N) r -= N; } } return r; }
/* inverse of a mod N (1 <= N < 2^127); returns 1 iff gcd(a, N) = 1 */
static int invmod128(u128 a, u128 N, u128 *res){
  if(N == 1){ *res = 0; return 1; }
  if(a >= N) a %= N;
  u128 r0 = N, r1 = a; i128 s0 = 0, s1 = 1;
  while(r1){
    u128 q, t;
    if(!(r0 >> 64)){ u64 x = (u64)r0, y = (u64)r1; u64 qq = x / y; q = qq; t = x - qq * y; }
    else { q = r0 / r1; t = r0 - q * r1; }
    r0 = r1; r1 = t; i128 ts = s0 - (i128)q * s1; s0 = s1; s1 = ts;
  }
  if(r0 != 1) return 0;
  i128 v = s0 % (i128)N; if(v < 0) v += (i128)N; *res = (u128)v; return 1; }
/* x == r1 (L1), x == r2 (L2), r_i < L_i, 1 <= L_i < 2^127. 0: incompatible; 1: ok (r, L < 2^127); 2: lcm >= 2^127 */
static int crt128(u128 r1, u128 L1, u128 r2, u128 L2, u128 *r, u128 *L){
  u128 g = gcd128(L1, L2);
  if(r1 % g != r2 % g) return 0;
  u128 L2g = L2 / g;
  if(L2g == 1){ *r = r1; *L = L1; return 1; }
  if(L1 > ((((u128)1) << 127) - 1) / L2g) return 2;
  u128 r1m = r1 % L2, dd = (r2 >= r1m) ? r2 - r1m : r2 + L2 - r1m, inv;
  if(!invmod128((L1 / g) % L2g, L2g, &inv)) die("crt128: internal inverse failure");
  u128 t = mulmod128((dd / g) % L2g, inv, L2g);
  *r = r1 + L1 * t; *L = L1 * L2g; return 1; }
/* GMP CRT: same contract, returns 0 if incompatible */
static int crt_mpz(mpz_t r, mpz_t L, const mpz_t r1, const mpz_t L1, const mpz_t r2, const mpz_t L2){
  mpz_t g, d, L1g, L2g, inv; mpz_inits(g, d, L1g, L2g, inv, NULL); int ok = 1;
  mpz_gcd(g, L1, L2); mpz_sub(d, r2, r1);
  if(!mpz_divisible_p(d, g)){ ok = 0; goto out; }
  mpz_divexact(L2g, L2, g); mpz_divexact(L1g, L1, g);
  if(mpz_cmp_ui(L2g, 1) == 0){ mpz_mod(r, r1, L1); mpz_set(L, L1); goto out; }
  if(!mpz_invert(inv, L1g, L2g)) die("crt_mpz: internal inverse failure");
  mpz_divexact(d, d, g); mpz_mul(d, d, inv); mpz_mod(d, d, L2g);
  mpz_mul(L, L1, L2g); mpz_mul(d, d, L1); mpz_add(r, r1, d); mpz_mod(r, r, L);
out: mpz_clears(g, d, L1g, L2g, inv, NULL); return ok; }
static u32 inv_u32(u32 a, u32 m){ /* gcd(a,m)=1 assumed, m >= 2 */
  int64_t r0 = m, r1 = a % m, s0 = 0, s1 = 1;
  while(r1){ int64_t q = r0 / r1, t = r0 - q * r1; r0 = r1; r1 = t; t = s0 - q * s1; s0 = s1; s1 = t; }
  if(r0 != 1) return 0; s0 %= (int64_t)m; if(s0 < 0) s0 += m; return (u32)s0; }

/* ================================ Montgomery, odd modulus < 2^63 ================================ */
typedef struct { u64 n, ninv, r1, r2; } m64;
static inline void m64_init_light(m64 *c, u64 n){ u64 inv = n; for(int i = 0; i < 5; i++) inv *= 2 - n * inv;
  c->n = n; c->ninv = (u64)0 - inv; c->r1 = (u64)((((u128)1) << 64) % n); c->r2 = 0; }
static inline void m64_init(m64 *c, u64 n){ m64_init_light(c, n); c->r2 = (u64)(((u128)c->r1 * c->r1) % n); }
static inline u64 m64_redc(u128 t, const m64 *c){ u64 m = (u64)t * c->ninv; u128 u = (t + (u128)m * c->n) >> 64; u64 r = (u64)u; return r >= c->n ? r - c->n : r; }
static inline u64 m64_mul(u64 a, u64 b, const m64 *c){ return m64_redc((u128)a * b, c); }
static inline u64 m64_add(u64 a, u64 b, const m64 *c){ u64 s = a + b; return s >= c->n ? s - c->n : s; }
static inline u64 m64_sub(u64 a, u64 b, const m64 *c){ return a >= b ? a - b : a + c->n - b; }
static inline u64 m64_half(u64 a, const m64 *c){ return (a & 1) ? (a >> 1) + (c->n >> 1) + 1 : (a >> 1); }
static inline u64 m64_to(u64 a, const m64 *c){ return m64_mul(a % c->n, c->r2, c); }
static inline u64 m64_from(u64 a, const m64 *c){ return m64_redc((u128)a, c); }
static u64 m64_pow(u64 am, u64 e, const m64 *c){ u64 r = c->r1; while(e){ if(e & 1) r = m64_mul(r, am, c); am = m64_mul(am, am, c); e >>= 1; } return r; }
/* 2^e == 1 mod n ? (single-limb exponent, context given) */
static int m64_pow2_is1_c(u64 e, const m64 *c){
  if(!e) return 1;
  int b = 63 - __builtin_clzll(e); u64 x = m64_add(c->r1, c->r1, c);
  for(b--; b >= 0; b--){ x = m64_mul(x, x, c); if((e >> b) & 1) x = m64_add(x, x, c); }
  return x == c->r1; }
/* 2^e == 1 mod n ? e = little-endian limbs e[0..nl-1]; n odd, 1 < n < 2^63 */
static int m64_pow2_is1(const u64 *e, int nl, u64 n){
  m64 c; m64_init_light(&c, n);
  int w = nl - 1; while(w > 0 && !e[w]) w--;
  u64 word = e[w]; if(!word) return 1;
  int b = 63 - __builtin_clzll(word); u64 x = m64_add(c.r1, c.r1, &c);
  for(b--;;){
    for(; b >= 0; b--){ x = m64_mul(x, x, &c); x = m64_add(x, x & ((u64)0 - ((word >> b) & 1)), &c); }
    if(--w < 0) break; word = e[w]; b = 63;
  }
  return x == c.r1; }
/* 4 independent lanes (ILP): res[i] = [2^E_i == 1 mod n_i], E_i = (e1[i]:e0[i]) < 2^128, n_i odd in (1, 2^63).
   Leading zero bits are harmless: x starts at Montgomery 1 and 1^2 = 1. */
static void m64_pow2_is1_x4(const u64 *nv, const u64 *e0, const u64 *e1, int *res){
  m64 c[4]; u64 x[4]; u128 E[4]; u128 Eor = 0;
  for(int i = 0; i < 4; i++){ m64_init_light(&c[i], nv[i]); x[i] = c[i].r1; E[i] = ((u128)e1[i] << 64) | e0[i]; Eor |= E[i]; }
  int nb = (Eor >> 64) ? 128 - __builtin_clzll((u64)(Eor >> 64)) : (Eor ? 64 - __builtin_clzll((u64)Eor) : 0);
  for(int b = nb - 1; b >= 0; b--){
    #pragma GCC unroll 4
    for(int i = 0; i < 4; i++){ u64 y = m64_mul(x[i], x[i], &c[i]); x[i] = m64_add(y, y & ((u64)0 - (u64)((E[i] >> b) & 1)), &c[i]); }
  }
  for(int i = 0; i < 4; i++) res[i] = (x[i] == c[i].r1); }
/* U_k(1,Q) == 0 mod n ?  Dm, Qm in Montgomery form */
static int m64_lucasU_is0(u64 k, u64 Dm, u64 Qm, const m64 *c){
  if(!k) return 1;
  u64 U = 0, V = m64_add(c->r1, c->r1, c), Qk = c->r1;
  for(int b = 63 - __builtin_clzll(k); b >= 0; b--){
    U = m64_mul(U, V, c);
    V = m64_sub(m64_mul(V, V, c), m64_add(Qk, Qk, c), c);
    Qk = m64_mul(Qk, Qk, c);
    if((k >> b) & 1){ u64 U1 = m64_half(m64_add(U, V, c), c); u64 V1 = m64_half(m64_add(m64_mul(Dm, U, c), V, c), c); U = U1; V = V1; Qk = m64_mul(Qk, Qm, c); }
  }
  return U == 0; }
static u64 s64mod(int64_t a, u64 n){ int64_t r = a % (int64_t)n; if(r < 0) r += (int64_t)n; return (u64)r; }

/* ================================ Montgomery, odd modulus in [2^63, 2^126) ================================ */
typedef struct { u128 n, ninv, r1; u64 n0, n1, ninv0; } m128;
static inline void mul128x128(u128 a, u128 b, u128 *lo, u128 *hi){
  u64 a0 = (u64)a, a1 = (u64)(a >> 64), b0 = (u64)b, b1 = (u64)(b >> 64);
  u128 p00 = (u128)a0 * b0, p01 = (u128)a0 * b1, p10 = (u128)a1 * b0, p11 = (u128)a1 * b1;
  u128 mid = (p00 >> 64) + (u64)p01 + (u64)p10;
  *lo = (mid << 64) | (u64)p00; *hi = p11 + (p01 >> 64) + (p10 >> 64) + (mid >> 64); }
static inline void m128_init(m128 *c, u128 n){ u128 inv = n; for(int i = 0; i < 6; i++) inv *= 2 - n * inv; c->n = n; c->ninv = (u128)0 - inv; c->r1 = ((u128)0 - n) % n;
  c->n0 = (u64)n; c->n1 = (u64)(n >> 64); c->ninv0 = (u64)c->ninv; }
/* reference Montgomery product (R = 2^128), full 128-bit REDC: used by the self-test only */
static inline u128 m128_mul_ref(u128 a, u128 b, const m128 *c){
  u128 tlo, thi, mlo, mhi; mul128x128(a, b, &tlo, &thi);
  u128 m = tlo * c->ninv; mul128x128(m, c->n, &mlo, &mhi);
  u128 r = thi + mhi + (tlo != 0); return r >= c->n ? r - c->n : r; }
/* production: word-by-word (CIOS) Montgomery product, same R = 2^128; a, b < n < 2^126 */
static inline u128 m128_mul(u128 a, u128 b, const m128 *c){
  u64 a0 = (u64)a, a1 = (u64)(a >> 64), b0 = (u64)b, b1 = (u64)(b >> 64), t0, t1, t2, u;
  u128 p = (u128)a0 * b0, q;
  t0 = (u64)p; q = (p >> 64) + (u128)a0 * b1; t1 = (u64)q; t2 = (u64)(q >> 64);
  u = t0 * c->ninv0; p = (u128)u * c->n0 + t0;
  q = (p >> 64) + (u128)u * c->n1 + t1; t0 = (u64)q; q = (q >> 64) + t2; t1 = (u64)q; t2 = (u64)(q >> 64);
  p = (u128)a1 * b0 + t0; t0 = (u64)p; q = (p >> 64) + (u128)a1 * b1 + t1; t1 = (u64)q; t2 += (u64)(q >> 64);
  u = t0 * c->ninv0; p = (u128)u * c->n0 + t0;
  q = (p >> 64) + (u128)u * c->n1 + t1; t0 = (u64)q; q = (q >> 64) + t2; t1 = (u64)q;
  u128 r = ((u128)t1 << 64) | t0; return r >= c->n ? r - c->n : r; }
static inline u128 m128_add(u128 a, u128 b, const m128 *c){ u128 s = a + b; return s >= c->n ? s - c->n : s; }
static int m128_pow2_is1(const u64 *e, int nl, u128 n){
  m128 c; m128_init(&c, n);
  int w = nl - 1; while(w > 0 && !e[w]) w--;
  u64 word = e[w]; if(!word) return 1;
  int b = 63 - __builtin_clzll(word); u128 x = m128_add(c.r1, c.r1, &c);
  for(b--;;){
    for(; b >= 0; b--){ x = m128_mul(x, x, &c); x = m128_add(x, x & ((u128)0 - ((word >> b) & 1)), &c); }
    if(--w < 0) break; word = e[w]; b = 63;
  }
  return x == c.r1; }
static void m128_pow2_is1_x4(const u128 *nv, const u64 *e0, const u64 *e1, int *res){
  m128 c[4]; u128 x[4]; u128 E[4]; u128 Eor = 0;
  for(int i = 0; i < 4; i++){ m128_init(&c[i], nv[i]); x[i] = c[i].r1; E[i] = ((u128)e1[i] << 64) | e0[i]; Eor |= E[i]; }
  int nb = (Eor >> 64) ? 128 - __builtin_clzll((u64)(Eor >> 64)) : (Eor ? 64 - __builtin_clzll((u64)Eor) : 0);
  for(int b = nb - 1; b >= 0; b--){
    #pragma GCC unroll 4
    for(int i = 0; i < 4; i++){ u128 y = m128_mul(x[i], x[i], &c[i]); x[i] = m128_add(y, y & ((u128)0 - (u128)((E[i] >> b) & 1)), &c[i]); }
  }
  for(int i = 0; i < 4; i++) res[i] = (x[i] == c[i].r1); }
/* Fermat-mod-m prefilter: 2^e == 1 (mod m), m odd, 15 <= m < 2^126 */
static int fermat_m(u128 m, const u64 *e, int nl){ if(!(m >> 63)) return m64_pow2_is1(e, nl, (u64)m); return m128_pow2_is1(e, nl, m); }

/* ================================ Jacobi symbol ================================ */
static int jacobi_s64(int64_t a, u64 n){ /* n odd > 0 */
  u64 x = s64mod(a, n), y = n; int r = 1;
  while(x){ while(!(x & 1)){ x >>= 1; u64 t = y & 7; if(t == 3 || t == 5) r = -r; } u64 t = x; x = y; y = t; if((x & 3) == 3 && (y & 3) == 3) r = -r; x %= y; }
  return y == 1 ? r : 0; }

/* ================================ globals ================================ */
enum { MODE_BPSW, MODE_WEAK, MODE_PLANT, MODE_FERMAT, MODE_LUCAS };
static const char *MODENAME[] = { "bpsw", "weak", "plant", "fermat", "lucas" };
static int gMode = MODE_BPSW;
static int64_t gD = 5, gQ = -1;
static u64 gX = 0, gX0 = 3;
static int gThreads = 4;
static long long gCap = -2;               /* -2 auto, -1 all A, 0 all C, >0 threshold */
static int gSieve = 1; static u32 gSieveB = 4096;
/* cost model (CPU microseconds, measured on the VM): branch-A candidate with m < 2^63 / m >= 2^63 (--sample windows at 1e9..1e10),
   branch-C pair (--cap 0 at X = 1e6) */
static double gCtest = 0.21, gCtestW = 0.63, gCpair = 0.07;
static u64 gChunkA = 1u << 15, gChunkC = 1u << 16;
static u64 gSeg = (u64)1 << 18;
static int gProgress = 0;
static int gQuiet = 0;

/* table of inert primes */
static u64 *Tp = NULL; static u32 *Ti1 = NULL, *Ti2 = NULL; static size_t NT = 0, TCap = 0;
typedef struct { size_t idx; u64 i1, i2; } ovent;
static ovent *Tov = NULL; static size_t NTov = 0, CapTov = 0;
static u32 *BP = NULL; static size_t NBP = 0; static u64 BPmax = 0;

static inline u64 tab_i1(size_t i){ u32 v = Ti1[i]; if(v) return v; size_t lo = 0, hi = NTov; while(lo < hi){ size_t m = (lo + hi) / 2; if(Tov[m].idx < i) lo = m + 1; else hi = m; } if(lo >= NTov || Tov[lo].idx != i) die("overflow table miss"); return Tov[lo].i1; }
static inline u64 tab_i2(size_t i){ u32 v = Ti2[i]; if(v) return v; size_t lo = 0, hi = NTov; while(lo < hi){ size_t m = (lo + hi) / 2; if(Tov[m].idx < i) lo = m + 1; else hi = m; } if(lo >= NTov || Tov[lo].idx != i) die("overflow table miss"); return Tov[lo].i2; }
static inline u64 tab_ord(size_t i){ return (Tp[i] - 1) / tab_i1(i); }
static inline u64 tab_z(size_t i){ return (Tp[i] + 1) / tab_i2(i); }
static int64_t tab_find(u64 p){ size_t lo = 0, hi = NT; while(lo < hi){ size_t m = (lo + hi) / 2; if(Tp[m] < p) lo = m + 1; else hi = m; } return (lo < NT && Tp[lo] == p) ? (int64_t)lo : -1; }

/* plant overrides */
typedef struct { u64 p; u128 M; } ovr;
static ovr *OV = NULL; static size_t NOV = 0;
static int ov_find(u64 p, u128 *M){ size_t lo = 0, hi = NOV; while(lo < hi){ size_t m = (lo + hi) / 2; if(OV[m].p < p) lo = m + 1; else hi = m; } if(lo < NOV && OV[lo].p == p){ *M = OV[lo].M; return 1; } return 0; }

/* the modulus M(p) used by (W1), per mode */
static u128 get_M(size_t i){
  if(gMode == MODE_PLANT && NOV){ u128 Mo; if(ov_find(Tp[i], &Mo)) return Mo; }
  u64 o = tab_ord(i), z = tab_z(i);
  if(gMode == MODE_FERMAT) return o;
  if(gMode == MODE_LUCAS) return z;
  return (u128)o * (z / gcd64(o, z)); }

/* ================================ pass context (split part) ================================ */
typedef struct {
  mpz_t S, A, B, W, nu; char Ss[512], As[512], Bs[512];
  int S1, W1, Wsmall; u64 Ws, nus;
  u64 *Sl; int nSl; u32 *S32; int nS32;
} passctx;
static passctx PC;
static int PCinit = 0;
static inline u128 S_mod(u128 q){ /* S mod q, 1 <= q < 2^95 */
  if(q == 1) return 0; if(PC.S1) return 1;
  u128 r = 0; for(int i = 0; i < PC.nS32; i++) r = ((r << 32) | PC.S32[i]) % q; return r; }
static int pass_setup(const char *Ss, const char *As, const char *Bs){
  if(!PCinit){ mpz_inits(PC.S, PC.A, PC.B, PC.W, PC.nu, NULL); PCinit = 1; PC.Sl = NULL; PC.S32 = NULL; }
  if(mpz_set_str(PC.S, Ss, 10) || mpz_set_str(PC.A, As, 10) || mpz_set_str(PC.B, Bs, 10)) die("bad S A B");
  if(mpz_sgn(PC.S) <= 0 || mpz_even_p(PC.S) || mpz_sgn(PC.A) <= 0 || mpz_sgn(PC.B) <= 0) die("S must be odd >= 1, A, B >= 1");
  mpz_get_str(PC.Ss, 10, PC.S); mpz_get_str(PC.As, 10, PC.A); mpz_get_str(PC.Bs, 10, PC.B);
  mpz_t g, one, r2; mpz_inits(g, one, r2, NULL);
  mpz_gcd(g, PC.A, PC.B);
  int compat = (mpz_cmp_ui(g, 1) == 0 || mpz_cmp_ui(g, 2) == 0);
  mpz_lcm(PC.W, PC.A, PC.B);
  if(compat){ mpz_set_ui(one, 1); mpz_mod(one, one, PC.A); mpz_sub_ui(r2, PC.B, 1); if(!crt_mpz(PC.nu, g, one, PC.A, r2, PC.B)) compat = 0; }
  mpz_clears(g, one, r2, NULL);
  if(!compat) return 0;
  PC.S1 = (mpz_cmp_ui(PC.S, 1) == 0); PC.W1 = (mpz_cmp_ui(PC.W, 1) == 0);
  PC.Wsmall = (mpz_sizeinbase(PC.W, 2) <= 62); PC.Ws = PC.Wsmall ? mpz_get_ui(PC.W) : 0; PC.nus = PC.Wsmall ? mpz_get_ui(PC.nu) : 0;
  free(PC.Sl); free(PC.S32);
  size_t cnt = 0; PC.nSl = (int)((mpz_sizeinbase(PC.S, 2) + 63) / 64); PC.Sl = calloc(PC.nSl + 1, 8); mpz_export(PC.Sl, &cnt, -1, 8, 0, 0, PC.S);
  PC.nS32 = (int)((mpz_sizeinbase(PC.S, 2) + 31) / 32); PC.S32 = calloc(PC.nS32 + 1, 4); mpz_export(PC.S32, &cnt, 1, 4, 0, 0, PC.S);
  if((int)cnt != PC.nS32) die("S export");
  return 1; }

/* e = p3 * m * S - 1 as little-endian limbs (buffer >= nSl + 3); returns limb count */
static int build_nm1(u64 p3, u128 m, u64 *e){
  u128 pm = (u128)p3 * m;
  if(PC.S1){ u128 t = pm - 1; e[0] = (u64)t; e[1] = (u64)(t >> 64); return e[1] ? 2 : 1; }
  u64 a[2] = { (u64)pm, (u64)(pm >> 64) }; int ns = PC.nSl;
  for(int i = 0; i < ns + 2; i++) e[i] = 0;
  for(int i = 0; i < ns; i++){ u128 carry = 0; for(int j = 0; j < 2; j++){ u128 t = (u128)a[j] * PC.Sl[i] + e[i + j] + carry; e[i + j] = (u64)t; carry = t >> 64; } e[i + 2] = (u64)carry; }
  for(int i = 0; i < ns + 2; i++){ if(e[i]){ e[i]--; break; } e[i] = ~(u64)0; }
  int nl = ns + 2; while(nl > 1 && !e[nl - 1]) nl--; return nl; }
static void make_n(mpz_t n, u64 a, u64 b, u64 c){ mpz_set_ui(n, a); mpz_mul_ui(n, n, b); mpz_mul_ui(n, n, c); if(!PC.S1) mpz_mul(n, n, PC.S); }
static void make_n_m(mpz_t n, u64 p3, u128 m){ mpz_set_u128(n, m); mpz_mul_ui(n, n, p3); if(!PC.S1) mpz_mul(n, n, PC.S); }
/* (W1): M | a*b*S - 1 */
static int w1_holds(u128 M, u64 a, u64 b){ if(M == 1) return 1; u128 x = mulmod128(((u128)a * b) % M, S_mod(M), M); return x == 1; }
/* (W2): a*b*c*S == 1 mod A and == -1 mod B */
static int w2_holds(u64 a, u64 b, u64 c, mpz_t t){
  if(PC.W1) return 1;
  make_n(t, a, b, c);
  if(PC.Wsmall){ u64 r = mpz_fdiv_ui(t, PC.Ws); return r == PC.nus; }
  mpz_mod(t, t, PC.W); return mpz_cmp(t, PC.nu) == 0; }

/* ================================ exact BPSW (GMP) ================================ */
typedef struct { mpz_t d, x, nm1, U, V, Qk, t, two, dl; } bctx;
static void bctx_init(bctx *b){ mpz_inits(b->d, b->x, b->nm1, b->U, b->V, b->Qk, b->t, b->two, b->dl, NULL); mpz_set_ui(b->two, 2); }
static void bctx_clear(bctx *b){ mpz_clears(b->d, b->x, b->nm1, b->U, b->V, b->Qk, b->t, b->two, b->dl, NULL); }
/* strong Fermat base 2, n odd >= 3 */
static int sprp2_mpz(bctx *b, const mpz_t n){
  mpz_sub_ui(b->nm1, n, 1); mp_bitcnt_t s = mpz_scan1(b->nm1, 0); mpz_tdiv_q_2exp(b->d, b->nm1, s);
  mpz_powm(b->x, b->two, b->d, n);
  if(mpz_cmp_ui(b->x, 1) == 0 || mpz_cmp(b->x, b->nm1) == 0) return 1;
  for(mp_bitcnt_t i = 1; i < s; i++){ mpz_mul(b->x, b->x, b->x); mpz_mod(b->x, b->x, n); if(mpz_cmp(b->x, b->nm1) == 0) return 1; }
  return 0; }
static inline void half_mod(mpz_t x, const mpz_t n){ if(mpz_odd_p(x)) mpz_add(x, x, n); mpz_tdiv_q_2exp(x, x, 1); }
/* strong Lucas, P = 1, Q = (1-D)/4, with delta = n - eps given in b->dl */
static int slprp_core(bctx *b, const mpz_t n, long D){
  long Q = (1 - D) / 4;
  mp_bitcnt_t s = mpz_scan1(b->dl, 0); mpz_tdiv_q_2exp(b->d, b->dl, s);
  mpz_set_ui(b->U, 0); mpz_set_ui(b->V, 2); mpz_set_ui(b->Qk, 1);
  for(long i = (long)mpz_sizeinbase(b->d, 2) - 1; i >= 0; i--){
    mpz_mul(b->U, b->U, b->V); mpz_mod(b->U, b->U, n);
    mpz_mul(b->V, b->V, b->V); mpz_submul_ui(b->V, b->Qk, 2); mpz_mod(b->V, b->V, n);
    mpz_mul(b->Qk, b->Qk, b->Qk); mpz_mod(b->Qk, b->Qk, n);
    if(mpz_tstbit(b->d, i)){
      mpz_mul_si(b->x, b->U, D); mpz_add(b->x, b->x, b->V); mpz_mod(b->x, b->x, n); half_mod(b->x, n);
      mpz_add(b->U, b->U, b->V); mpz_mod(b->U, b->U, n); half_mod(b->U, n);
      mpz_swap(b->V, b->x);
      mpz_mul_si(b->Qk, b->Qk, Q); mpz_mod(b->Qk, b->Qk, n);
    }
  }
  if(!mpz_sgn(b->U) || !mpz_sgn(b->V)) return 1;
  for(mp_bitcnt_t i = 1; i < s; i++){
    mpz_mul(b->V, b->V, b->V); mpz_submul_ui(b->V, b->Qk, 2); mpz_mod(b->V, b->V, n);
    mpz_mul(b->Qk, b->Qk, b->Qk); mpz_mod(b->Qk, b->Qk, n);
    if(!mpz_sgn(b->V)) return 1;
  }
  return 0; }
static int selfridge_mpz(const mpz_t n, long *Dout){
  if(mpz_perfect_square_p(n)) return 0;
  long D = 5;
  for(;;){ if(mpz_si_kronecker(D, n) == -1){ *Dout = D; return 1; } D = D > 0 ? -D - 2 : -D + 2; } }
/* exact strong BPSW as src/bpsw.py: n > 1, odd, strong Fermat 2, Selfridge A strong Lucas */
static int bpsw_mpz(bctx *b, const mpz_t n){
  if(mpz_cmp_ui(n, 1) <= 0 || mpz_even_p(n)) return 0;
  if(!sprp2_mpz(b, n)) return 0;
  long D; if(!selfridge_mpz(n, &D)) return 0;
  mpz_add_ui(b->dl, n, 1); return slprp_core(b, n, D); }
/* strong Lucas with FIXED D (lucas control mode): delta = n - (D/n) */
static int slprp_fixed_mpz(bctx *b, const mpz_t n, long D){
  if(mpz_cmp_ui(n, 1) <= 0 || mpz_even_p(n)) return 0;
  int eps = mpz_si_kronecker(D, n); if(!eps) return 0;
  if(eps < 0) mpz_add_ui(b->dl, n, 1); else mpz_sub_ui(b->dl, n, 1);
  return slprp_core(b, n, D); }

/* ================================ u64 primality / factoring (control modes) ================================ */
static int is_prime_u64(u64 n){ /* n < 2^63 */
  if(n < 2) return 0;
  static const u32 sp[] = {2,3,5,7,11,13,17,19,23,29,31,37};
  for(int i = 0; i < 12; i++){ if(n == sp[i]) return 1; if(n % sp[i] == 0) return 0; }
  if(n < 1369) return 1;
  m64 c; m64_init(&c, n); u64 d = n - 1; int s = __builtin_ctzll(d); d >>= s;
  u64 one = c.r1, mone = c.n - c.r1;
  static const u64 bases[] = {2, 325, 9375, 28178, 450775, 9780504, 1795265022};
  for(int i = 0; i < 7; i++){
    u64 a = bases[i] % n; if(!a) continue;
    u64 x = m64_pow(m64_to(a, &c), d, &c);
    if(x == one || x == mone) continue;
    int ok = 0; for(int r = 1; r < s; r++){ x = m64_mul(x, x, &c); if(x == mone){ ok = 1; break; } }
    if(!ok) return 0;
  }
  return 1; }
static u64 rho_u64(u64 n){ /* nontrivial factor of odd composite n < 2^63, 0 on failure */
  m64 c; m64_init(&c, n);
  for(u64 c0 = 1; c0 < 2000; c0++){
    u64 cm = m64_to(c0, &c), y = m64_to(c0 + 2, &c), x = y, ys = y, q = c.r1, g = 1, r = 1;
    while(g == 1){
      x = y; for(u64 i = 0; i < r; i++) y = m64_add(m64_mul(y, y, &c), cm, &c);
      for(u64 k = 0; k < r && g == 1; k += 128){
        ys = y; u64 lim = (r - k < 128) ? r - k : 128;
        for(u64 i = 0; i < lim; i++){ y = m64_add(m64_mul(y, y, &c), cm, &c); q = m64_mul(q, m64_sub(x, y, &c), &c); }
        g = gcd64(q, n);
      }
      r <<= 1; if(r > ((u64)1 << 36)) break;
    }
    if(g == n || g == 0){ g = 1; for(u64 it = 0; it < ((u64)1 << 36) && g == 1; it++){ ys = m64_add(m64_mul(ys, ys, &c), cm, &c); g = gcd64(m64_sub(x, ys, &c), n); } }
    if(g != n && g > 1) return g;
  }
  return 0; }
static const u32 SMALLP[] = {3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97};
/* m = a*b with a < b both inert table primes?  1 yes (sets ia, ib), 0 no, dies on factoring failure */
static int split_m(u64 m, int64_t *ia, int64_t *ib){
  if(is_prime_u64(m)) return 0;
  u64 f = 0;
  for(size_t i = 0; i < sizeof(SMALLP) / sizeof(SMALLP[0]); i++) if(m % SMALLP[i] == 0){ f = SMALLP[i]; break; }
  if(!f) f = rho_u64(m);
  if(!f || f == 1 || f == m || m % f) die("split_m: factoring failure (completeness would be at risk)");
  u64 g = m / f, x = f < g ? f : g, y = f < g ? g : f;
  if(x == y) return 0;
  int64_t i1 = tab_find(x); if(i1 < 0) return 0;
  int64_t i2 = tab_find(y); if(i2 < 0) return 0;
  *ia = i1; *ib = i2; return 1; }

/* factorisation for printing (mpz) */
static void mpz_rho_factor(const mpz_t n, char *out, size_t cap){
  if(mpz_cmp_ui(n, 1) == 0) return;
  if(mpz_probab_prime_p(n, 30)){ size_t l = strlen(out); if(l) { strncat(out, " * ", cap - l - 1); } l = strlen(out); mpz_get_str(out + l, 10, n); return; }
  mpz_t x, y, g, c, t; mpz_inits(x, y, g, c, t, NULL);
  for(unsigned long cc = 1;; cc++){
    mpz_set_ui(x, 2); mpz_set_ui(y, 2); mpz_set_ui(g, 1); mpz_set_ui(c, cc);
    while(mpz_cmp_ui(g, 1) == 0){
      mpz_mul(x, x, x); mpz_add(x, x, c); mpz_mod(x, x, n);
      mpz_mul(y, y, y); mpz_add(y, y, c); mpz_mod(y, y, n); mpz_mul(y, y, y); mpz_add(y, y, c); mpz_mod(y, y, n);
      mpz_sub(t, x, y); mpz_abs(t, t); mpz_gcd(g, t, n);
    }
    if(mpz_cmp(g, n) != 0) break;
  }
  mpz_divexact(t, n, g); mpz_rho_factor(g, out, cap); mpz_rho_factor(t, out, cap);
  mpz_clears(x, y, g, c, t, NULL); }
static void factor_string(const mpz_t n0, char *out, size_t cap){
  out[0] = 0; mpz_t n; mpz_init_set(n, n0);
  for(unsigned long p = 2; p < 100000 && mpz_cmp_ui(n, 1) > 0; p++){
    while(mpz_divisible_ui_p(n, p)){ size_t l = strlen(out); snprintf(out + l, cap - l, "%s%lu", l ? " * " : "", p); mpz_divexact_ui(n, n, p); }
  }
  if(mpz_cmp_ui(n, 1) > 0) mpz_rho_factor(n, out, cap);
  mpz_clear(n); }

/* ================================ table construction ================================ */
#define MAXF 16
typedef struct { u8 *comp, *need, *nf; u64 *rem; u32 *fac; u64 *op; u32 *o1, *o2; size_t nout; size_t *ovp; u64 *ov1, *ov2; size_t nov, capov; } segws;
static void init_base_primes(u64 X){
  u64 sq = (u64)sqrtl((long double)(X + 1)) + 2; while(sq * sq < X + 2) sq++;
  if(BP && BPmax >= sq) return;
  free(BP); u8 *s = calloc(sq + 1, 1); BP = malloc(sizeof(u32) * (sq / 2 + 16)); NBP = 0;
  for(u64 i = 2; i <= sq; i++){ if(!s[i]){ BP[NBP++] = (u32)i; for(u64 j = i * i; j <= sq; j += i) s[j] = 1; } }
  free(s); BPmax = sq; }
static void segws_alloc(segws *w){
  u64 L = gSeg + 4;
  w->comp = malloc(L); w->need = malloc(L); w->nf = malloc(L); w->rem = malloc(8 * L); w->fac = malloc(4 * MAXF * L);
  w->op = malloc(8 * (L / 2 + 4)); w->o1 = malloc(4 * (L / 2 + 4)); w->o2 = malloc(4 * (L / 2 + 4));
  w->capov = 16; w->ovp = malloc(sizeof(size_t) * 16); w->ov1 = malloc(8 * 16); w->ov2 = malloc(8 * 16); w->nov = 0; }
static void segws_free(segws *w){ free(w->comp); free(w->need); free(w->nf); free(w->rem); free(w->fac); free(w->op); free(w->o1); free(w->o2); free(w->ovp); free(w->ov1); free(w->ov2); }
/* all inert primes p in [lo, hi) (lo >= 3), with exact ord_p(2) and z_D(p). Requires base primes >= sqrt(hi+1). */
static void seg_process(u64 lo, u64 hi, segws *w){
  u64 len = hi - lo, base = lo - 1;   /* factor-sieve index t <-> v = base + t, t in [0, len+1] */
  w->nout = 0; w->nov = 0;
  memset(w->comp, 0, len);
  for(size_t i = 0; i < NBP; i++){ u64 r = BP[i]; if(r * r >= hi) break; u64 st = ((lo + r - 1) / r) * r; if(st < r * r) st = r * r; for(u64 j = st; j < hi; j += r) w->comp[j - lo] = 1; }
  memset(w->need, 0, len + 2);
  u64 aQ = (u64)(gQ < 0 ? -gQ : gQ);
  for(u64 t = 0; t < len; t++){
    u64 p = lo + t;
    if(w->comp[t] || !(p & 1)){ w->comp[t] = 1; continue; }
    if(jacobi_s64(gD, p) != -1 || (aQ % p) == 0){ w->comp[t] = 1; continue; }   /* [T-INERT] keep (D/p) = -1, p not | Q */
    w->need[t] = 1; w->need[t + 2] = 1;
  }
  for(u64 t = 0; t <= len + 1; t++) if(w->need[t]){ w->rem[t] = base + t; w->nf[t] = 0; }
  for(size_t i = 0; i < NBP; i++){
    u64 r = BP[i]; if(r > hi) break;
    u64 st = ((base + r - 1) / r) * r;
    for(u64 v = st; v <= hi; v += r){
      u64 t = v - base; if(!w->need[t]) continue;
      u64 x = w->rem[t]; u32 e = 0; do { x /= r; e++; } while(x % r == 0);
      w->rem[t] = x; if(w->nf[t] >= MAXF) die("MAXF"); w->fac[t * MAXF + w->nf[t]++] = (u32)(r << 6) | e;
    }
  }
  u64 Dabs = (u64)(gD < 0 ? -gD : gD);
  for(u64 t = 0; t < len; t++){
    if(w->comp[t]) continue;
    u64 p = lo + t; m64 c; m64_init(&c, p);
    /* ord_p(2): strip prime factors of p - 1 */
    u64 o = p - 1; u64 tt = t;
    for(int k = 0; k <= w->nf[tt]; k++){
      u64 q; u32 e;
      if(k < w->nf[tt]){ q = w->fac[tt * MAXF + k] >> 6; e = w->fac[tt * MAXF + k] & 63; } else { if(w->rem[tt] <= 1) break; q = w->rem[tt]; e = 1; }
      for(u32 j = 0; j < e; j++){ if(m64_pow2_is1_c(o / q, &c)) o /= q; else break; }   /* [T-ORD] */
    }
    /* z_D(p): strip prime factors of p + 1 */
    u64 Dm = m64_to(gD < 0 ? p - (Dabs % p) : Dabs % p, &c), Qm = m64_to(s64mod(gQ, p), &c);
    u64 z = p + 1; tt = t + 2;
    if(!m64_lucasU_is0(z, Dm, Qm, &c)) die("U_(p+1) != 0 mod inert p");
    for(int k = 0; k <= w->nf[tt]; k++){
      u64 q; u32 e;
      if(k < w->nf[tt]){ q = w->fac[tt * MAXF + k] >> 6; e = w->fac[tt * MAXF + k] & 63; } else { if(w->rem[tt] <= 1) break; q = w->rem[tt]; e = 1; }
      for(u32 j = 0; j < e; j++){ if(m64_lucasU_is0(z / q, Dm, Qm, &c)) z /= q; else break; }   /* [T-Z] */
    }
    u64 i1 = (p - 1) / o, i2 = (p + 1) / z;
    w->op[w->nout] = p;
    if(i1 >> 32 || i2 >> 32){
      if(w->nov == w->capov){ w->capov *= 2; w->ovp = realloc(w->ovp, sizeof(size_t) * w->capov); w->ov1 = realloc(w->ov1, 8 * w->capov); w->ov2 = realloc(w->ov2, 8 * w->capov); }
      w->ovp[w->nov] = w->nout; w->ov1[w->nov] = i1; w->ov2[w->nov] = i2; w->nov++;
      w->o1[w->nout] = 0; w->o2[w->nout] = 0;
    } else { w->o1[w->nout] = (u32)i1; w->o2[w->nout] = (u32)i2; }
    w->nout++;
  }
}
static void build_table(u64 X){
  free(Tp); free(Ti1); free(Ti2); free(Tov); Tp = NULL; Ti1 = Ti2 = NULL; Tov = NULL; NT = 0; NTov = 0; CapTov = 0;
  init_base_primes(X);
  TCap = X < 1000 ? X / 2 + 64 : (size_t)(0.65 * (double)X / (log((double)X) - 1.1)) + 4096;
  Tp = malloc(8 * TCap); Ti1 = malloc(4 * TCap); Ti2 = malloc(4 * TCap);
  if(X < 3) return;
  u64 lo0 = gX0 < 3 ? 3 : (gX0 | 1);
  if(lo0 > X) return;
  u64 nseg = (X - lo0) / gSeg + 1;
  int R = 2 * gThreads; segws *ws = malloc(sizeof(segws) * R); for(int i = 0; i < R; i++) segws_alloc(&ws[i]);
  for(u64 s0 = 0; s0 < nseg; s0 += R){
    u64 s1 = s0 + R < nseg ? s0 + R : nseg;
    #pragma omp parallel for schedule(dynamic, 1) num_threads(gThreads)
    for(u64 s = s0; s < s1; s++){ u64 lo = lo0 + s * gSeg, hi = lo + gSeg; if(hi > X + 1) hi = X + 1; seg_process(lo, hi, &ws[s - s0]); }
    for(u64 s = s0; s < s1; s++){
      segws *w = &ws[s - s0];
      if(NT + w->nout > TCap){ TCap = (NT + w->nout) * 5 / 4 + 4096; Tp = realloc(Tp, 8 * TCap); Ti1 = realloc(Ti1, 4 * TCap); Ti2 = realloc(Ti2, 4 * TCap); if(!Tp || !Ti1 || !Ti2) die("table realloc"); }
      memcpy(Tp + NT, w->op, 8 * w->nout); memcpy(Ti1 + NT, w->o1, 4 * w->nout); memcpy(Ti2 + NT, w->o2, 4 * w->nout);
      for(size_t k = 0; k < w->nov; k++){
        if(NTov == CapTov){ CapTov = CapTov ? 2 * CapTov : 64; Tov = realloc(Tov, sizeof(ovent) * CapTov); }
        Tov[NTov].idx = NT + w->ovp[k]; Tov[NTov].i1 = w->ov1[k]; Tov[NTov].i2 = w->ov2[k]; NTov++;
      }
      NT += w->nout;
    }
  }
  for(int i = 0; i < R; i++) segws_free(&ws[i]); free(ws);
}

/* ================================ small-prime sieve for the Fermat-mod-m prefilter ================================ */
/* Lemire fastmod: a mod d for a, d < 2^32, M = floor((2^64-1)/d) + 1 */
static inline u32 fastmod32(u32 a, u64 M, u32 d){ u64 low = M * a; return (u32)(((u128)low * d) >> 64); }
typedef struct { u32 l, o; u64 Ml, Mo; u32 t32l, t64l, t32o, t64o; u16 *invl, *invo; } svprime;
static svprime *SV = NULL; static int NSV = 0; static u32 *PSo = NULL;   /* PSo[q] = S mod o_q for the current pass */
/* x mod d for x < 2^96, d < 2^12 (all partial sums < 2^26) */
static inline u32 modsmall(u128 x, u64 M, u32 d, u32 t32, u32 t64){
  u32 s = fastmod32((u32)(x >> 64), M, d) * t64 + fastmod32((u32)(x >> 32), M, d) * t32 + fastmod32((u32)x, M, d);
  return fastmod32(s, M, d); }
static void init_sieve_primes(u32 B){
  if(B > 4096) B = 4096;     /* modsmall / u16 tables assume l <= 4096 */
  for(int q = 0; q < NSV; q++){ free(SV[q].invl); free(SV[q].invo); } free(SV); free(PSo);
  SV = calloc(B / 2 + 8, sizeof(svprime)); PSo = calloc(B / 2 + 8, sizeof(u32)); NSV = 0;
  u8 *s = calloc(B + 1, 1);
  for(u32 i = 2; i <= B; i++){ if(s[i]) continue; for(u64 j = (u64)i * i; j <= B; j += i) s[j] = 1;
    if(i == 2) continue;
    u32 o = 1; u64 x = 2 % i; while(x != 1){ x = x * 2 % i; o++; }
    svprime *P = &SV[NSV]; P->l = i; P->o = o; P->Ml = UINT64_MAX / i + 1; P->Mo = UINT64_MAX / o + 1;
    P->t32l = (u32)((((u64)1) << 32) % i); P->t64l = (u32)((((u128)1) << 64) % i); P->t32o = (u32)((((u64)1) << 32) % o); P->t64o = (u32)((((u128)1) << 64) % o);
    P->invl = malloc(sizeof(u16) * i); P->invo = malloc(sizeof(u16) * o);
    for(u32 a = 0; a < i; a++) P->invl[a] = a ? (u16)inv_u32(a, i) : 0xFFFF;
    for(u32 a = 0; a < o; a++) P->invo[a] = (o == 1) ? 0 : (gcd64(a, o) == 1 ? (u16)inv_u32(a, o) : 0xFFFF);
    PSo[NSV] = 1 % o; NSV++; }
  free(s); }
static void sieve_pass_init(void){ for(int q = 0; q < NSV; q++) PSo[q] = (u32)S_mod(SV[q].o); }
/* mark dead[d] (d < len) for m_d = m + d*step when some small l | m_d and ord_l(2) does not divide n_d - 1 */
static u64 sieveA(u64 p3, u128 mstart, u128 step, u64 len, u8 *dead){
  memset(dead, 0, len);
  if((mstart >> 90) || (step >> 90)) return 0;        /* outside modsmall range: no marks (always safe) */
  int wide = ((mstart + (u128)(len - 1) * step) >> 63) != 0;   /* m128 prefilter is ~5x dearer: sieve deeper */
  u64 lim = (wide ? 30 : 6) * len + 64; if(lim > gSieveB) lim = gSieveB;
  for(int q = 0; q < NSV; q++){
    const svprime *P = &SV[q]; u32 l = P->l; if(l > lim) break;
    u32 o = P->o;
    u32 Lm = modsmall(step, P->Ml, l, P->t32l, P->t64l); if(!Lm) continue;      /* l | step: no marks (always safe) */
    u32 mk = modsmall(mstart, P->Ml, l, P->t32l, P->t64l);
    u32 d0 = fastmod32((l - mk) * (u32)P->invl[Lm], P->Ml, l);                  /* first d with l | m_d */
    if(d0 >= len) continue;
    u32 so = modsmall(step, P->Mo, o, P->t32o, P->t64o), mo = modsmall(mstart, P->Mo, o, P->t32o, P->t64o);
    u32 x = fastmod32(mo + fastmod32(d0, P->Mo, o) * so, P->Mo, o);              /* m_d0 mod o */
    u32 dx = fastmod32(fastmod32(l, P->Mo, o) * so, P->Mo, o);                   /* (l * step) mod o */
    u32 ps = fastmod32(modsmall((u128)p3, P->Mo, o, P->t32o, P->t64o) * PSo[q], P->Mo, o);
    u32 t = P->invo[ps];   /* x == t  <=>  ord_l | p3 S m - 1 ; t = 0xFFFF (never matched) if p3 S not invertible mod o */
    for(u64 d = d0; d < len; d += l){ if(x != t) dead[d] = 1; x += dx; if(x >= o) x -= o; }
  }
  u64 cnt = 0; for(u64 d = 0; d < len; d++) cnt += dead[d]; return cnt; }

/* ================================ plants ================================ */
typedef struct { u64 p1, p2, p3; char S[512], A[512], B[512]; int pass; int valid; int found; char why[64]; } plant;
static plant *PL = NULL; static size_t NPL = 0;
typedef struct { u64 p3; u128 m; size_t idx; } pkey;
static pkey *PK = NULL; static size_t NPK = 0;
static int pkey_cmp(const void *a, const void *b){ const pkey *x = a, *y = b; if(x->p3 != y->p3) return x->p3 < y->p3 ? -1 : 1; if(x->m != y->m) return x->m < y->m ? -1 : 1; return 0; }
static int ovr_cmp(const void *a, const void *b){ const ovr *x = a, *y = b; return x->p < y->p ? -1 : x->p > y->p; }
static int64_t pkey_find(u64 p3, u128 m){ size_t lo = 0, hi = NPK; pkey k = { p3, m, 0 };
  while(lo < hi){ size_t mid = (lo + hi) / 2; if(pkey_cmp(&PK[mid], &k) < 0) lo = mid + 1; else hi = mid; }
  return (lo < NPK && PK[lo].p3 == p3 && PK[lo].m == m) ? (int64_t)PK[lo].idx : -1; }
static void load_plants(const char *fn){
  FILE *f = fopen(fn, "r"); if(!f) die("cannot open plant file");
  char line[2048]; size_t capP = 0, capO = 0;
  mpz_t a, b, c; mpz_inits(a, b, c, NULL);
  while(fgets(line, sizeof line, f)){
    char *s = line; while(*s == ' ' || *s == '\t') s++;
    if(*s == 'O'){ char ps[64], Ms[256]; if(sscanf(s + 1, "%63s %255s", ps, Ms) != 2) die("bad O line");
      if(NOV == capO){ capO = capO ? 2 * capO : 64; OV = realloc(OV, sizeof(ovr) * capO); }
      mpz_set_str(a, ps, 10); mpz_set_str(b, Ms, 10); u128 M; if(!mpz_get_u128(&M, b) || mpz_sgn(b) <= 0) die("bad override M");
      u64 p = mpz_get_ui(a); mpz_mul(c, a, a); mpz_sub_ui(c, c, 1); if(!mpz_divisible_p(c, b)) die("override M does not divide p^2 - 1");
      OV[NOV].p = p; OV[NOV].M = M; NOV++; }
    else if(*s == 'T'){ plant P; memset(&P, 0, sizeof P); unsigned long long x1, x2, x3;
      if(sscanf(s + 1, "%llu %llu %llu %511s %511s %511s", &x1, &x2, &x3, P.S, P.A, P.B) != 6) die("bad T line");
      P.p1 = x1; P.p2 = x2; P.p3 = x3;
      mpz_set_str(a, P.S, 10); mpz_get_str(P.S, 10, a); mpz_set_str(a, P.A, 10); mpz_get_str(P.A, 10, a); mpz_set_str(a, P.B, 10); mpz_get_str(P.B, 10, a);
      if(NPL == capP){ capP = capP ? 2 * capP : 64; PL = realloc(PL, sizeof(plant) * capP); }
      PL[NPL++] = P; }
  }
  fclose(f); mpz_clears(a, b, c, NULL);
  qsort(OV, NOV, sizeof(ovr), ovr_cmp);
  for(size_t i = 1; i < NOV; i++) if(OV[i].p == OV[i - 1].p && OV[i].M != OV[i - 1].M) die("conflicting overrides");
}

/* ================================ counters / thread state / hits ================================ */
typedef struct {
  u64 Aprimes, Acands, Asieved, Apre, Ashape, Afinal, Aunsolv, Cprimes, Cpairs, Cnoinv, Cincompat, Cinner, Cp1, Cw1, Cfinal, Cunsolv;
  u64 hits, longpairs, Abig;
  double tA, tC, cpu; char pad[64];
} ctrs;
typedef struct { ctrs c; bctx b; mpz_t n, t1; u8 *dead; u64 *e; } tstate;
typedef struct { u64 p1, p2, p3; u128 m; } hitrec;
static hitrec *HITS = NULL; static size_t NHITS = 0, CapHITS = 0;
static void record_hit(const char *tag, tstate *ts, u64 p1, u64 p2, u64 p3, u128 m, int branch){
  char buf[4096], mb[48];
  #pragma omp critical(hits)
  {
    if(NHITS == CapHITS){ CapHITS = CapHITS ? 2 * CapHITS : 256; HITS = realloc(HITS, sizeof(hitrec) * CapHITS); }
    memset(&HITS[NHITS], 0, sizeof(hitrec)); HITS[NHITS].p1 = p1; HITS[NHITS].p2 = p2; HITS[NHITS].p3 = p3; HITS[NHITS].m = m; NHITS++;
    if(!gQuiet || gMode == MODE_BPSW){   /* FIX F1: a bpsw HIT (counterexample) is never silenced by --quiet */
      char *ns = mpz_get_str(NULL, 10, ts->n);
      if(p1){ printf("%s n = %s = %llu * %llu * %llu%s%s  [branch %c]\n", tag, ns, (unsigned long long)p1, (unsigned long long)p2, (unsigned long long)p3, PC.S1 ? "" : " * ", PC.S1 ? "" : PC.Ss, branch); }
      else { mpz_t mm; mpz_init(mm); mpz_set_u128(mm, m); factor_string(mm, buf, sizeof buf); mpz_clear(mm);
        printf("%s n = %s = (m = %s = %s) * %llu%s%s  [branch %c]\n", tag, ns, u128str(m, mb), buf, (unsigned long long)p3, PC.S1 ? "" : " * ", PC.S1 ? "" : PC.Ss, branch);
        mpz_t nn; mpz_init_set(nn, ts->n); factor_string(nn, buf, sizeof buf); printf("     full factorisation: %s\n", buf); mpz_clear(nn); }
      fflush(stdout); free(ns);
    }
  }
  ts->c.hits++; }

/* final predicate for a triple of table indices (i1 < i2 < i3), per mode */
static void final_triple(tstate *ts, int64_t i1, int64_t i2, size_t i3, int branch){
  u64 p1 = Tp[i1], p2 = Tp[i2], p3 = Tp[i3];
  if(branch == 'A') ts->c.Afinal++; else ts->c.Cfinal++;
  switch(gMode){
    case MODE_BPSW: make_n(ts->n, p1, p2, p3); if(bpsw_mpz(&ts->b, ts->n)) record_hit("HIT", ts, p1, p2, p3, 0, branch); break;
    case MODE_WEAK:
      if(w1_holds(get_M(i1), p2, p3) && w1_holds(get_M(i2), p1, p3) && w1_holds(get_M(i3), p1, p2) && w2_holds(p1, p2, p3, ts->t1)){ make_n(ts->n, p1, p2, p3); record_hit("WEAK", ts, p1, p2, p3, 0, branch); }
      break;
    case MODE_PLANT: { int64_t k = pkey_find(p3, (u128)p1 * p2); if(k >= 0){ make_n(ts->n, p1, p2, p3);
        #pragma omp atomic
        PL[k].found++;
        record_hit("PLANTFOUND", ts, p1, p2, p3, 0, branch); } } break;
    case MODE_FERMAT: make_n(ts->n, p1, p2, p3); if(sprp2_mpz(&ts->b, ts->n)) record_hit("FERMAT", ts, p1, p2, p3, 0, branch); break;
    case MODE_LUCAS: make_n(ts->n, p1, p2, p3); if(slprp_fixed_mpz(&ts->b, ts->n, (long)gD)) record_hit("LUCAS", ts, p1, p2, p3, 0, branch); break;
  }
}

/* ================================ branch A ================================ */
typedef struct { u128 m0, step; u64 count; } progA;
/* the class of m for anchor p3; 0 = unsolvable */
static int setupA(size_t i3, progA *pa){
  u64 p3 = Tp[i3]; u128 M3 = get_M(i3), bound = (u128)p3 * p3, r, L;
  pa->count = 0; pa->m0 = 0; pa->step = 1;
  if(!invmod128(S_mod(M3), M3, &r)) return 0;                                   /* [JA1] */
  L = M3;
  if(!PC.W1){
    int st = 2;
    if(PC.Wsmall){
      u128 W = PC.Ws, ti; u128 t = mulmod128(p3 % W, S_mod(W), W);
      if(!invmod128(t, W, &ti)) return 0;                                       /* [JA2] */
      st = crt128(r, L, mulmod128(ti, PC.nus, W), W, &r, &L);
      if(st == 0) return 0;                                                     /* [JA3] */
    }
    if(st == 2){   /* big moduli: GMP */
      mpz_t t, ti, aw, rz, Lz, r1, L1; mpz_inits(t, ti, aw, rz, Lz, r1, L1, NULL); int ok = 1;
      mpz_set_ui(t, p3); mpz_mul(t, t, PC.S); mpz_mod(t, t, PC.W);
      if(!mpz_invert(ti, t, PC.W)) ok = 0;                                      /* [JA2] */
      if(ok){ mpz_mul(aw, ti, PC.nu); mpz_mod(aw, aw, PC.W); mpz_set_u128(r1, r); mpz_set_u128(L1, M3);
        if(!crt_mpz(rz, Lz, r1, L1, aw, PC.W)) ok = 0; }                        /* [JA3] */
      if(ok){
        if(mpz_sizeinbase(Lz, 2) <= 126){ mpz_get_u128(&r, rz); mpz_get_u128(&L, Lz); }
        else { u128 rr; int small = mpz_get_u128(&rr, rz);   /* modulus > 2^126 > bound: at most m = rz */
          pa->step = 0; pa->m0 = rr; pa->count = (small && rr >= 15 && rr < bound && (rr & 1)) ? 1 : 0; ok = 2; }
      }
      mpz_clears(t, ti, aw, rz, Lz, r1, L1, NULL);
      if(ok == 0) return 0;                                                     /* [JA2]/[JA3] (GMP path) */
      if(ok == 2) return 1;
    }
  }
  u128 m = r;
  if(m < 15) m += ((15 - m + L - 1) / L) * L;                                    /* [A-M15] */
  u128 step = L;
  if(L & 1){ if(!(m & 1)) m += L; step = 2 * L; }                                /* [A-ODD] */
  else if(!(m & 1)) return 0;                                                   /* [JA4] */
  pa->m0 = m; pa->step = step;
  if(m < bound){ u128 cnt = (bound - 1 - m) / step + 1; if(cnt >> 63) die("count overflow"); pa->count = (u64)cnt; }  /* [A-BOUND] */
  else pa->count = 0;
  return 1; }

/* branch decision: 0 skip (unsolvable), 1 = A, 2 = C */
static int decide(size_t i3, progA *pa){
  int ok = setupA(i3, pa);
  if(gCap == 0) return 2;                                                       /* [P-DECIDE] */
  if(!ok) return 0;                                                             /* [JP0] */
  if(gCap == -1) return 1;
  if(gCap > 0) return pa->count < (u64)gCap ? 1 : 2;
  double p3 = (double)Tp[i3], fw = 1.0 - 9.223372036854775808e18 / (p3 * p3); if(fw < 0) fw = 0;   /* share of m >= 2^63 */
  return (double)pa->count * (gCtest * (1 - fw) + gCtestW * fw) <= (double)i3 * gCpair ? 1 : 2; }

typedef struct { u8 type; size_t i3; u128 m0, step; u64 k0, k1; } item;

static void runA(tstate *ts, const item *it){
  size_t i3 = it->i3; u64 p3 = Tp[i3]; u64 len = it->k1 - it->k0;
  u128 step = it->step, m = it->m0 + (u128)it->k0 * step;
  ts->c.Acands += len;
  int usesieve = gSieve && (gMode == MODE_BPSW || gMode == MODE_FERMAT) && len >= 4 && step;
  if(usesieve) ts->c.Asieved += sieveA(p3, m, step, len, ts->dead);
  int lanes = PC.S1 && (gMode == MODE_BPSW || gMode == MODE_FERMAT);
  u128 qm[2][4]; u64 qe0[2][4], qe1[2][4]; int qn[2] = {0, 0}, qres[4];   /* lane queues: [0] m < 2^63, [1] m >= 2^63 */
  for(u64 d = 0; d <= len; d++, m += step){
    if(lanes){   /* queue sieve survivors; flush full queues (and everything at the end) through the 4-lane prefilter */
      if(d < len && !(usesieve && ts->dead[d])){                                 /* [JA5] dead -> not queued */
        int w = (m >> 63) ? 1 : 0; u128 t = (u128)p3 * m - 1; qm[w][qn[w]] = m; qe0[w][qn[w]] = (u64)t; qe1[w][qn[w]] = (u64)(t >> 64); qn[w]++; }
      for(int w = 0; w < 2; w++){
        if(!qn[w] || (qn[w] < 4 && d < len)) continue;                           /* not full yet */
        int k = qn[w]; for(int i = k; i < 4; i++){ qm[w][i] = qm[w][0]; qe0[w][i] = qe0[w][0]; qe1[w][i] = qe1[w][0]; }
        if(w == 0){ u64 nv[4]; for(int i = 0; i < 4; i++) nv[i] = (u64)qm[0][i]; m64_pow2_is1_x4(nv, qe0[0], qe1[0], qres); }
        else m128_pow2_is1_x4(qm[1], qe0[1], qe1[1], qres);
        qn[w] = 0;
        for(int i = 0; i < k; i++){
          if(!qres[i]) continue;                                                 /* [JA6] */
          ts->c.Apre++; u128 mm = qm[w][i];
          if(gMode == MODE_BPSW){ ts->c.Afinal++; make_n_m(ts->n, p3, mm); if(bpsw_mpz(&ts->b, ts->n)) record_hit("HIT", ts, 0, 0, p3, mm, 'A'); }
          else { int64_t ia, ib; if(!split_m((u64)mm, &ia, &ib) || Tp[ib] >= p3) continue; ts->c.Ashape++; final_triple(ts, ia, ib, i3, 'A'); }   /* [JA7] */
        }
      }
      continue;
    }
    if(d == len) break;
    if(usesieve && ts->dead[d]) continue;                                        /* [JA5] */
    if(gMode == MODE_BPSW || gMode == MODE_FERMAT){
      int nl = build_nm1(p3, m, ts->e);
      if(!fermat_m(m, ts->e, nl)) continue;                                      /* [JA6] */
      ts->c.Apre++;
      if(gMode == MODE_BPSW){ ts->c.Afinal++; make_n_m(ts->n, p3, m); if(bpsw_mpz(&ts->b, ts->n)) record_hit("HIT", ts, 0, 0, p3, m, 'A'); continue; }
    }
    if(gMode == MODE_PLANT){ ts->c.Afinal++; int64_t k = pkey_find(p3, m);
      if(k >= 0){
        #pragma omp atomic
        PL[k].found++;
        make_n(ts->n, PL[k].p1, PL[k].p2, p3); record_hit("PLANTFOUND", ts, PL[k].p1, PL[k].p2, p3, 0, 'A'); }
      continue; }
    int64_t ia, ib;
    if(!split_m((u64)m, &ia, &ib) || Tp[ib] >= p3) continue;                     /* [JA7] */
    ts->c.Ashape++;
    final_triple(ts, ia, ib, i3, 'A');
  }
}

/* ================================ branch C ================================ */
static void runC(tstate *ts, const item *it){
  size_t i3 = it->i3; u64 p3 = Tp[i3]; u128 M3 = get_M(i3), s3, c3 = 0;
  if(!invmod128(S_mod(M3), M3, &s3)){ if(it->k0 == 0) ts->c.Cunsolv++; return; }                     /* [JC1] */
  if(!PC.W1 && PC.Wsmall){ u128 W = PC.Ws, t = mulmod128(p3 % W, S_mod(W), W), ti;
    if(!invmod128(t, W, &ti)){ if(it->k0 == 0) ts->c.Cunsolv++; return; }                           /* [JC2] */
    c3 = mulmod128(ti, PC.nus, W); }
  for(size_t j = it->k0; j < it->k1; j++){
    u64 p2 = Tp[j]; ts->c.Cpairs++;
    u128 r, L = M3, inv, M2 = 0;
    if(!invmod128(p2, M3, &inv)){ ts->c.Cnoinv++; continue; }                                     /* [JC3] */
    r = PC.S1 ? inv : mulmod128(inv, s3, M3);
    int haveM2 = 0, haveW = PC.W1;
    if(L < p2){   /* combine (W1) at p2 */
      M2 = get_M(j); u128 a2, t = mulmod128(p3 % M2, S_mod(M2), M2);
      if(!invmod128(t, M2, &a2)){ ts->c.Cnoinv++; continue; }                                     /* [JC4] */
      int st = crt128(r, L, a2, M2, &r, &L);
      if(st == 0){ ts->c.Cincompat++; continue; }                                                 /* [JC5] */
      haveM2 = (st == 1);
    }
    if(L < p2 && !haveW && PC.Wsmall){   /* combine (W2) */
      u128 W = PC.Ws, ip;
      if(!invmod128(p2 % W, W, &ip)){ ts->c.Cnoinv++; continue; }                                 /* [JC6] */
      int st = crt128(r, L, mulmod128(c3, ip, W), W, &r, &L);
      if(st == 0){ ts->c.Cincompat++; continue; }                                                 /* [JC7] */
      haveW = (st == 1);
    }
    u128 p1 = r;
    if(p1 < 3) p1 += ((3 - p1 + L - 1) / L) * L;                                                  /* [C-P1GE3] */
    u64 cnt = 0;
    for(; p1 < p2; p1 += L){                                                                      /* [C-P1LTP2] */
      cnt++; u64 q1 = (u64)p1;
      if(!haveM2){ if(!M2) M2 = get_M(j); if(!w1_holds(M2, q1, p3)) continue; }                   /* [C-DEF2] */
      if(!haveW){ if(!w2_holds(q1, p2, p3, ts->t1)) continue; }                                   /* [C-DEFW] */
      int64_t i1 = tab_find(q1);
      if(i1 < 0) continue;                                                                        /* [C-P1TAB] */
      ts->c.Cp1++;
      if(!w1_holds(get_M(i1), p2, p3)) continue;                                                  /* [C-W1P1] */
      ts->c.Cw1++;
      final_triple(ts, i1, (int64_t)j, i3, 'C');
    }
    ts->c.Cinner += cnt;
    if(cnt > 1000000){ ts->c.longpairs++;
      #pragma omp critical(log)
      { printf("LONGPAIR p2 = %llu p3 = %llu inner = %llu\n", (unsigned long long)p2, (unsigned long long)p3, (unsigned long long)cnt); fflush(stdout); } }
  }
}

/* ================================ one pass ================================ */
static void sum_ctrs(ctrs *t, const ctrs *a){
  t->Aprimes += a->Aprimes; t->Acands += a->Acands; t->Asieved += a->Asieved; t->Apre += a->Apre; t->Ashape += a->Ashape; t->Afinal += a->Afinal; t->Aunsolv += a->Aunsolv;
  t->Cprimes += a->Cprimes; t->Cpairs += a->Cpairs; t->Cnoinv += a->Cnoinv; t->Cincompat += a->Cincompat; t->Cinner += a->Cinner; t->Cp1 += a->Cp1; t->Cw1 += a->Cw1;
  t->Cfinal += a->Cfinal; t->Cunsolv += a->Cunsolv; t->hits += a->hits; t->longpairs += a->longpairs; t->Abig += a->Abig; t->tA += a->tA; t->tC += a->tC; }

static ctrs run_pass(double *wt){
  double t0 = wall(), c0 = pcpu();
  sieve_pass_init();
  int T = gThreads;
  tstate *TS = calloc(T, sizeof(tstate));
  for(int i = 0; i < T; i++){ bctx_init(&TS[i].b); mpz_inits(TS[i].n, TS[i].t1, NULL); TS[i].dead = malloc(gChunkA + 8); TS[i].e = malloc(8 * (PC.nSl + 8)); }
  ctrs plan_c; memset(&plan_c, 0, sizeof plan_c);
  const size_t BLK = 1 << 14;
  progA *PA = malloc(sizeof(progA) * BLK); u8 *BR = malloc(BLK);
  size_t capI = 1 << 16; item *IT = malloc(sizeof(item) * capI);
  double tlast = t0;
  for(size_t b0 = 0; b0 < NT; b0 += BLK){                                                         /* [P-ALL] */
    size_t b1 = b0 + BLK < NT ? b0 + BLK : NT;
    #pragma omp parallel for schedule(dynamic, 256) num_threads(T)
    for(size_t i = b0; i < b1; i++) BR[i - b0] = (u8)decide(i, &PA[i - b0]);
    size_t ni = 0;
    for(size_t i = b0; i < b1; i++){
      int br = BR[i - b0];
      if(br == 0){ plan_c.Aunsolv++; continue; }                                                /* [JP0] */
      if(br == 1){
        plan_c.Aprimes++; u64 cnt = PA[i - b0].count; if(cnt > 100000000ULL) plan_c.Abig++;
        for(u64 k = 0; k < cnt; k += gChunkA){                                                    /* [A-CHUNK] */
          if(ni == capI){ capI *= 2; IT = realloc(IT, sizeof(item) * capI); }
          IT[ni].type = 1; IT[ni].i3 = i; IT[ni].m0 = PA[i - b0].m0; IT[ni].step = PA[i - b0].step; IT[ni].k0 = k; IT[ni].k1 = (k + gChunkA < cnt) ? k + gChunkA : cnt; ni++;
        }
      } else {
        plan_c.Cprimes++;
        for(u64 j = 0; j < i; j += gChunkC){                                                      /* [C-ALLP2] */
          if(ni == capI){ capI *= 2; IT = realloc(IT, sizeof(item) * capI); }
          IT[ni].type = 2; IT[ni].i3 = i; IT[ni].k0 = j; IT[ni].k1 = (j + gChunkC < i) ? j + gChunkC : i; ni++;
        }
      }
    }
    #pragma omp parallel for schedule(dynamic, 1) num_threads(T)
    for(size_t q = 0; q < ni; q++){
      tstate *ts = &TS[omp_get_thread_num()]; double ta = wall();
      if(IT[q].type == 1){ runA(ts, &IT[q]); ts->c.tA += wall() - ta; } else { runC(ts, &IT[q]); ts->c.tC += wall() - ta; }
    }
    if(gProgress && wall() - tlast > 30){ tlast = wall(); fprintf(stderr, "progress p3 <= %llu (%.1f%% of table) elapsed %.0fs\n", (unsigned long long)Tp[b1 - 1], 100.0 * b1 / NT, tlast - t0); }
  }
  ctrs tot = plan_c;
  for(int i = 0; i < T; i++){ sum_ctrs(&tot, &TS[i].c); bctx_clear(&TS[i].b); mpz_clears(TS[i].n, TS[i].t1, NULL); free(TS[i].dead); free(TS[i].e); }
  free(TS); free(PA); free(BR); free(IT);
  tot.cpu = pcpu() - c0; *wt = wall() - t0; return tot; }

static void print_pass(const ctrs *c, double wt){
  char cap[32]; if(gCap == -2) strcpy(cap, "auto"); else snprintf(cap, sizeof cap, "%lld", gCap);
  printf("PASS D=%lld S=%s A=%s B=%s mode=%s cap=%s X=%llu inert=%zu Aprimes=%llu Acands=%llu Cprimes=%llu Cpairs=%llu Cinner=%llu prefilter=%llu hits=%llu wall=%.2fs\n",
    (long long)gD, PC.Ss, PC.As, PC.Bs, MODENAME[gMode], cap, (unsigned long long)gX, NT, (unsigned long long)c->Aprimes, (unsigned long long)c->Acands,
    (unsigned long long)c->Cprimes, (unsigned long long)c->Cpairs, (unsigned long long)c->Cinner, (unsigned long long)c->Apre, (unsigned long long)c->hits, wt);
  printf("DIAG Aunsolv=%llu Asieved=%llu Ashape=%llu Afinal=%llu Abig=%llu Cunsolv=%llu Cnoinv=%llu Cincompat=%llu Cp1=%llu Cw1=%llu Cfinal=%llu longpairs=%llu "
         "threadsec_A=%.2f threadsec_C=%.2f us_per_Acand=%.4f us_per_Cpair=%.4f sieve=%s cpu=%.2fs\n",
    (unsigned long long)c->Aunsolv, (unsigned long long)c->Asieved, (unsigned long long)c->Ashape, (unsigned long long)c->Afinal, (unsigned long long)c->Abig,
    (unsigned long long)c->Cunsolv, (unsigned long long)c->Cnoinv, (unsigned long long)c->Cincompat, (unsigned long long)c->Cp1, (unsigned long long)c->Cw1,
    (unsigned long long)c->Cfinal, (unsigned long long)c->longpairs, c->tA, c->tC, c->Acands ? 1e6 * c->tA / c->Acands : 0.0, c->Cpairs ? 1e6 * c->tC / c->Cpairs : 0.0,
    gSieve ? "on" : "off", c->cpu);
  fflush(stdout); }

/* plant pass bookkeeping */
static void plants_prepare(int pass){
  free(PK); PK = malloc(sizeof(pkey) * (NPL + 1)); NPK = 0;
  mpz_t t; mpz_init(t);
  for(size_t k = 0; k < NPL; k++){
    plant *P = &PL[k]; if(P->pass != pass) continue;
    P->found = 0; P->valid = 1; P->why[0] = 0;
    int64_t a = tab_find(P->p1), b = tab_find(P->p2), c = tab_find(P->p3);
    if(!(P->p1 < P->p2 && P->p2 < P->p3 && P->p3 <= gX)){ P->valid = 0; strcpy(P->why, "order/bound"); }
    else if(a < 0 || b < 0 || c < 0){ P->valid = 0; strcpy(P->why, "not inert table primes"); }
    else if(!w1_holds(get_M(a), P->p2, P->p3)){ P->valid = 0; strcpy(P->why, "W1 fails at p1"); }
    else if(!w1_holds(get_M(b), P->p1, P->p3)){ P->valid = 0; strcpy(P->why, "W1 fails at p2"); }
    else if(!w1_holds(get_M(c), P->p1, P->p2)){ P->valid = 0; strcpy(P->why, "W1 fails at p3"); }
    else if(!w2_holds(P->p1, P->p2, P->p3, t)){ P->valid = 0; strcpy(P->why, "W2 fails"); }
    PK[NPK].p3 = P->p3; PK[NPK].m = (u128)P->p1 * P->p2; PK[NPK].idx = k; NPK++;
  }
  mpz_clear(t);
  qsort(PK, NPK, sizeof(pkey), pkey_cmp); }
static void plants_report(int pass){
  size_t nv = 0, nf = 0, nmiss = 0, ninv = 0, ninvfound = 0, ndup = 0;
  for(size_t k = 0; k < NPL; k++){ plant *P = &PL[k]; if(P->pass != pass) continue;
    if(P->valid){ nv++; if(P->found) nf++; else nmiss++; if(P->found > 1) ndup++; } else { ninv++; if(P->found) ninvfound++; }
    printf("PLANT %zu %llu %llu %llu S=%s A=%s B=%s valid=%d found=%d%s%s\n", k, (unsigned long long)P->p1, (unsigned long long)P->p2, (unsigned long long)P->p3,
      P->S, P->A, P->B, P->valid, P->found, P->valid ? "" : " why=", P->why); }
  printf("PLANTS pass=%d valid=%zu found=%zu MISSED=%zu duplicates=%zu invalid=%zu invalid_found=%zu\n", pass, nv, nf, nmiss, ndup, ninv, ninvfound); fflush(stdout); }

/* ================================ self-test ================================ */
static u64 rs = 0x9E3779B97F4A7C15ULL;
static u64 rnd(void){ u64 z = (rs += 0x9E3779B97F4A7C15ULL); z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL; z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL; return z ^ (z >> 31); }
static u128 rnd128(int bits){ u128 x = ((u128)rnd() << 64) | rnd(); if(bits < 128) x &= (((u128)1) << bits) - 1; return x; }
static long st_fail = 0;
#define CHECK(cond, ...) do { if(!(cond)){ st_fail++; if(st_fail < 50){ printf("SELFTEST FAIL: "); printf(__VA_ARGS__); printf("\n"); } } } while(0)

/* direct reference BPSW (u64, literal transcription of src/bpsw.py; n < 2^31) */
static int ref_jacobi(int64_t a, int64_t n){ a %= n; if(a < 0) a += n; int r = 1;
  while(a){ while(a % 2 == 0){ a /= 2; if(n % 8 == 3 || n % 8 == 5) r = -r; } int64_t t = a; a = n; n = t; if(a % 4 == 3 && n % 4 == 3) r = -r; a %= n; }
  return n == 1 ? r : 0; }
static int64_t pmod(i128 a, int64_t n){ i128 r = a % n; if(r < 0) r += n; return (int64_t)r; }
static int ref_sf2(int64_t n){ int64_t d = n - 1, s = 0; while(d % 2 == 0){ d /= 2; s++; }
  int64_t x = 1, b = 2 % n, e = d; while(e){ if(e & 1) x = pmod((i128)x * b, n); b = pmod((i128)b * b, n); e >>= 1; }
  if(x == 1 || x == n - 1) return 1; for(int64_t i = 0; i < s - 1; i++){ x = pmod((i128)x * x, n); if(x == n - 1) return 1; } return 0; }
static int ref_sl(int64_t n){
  int64_t rt = (int64_t)sqrtl((long double)n); while(rt * rt > n) rt--; while((rt + 1) * (rt + 1) <= n) rt++;
  if(rt * rt == n) return 0;
  int64_t D = 5; while(ref_jacobi(D, n) != -1) D = D > 0 ? -D - 2 : -D + 2;
  int64_t P = 1, Q = (1 - D) / 4, e = n + 1, t = 0; while(e % 2 == 0){ e /= 2; t++; }
  int64_t U = 0, V = 2, Qk = 1, inv2 = (n + 1) / 2;
  int nb = 63 - __builtin_clzll((u64)e);
  for(int b = nb; b >= 0; b--){
    U = pmod((i128)U * V, n); V = pmod((i128)V * V - 2 * (i128)Qk, n); Qk = pmod((i128)Qk * Qk, n);
    if((e >> b) & 1){ int64_t U1 = pmod(((i128)P * U + V) * inv2, n), V1 = pmod(((i128)D * U + (i128)P * V) * inv2, n); U = U1; V = V1; Qk = pmod((i128)Qk * Q, n); }
  }
  if(U == 0 || V == 0) return 1;
  for(int64_t i = 0; i < t - 1; i++){ V = pmod((i128)V * V - 2 * (i128)Qk, n); Qk = pmod((i128)Qk * Qk, n); if(V == 0) return 1; }
  return 0; }
static int ref_bpsw(int64_t n){ return n > 1 && n % 2 == 1 && ref_sf2(n) && ref_sl(n); }
/* literal mpz reference (for big n) */
static int refz_jacobi(const mpz_t a0, const mpz_t n0){ mpz_t a, n, t; mpz_inits(a, n, t, NULL); mpz_mod(a, a0, n0); mpz_set(n, n0); int r = 1;
  while(mpz_sgn(a)){ while(mpz_even_p(a)){ mpz_tdiv_q_2exp(a, a, 1); unsigned long m8 = mpz_fdiv_ui(n, 8); if(m8 == 3 || m8 == 5) r = -r; }
    mpz_swap(a, n); if(mpz_fdiv_ui(a, 4) == 3 && mpz_fdiv_ui(n, 4) == 3) r = -r; mpz_mod(a, a, n); }
  int res = mpz_cmp_ui(n, 1) == 0 ? r : 0; mpz_clears(a, n, t, NULL); return res; }
static int refz_bpsw(const mpz_t n){
  if(mpz_cmp_ui(n, 1) <= 0 || mpz_even_p(n)) return 0;
  mpz_t d, x, nm1, two, U, V, Qk, inv2, t, Dz, e, U1, V1; mpz_inits(d, x, nm1, two, U, V, Qk, inv2, t, Dz, e, U1, V1, NULL); int res = 0;
  mpz_sub_ui(nm1, n, 1); mpz_set(d, nm1); long s = 0; while(mpz_even_p(d)){ mpz_tdiv_q_2exp(d, d, 1); s++; }
  mpz_set_ui(two, 2); mpz_powm(x, two, d, n); int sf = 0;
  if(mpz_cmp_ui(x, 1) == 0 || mpz_cmp(x, nm1) == 0) sf = 1;
  for(long i = 0; i < s - 1 && !sf; i++){ mpz_mul(x, x, x); mpz_mod(x, x, n); if(mpz_cmp(x, nm1) == 0) sf = 1; }
  if(!sf) goto out;
  mpz_sqrt(t, n); mpz_mul(t, t, t); if(mpz_cmp(t, n) == 0) goto out;
  long D = 5; for(;;){ mpz_set_si(Dz, D); if(refz_jacobi(Dz, n) == -1) break; D = D > 0 ? -D - 2 : -D + 2; }
  long Q = (1 - D) / 4; mpz_add_ui(e, n, 1); long tt = 0; while(mpz_even_p(e)){ mpz_tdiv_q_2exp(e, e, 1); tt++; }
  mpz_set_ui(U, 0); mpz_set_ui(V, 2); mpz_set_ui(Qk, 1); mpz_add_ui(inv2, n, 1); mpz_tdiv_q_2exp(inv2, inv2, 1);
  for(long b = (long)mpz_sizeinbase(e, 2) - 1; b >= 0; b--){
    mpz_mul(U1, U, V); mpz_mod(U1, U1, n); mpz_mul(V1, V, V); mpz_submul_ui(V1, Qk, 2); mpz_mod(V1, V1, n); mpz_set(U, U1); mpz_set(V, V1); mpz_mul(Qk, Qk, Qk); mpz_mod(Qk, Qk, n);
    if(mpz_tstbit(e, b)){ mpz_add(U1, U, V); mpz_mul(U1, U1, inv2); mpz_mod(U1, U1, n); mpz_mul_si(V1, U, D); mpz_add(V1, V1, V); mpz_mul(V1, V1, inv2); mpz_mod(V1, V1, n);
      mpz_set(U, U1); mpz_set(V, V1); mpz_mul_si(Qk, Qk, Q); mpz_mod(Qk, Qk, n); }
  }
  if(!mpz_sgn(U) || !mpz_sgn(V)){ res = 1; goto out; }
  for(long i = 0; i < tt - 1; i++){ mpz_mul(V, V, V); mpz_submul_ui(V, Qk, 2); mpz_mod(V, V, n); mpz_mul(Qk, Qk, Qk); mpz_mod(Qk, Qk, n); if(!mpz_sgn(V)){ res = 1; break; } }
out: mpz_clears(d, x, nm1, two, U, V, Qk, inv2, t, Dz, e, U1, V1, NULL); return res; }
/* naive U_k(1,Q) mod p (mpz, binary method via matrix-free doubling as in bpsw.py, independent of m64) */
static int refz_lucasU_is0(u64 k, long D, u64 p){
  mpz_t U, V, Qk, n, inv2, U1, V1; mpz_inits(U, V, Qk, n, inv2, U1, V1, NULL); mpz_set_ui(n, p); long Q = (1 - D) / 4;
  mpz_set_ui(U, 0); mpz_set_ui(V, 2); mpz_set_ui(Qk, 1); mpz_set_ui(inv2, (p + 1) / 2);
  for(int b = 63 - __builtin_clzll(k); b >= 0; b--){
    mpz_mul(U1, U, V); mpz_mod(U1, U1, n); mpz_mul(V1, V, V); mpz_submul_ui(V1, Qk, 2); mpz_mod(V1, V1, n); mpz_set(U, U1); mpz_set(V, V1); mpz_mul(Qk, Qk, Qk); mpz_mod(Qk, Qk, n);
    if((k >> b) & 1){ mpz_add(U1, U, V); mpz_mul(U1, U1, inv2); mpz_mod(U1, U1, n); mpz_mul_si(V1, U, D); mpz_add(V1, V1, V); mpz_mul(V1, V1, inv2); mpz_mod(V1, V1, n); mpz_set(U, U1); mpz_set(V, V1); mpz_mul_si(Qk, Qk, Q); mpz_mod(Qk, Qk, n); }
  }
  int r = mpz_sgn(U) == 0; mpz_clears(U, V, Qk, n, inv2, U1, V1, NULL); return r; }
static void u64_prime_factors(u64 x, u64 *f, int *nf){ *nf = 0; for(u64 q = 2; q * q <= x; q++) if(x % q == 0){ f[(*nf)++] = q; while(x % q == 0) x /= q; } if(x > 1) f[(*nf)++] = x; }
/* certificate check of ord and z for table entry index i (independent of the sieve) */
static int certify_entry(u64 p, u64 o, u64 z){
  mpz_t b, e, m, r; mpz_inits(b, e, m, r, NULL); int ok = 1; u64 f[64]; int nf;
  mpz_set_ui(b, 2); mpz_set_ui(m, p);
  if((p - 1) % o) ok = 0;
  mpz_set_ui(e, o); mpz_powm(r, b, e, m); if(mpz_cmp_ui(r, 1)) ok = 0;
  u64_prime_factors(o, f, &nf); for(int k = 0; k < nf; k++){ mpz_set_ui(e, o / f[k]); mpz_powm(r, b, e, m); if(!mpz_cmp_ui(r, 1)) ok = 0; }
  if((p + 1) % z) ok = 0;
  if(!refz_lucasU_is0(z, (long)gD, p)) ok = 0;
  u64_prime_factors(z, f, &nf); for(int k = 0; k < nf; k++) if(refz_lucasU_is0(z / f[k], (long)gD, p)) ok = 0;
  mpz_clears(b, e, m, r, NULL); return ok; }

static void selftest_arith(void){
  mpz_t a, b, n, r, s, e; mpz_inits(a, b, n, r, s, e, NULL);
  long N = 200000;
  /* mulmod128 in all three regimes */
  for(long it = 0; it < N; it++){
    int bits = (int)(it % 3 == 0 ? 1 + rnd() % 64 : it % 3 == 1 ? 65 + rnd() % 30 : 96 + rnd() % 31);
    u128 Nm = rnd128(bits) | 1; if(Nm < 2) Nm = 3; u128 x = rnd128(127) % Nm, y = rnd128(127) % Nm;
    u128 got = mulmod128(x, y, Nm); mpz_set_u128(a, x); mpz_set_u128(b, y); mpz_set_u128(n, Nm); mpz_mul(r, a, b); mpz_mod(r, r, n);
    u128 want; mpz_get_u128(&want, r); CHECK(got == want, "mulmod128 bits=%d", bits);
  }
  /* invmod128, including even and non-invertible cases */
  for(long it = 0; it < N; it++){
    int bits = 2 + (int)(rnd() % 125); u128 Nm = rnd128(bits); if(Nm < 1) Nm = 1; u128 x = rnd128(127);
    u128 got; int ok = invmod128(x, Nm, &got);
    mpz_set_u128(a, x); mpz_set_u128(n, Nm); mpz_mod(a, a, n); mpz_gcd(r, a, n);
    int want_ok = (mpz_cmp_ui(r, 1) == 0) || (Nm == 1);
    CHECK(ok == want_ok, "invmod128 existence N bits %d", bits);
    if(ok && want_ok && Nm > 1){ mpz_invert(r, a, n); u128 w; mpz_get_u128(&w, r); CHECK(got == w, "invmod128 value"); }
  }
  /* crt128 vs crt_mpz */
  for(long it = 0; it < N; it++){
    int b1 = 1 + (int)(rnd() % 70), b2 = 1 + (int)(rnd() % 70);
    u128 L1 = rnd128(b1); if(!L1) L1 = 1; u128 L2 = rnd128(b2); if(!L2) L2 = 1;
    if(it % 3 == 0){ u128 g = 1 + rnd() % 1000; if(L1 < ((u128)1 << 100)) L1 *= g; if(L2 < ((u128)1 << 100)) L2 *= g; }
    u128 r1 = rnd128(127) % L1, r2 = rnd128(127) % L2;
    if(it % 2 == 0){ u128 x = rnd128(120); r1 = x % L1; r2 = x % L2; }
    u128 rr, LL; int st = crt128(r1, L1, r2, L2, &rr, &LL);
    mpz_t z1, z2, zr, zL; mpz_inits(z1, z2, zr, zL, NULL); mpz_set_u128(a, r1); mpz_set_u128(b, L1); mpz_set_u128(z1, r2); mpz_set_u128(z2, L2);
    int okz = crt_mpz(zr, zL, a, b, z1, z2);
    CHECK((st != 0) == okz, "crt existence");
    if(st == 1 && okz){ u128 w1, w2; mpz_get_u128(&w1, zr); mpz_get_u128(&w2, zL); CHECK(rr == w1 && LL == w2, "crt value"); CHECK(rr % L1 == r1 && rr % L2 == r2, "crt congruences"); }
    if(st == 2 && okz) CHECK(mpz_sizeinbase(zL, 2) >= 127, "crt overflow flag");
    mpz_clears(z1, z2, zr, zL, NULL);
  }
  /* Fermat prefilter, both Montgomery widths, multi-limb exponents; also exponents that make the result 1 */
  for(long it = 0; it < N; it++){
    int bits = (it & 1) ? 5 + (int)(rnd() % 58) : 64 + (int)(rnd() % 62);
    u128 m = rnd128(bits) | 1; if(m < 15) m = 15; if(bits >= 64) m |= ((u128)1) << 63;
    int nl = 1 + (int)(rnd() % 4); u64 ev[6]; for(int i = 0; i < nl; i++) ev[i] = rnd(); if(it % 7 == 0) ev[nl - 1] = 0;
    mpz_set_u128(n, m); mpz_import(e, nl, -1, 8, 0, 0, ev);
    if(it % 5 == 0){ /* force e = multiple of ord: use e = lambda-like multiple: e = (m-1)*k for prime m, else random */
      if(mpz_probab_prime_p(n, 20)){ mpz_sub_ui(e, n, 1); mpz_mul_ui(e, e, 1 + rnd() % 1000); size_t c2 = 0; memset(ev, 0, sizeof ev); mpz_export(ev, &c2, -1, 8, 0, 0, e); nl = (int)(c2 ? c2 : 1); } }
    mpz_set_ui(b, 2); mpz_powm(r, b, e, n); int want = mpz_cmp_ui(r, 1) == 0; int got = fermat_m(m, ev, nl);
    CHECK(got == want, "fermat_m bits=%d", bits);
  }
  /* m128 CIOS product vs reference REDC and vs mpz (a*b*2^-128 mod n) */
  for(long it = 0; it < N; it++){
    u128 nn = rnd128(64 + (int)(rnd() % 62)) | 1 | (((u128)1) << 63); m128 c; m128_init(&c, nn);
    u128 x = rnd128(127) % nn, y = rnd128(127) % nn; if(it % 5 == 0) x = nn - 1; if(it % 7 == 0) y = nn - 1 - (rnd() % 3);
    u128 g1 = m128_mul(x, y, &c), g2 = m128_mul_ref(x, y, &c);
    mpz_set_u128(a, x); mpz_set_u128(b, y); mpz_set_u128(n, nn); mpz_mul(r, a, b); mpz_set_ui(s, 1); mpz_mul_2exp(s, s, 128); mpz_invert(s, s, n); mpz_mul(r, r, s); mpz_mod(r, r, n);
    u128 w; mpz_get_u128(&w, r); CHECK(g1 == w && g2 == w, "m128 product");
  }
  /* 4-lane prefilters vs mpz (random lanes, plus lanes that must give 1) */
  for(long it = 0; it < N / 4; it++){
    int wide = (int)(it & 1); u128 mv[4]; u64 nv[4], e0[4], e1[4]; int got[4];
    for(int i = 0; i < 4; i++){
      int bits = wide ? 64 + (int)(rnd() % 62) : 5 + (int)(rnd() % 58);
      u128 m = rnd128(bits) | 1; if(m < 15) m = 15; if(wide) m |= ((u128)1) << 63; mv[i] = m; nv[i] = (u64)m;
      u128 E = rnd128(1 + (int)(rnd() % 120)); if(rnd() % 8 == 0) E = 0;
      mpz_set_u128(n, m); if(rnd() % 3 == 0 && mpz_probab_prime_p(n, 20)){ E = (m - 1) * (1 + rnd() % 64); }
      e0[i] = (u64)E; e1[i] = (u64)(E >> 64);
    }
    if(wide) m128_pow2_is1_x4(mv, e0, e1, got); else m64_pow2_is1_x4(nv, e0, e1, got);
    for(int i = 0; i < 4; i++){ mpz_set_u128(n, mv[i]); mpz_set_u128(e, ((u128)e1[i] << 64) | e0[i]); mpz_set_ui(b, 2); mpz_powm(r, b, e, n);
      CHECK(got[i] == (mpz_cmp_ui(r, 1) == 0), "4-lane fermat wide=%d", wide); }
  }
  /* fastmod32 / modsmall vs % */
  for(long it = 0; it < N; it++){
    const svprime *P = &SV[rnd() % NSV]; u128 x = rnd128(90); u32 a = (u32)rnd();
    CHECK(fastmod32(a, P->Ml, P->l) == a % P->l && fastmod32(a, P->Mo, P->o) == a % P->o, "fastmod32");
    CHECK(modsmall(x, P->Ml, P->l, P->t32l, P->t64l) == (u32)(x % P->l) && modsmall(x, P->Mo, P->o, P->t32o, P->t64o) == (u32)(x % P->o), "modsmall");
    u32 r0 = (u32)(rnd() % P->l); if(r0) CHECK((u64)r0 * P->invl[r0] % P->l == 1, "invl");
    u32 r1 = (u32)(rnd() % P->o); if(P->invo[r1] != 0xFFFF) CHECK(P->o == 1 || (u64)r1 * P->invo[r1] % P->o == 1, "invo"); else CHECK(gcd64(r1, P->o) != 1, "invo sentinel");
  }
  /* specifically: primes m (both widths), e = m - 1 -> must be 1 */
  for(long it = 0, cnt = 0; cnt < 2000 && it < 200000; it++){
    int bits = (it & 1) ? 20 + (int)(rnd() % 42) : 64 + (int)(rnd() % 60);
    u128 m = rnd128(bits) | 1; if(bits >= 64) m |= ((u128)1) << 63; mpz_set_u128(n, m); if(!mpz_probab_prime_p(n, 20)) continue; cnt++;
    u64 ev[3] = {0, 0, 0}; size_t c2; mpz_sub_ui(e, n, 1); mpz_export(ev, &c2, -1, 8, 0, 0, e);
    CHECK(fermat_m(m, ev, (int)c2) == 1, "fermat_m on prime");
  }
  /* build_nm1 with random S */
  {
    const char *Ss[] = {"1", "3", "1000000007", "340282366920938463463374607431768211457", "99999999999999999999999999999999999999999999999999999999999"};
    for(int si = 0; si < 5; si++){
      if(!pass_setup(Ss[si], "1", "1")) die("pass_setup");
      for(long it = 0; it < 20000; it++){
        u64 p3 = 3 + rnd() % ((u64)1 << 34); u128 m = rnd128(68);
        u64 ev[80]; int nl = build_nm1(p3, m, ev);
        mpz_import(e, nl, -1, 8, 0, 0, ev); make_n_m(n, p3, m); mpz_sub_ui(n, n, 1); CHECK(mpz_cmp(e, n) == 0, "build_nm1 S=%s", Ss[si]);
        u128 q = rnd128(1 + (int)(rnd() % 94)); if(!q) q = 1; u128 sm = S_mod(q); mpz_set_u128(a, q); mpz_mod(b, PC.S, a); u128 w; mpz_get_u128(&w, b); CHECK(sm == w, "S_mod");
      }
    }
    pass_setup("1", "1", "1");
  }
  /* Lucas U (Montgomery) vs mpz reference */
  for(long it = 0, cnt = 0; cnt < 20000 && it < 1000000; it++){
    u64 p = (rnd() % ((u64)1 << 40)) | 1; if(p < 5) continue; mpz_set_ui(n, p); if(!mpz_probab_prime_p(n, 20)) continue;
    long Ds[] = {5, -7, 13, -11, 17, -15, 21}; long D = Ds[rnd() % 7]; if(((u64)labs(D)) % p == 0) continue; long Q = (1 - D) / 4; if(Q == 0 || (u64)labs(Q) % p == 0) continue;
    cnt++; m64 c; m64_init(&c, p); u64 Dm = m64_to(s64mod(D, p), &c), Qm = m64_to(s64mod(Q, p), &c);
    u64 k = 1 + rnd() % (2 * p); if(it % 3 == 0) k = p - jacobi_s64(D, p);
    CHECK(m64_lucasU_is0(k, Dm, Qm, &c) == refz_lucasU_is0(k, D, p), "lucasU p=%llu D=%ld k=%llu", (unsigned long long)p, D, (unsigned long long)k);
  }
  /* is_prime_u64, rho, split helpers */
  for(long it = 0; it < 100000; it++){
    u64 x = (rnd() % ((u64)1 << 50)) | 1; mpz_set_ui(n, x); CHECK(is_prime_u64(x) == (mpz_probab_prime_p(n, 30) != 0), "is_prime_u64 %llu", (unsigned long long)x);
    if(!is_prime_u64(x) && x > 1 && it % 10 == 0){ u64 f = 0; for(int i = 0; i < 24; i++) if(x % SMALLP[i] == 0){ f = SMALLP[i]; break; } if(!f) f = rho_u64(x); CHECK(f > 1 && f < x && x % f == 0, "rho"); }
  }
  mpz_clears(a, b, n, r, s, e, NULL);
}

static void selftest_bpsw(void){
  /* all odd n < 1e6: production GMP vs u64 literal reference vs sieve primality */
  const int64_t LIM = 1000000; u8 *comp = calloc(LIM + 1, 1); for(int64_t i = 2; i * i <= LIM; i++) if(!comp[i]) for(int64_t j = i * i; j <= LIM; j += i) comp[j] = 1;
  bctx b; bctx_init(&b); mpz_t n; mpz_init(n);
  long nsp = 0, nsl = 0, disagree = 0, nonprime_pass = 0; int64_t spl[200], sll[400];
  for(int64_t x = 1; x < LIM; x += 2){
    mpz_set_si(n, x); int g = bpsw_mpz(&b, n), r = ref_bpsw(x), pr = (x > 1 && !comp[x]);
    if(g != r) disagree++; if(g != pr) nonprime_pass++;
    if(x > 1){ int sf = sprp2_mpz(&b, n); long D; int sl = 0; if(selfridge_mpz(n, &D)){ mpz_add_ui(b.dl, n, 1); sl = slprp_core(&b, n, D); }
      int rsf = ref_sf2(x), rsl = ref_sl(x); if(sf != rsf || sl != rsl) disagree++;
      if(comp[x] && sf){ if(nsp < 200) spl[nsp] = x; nsp++; }
      if(comp[x] && sl){ if(nsl < 400) sll[nsl] = x; nsl++; } }
  }
  CHECK(disagree == 0, "BPSW production vs reference disagreements: %ld", disagree);
  CHECK(nonprime_pass == 0, "BPSW vs primality below 1e6: %ld", nonprime_pass);
  printf("SELFTEST bpsw: odd n < 1e6: disagreements(prod vs ref, incl. each half) = %ld, BPSW != isprime: %ld\n", disagree, nonprime_pass);
  printf("SELFTEST bpsw: strong 2-psp below 1e6: %ld (first: ", nsp); for(int i = 0; i < nsp && i < 12; i++) printf("%lld ", (long long)spl[i]); printf(")\n");
  printf("SELFTEST bpsw: strong Lucas (Selfridge) psp below 1e6: %ld (first: ", nsl); for(int i = 0; i < nsl && i < 12; i++) printf("%lld ", (long long)sll[i]); printf(")\n");
  CHECK(nsp == 46, "expected 46 strong 2-psp below 1e6 (OEIS A001262), got %ld", nsp);
  const int64_t knownsl[] = {5459, 5777, 10877, 16109, 18971, 22499, 24569, 25199, 40309, 58519, 75077, 97439};
  for(int i = 0; i < 12; i++){ int f = 0; for(int k = 0; k < nsl && k < 400; k++) if(sll[k] == knownsl[i]) f = 1; CHECK(f, "known strong Lucas psp %lld not detected", (long long)knownsl[i]);
    mpz_set_si(n, knownsl[i]); CHECK(bpsw_mpz(&b, n) == 0, "known slpsp passes BPSW"); }
  for(int i = 0; i < nsp && i < 200; i++){ mpz_set_si(n, spl[i]); CHECK(bpsw_mpz(&b, n) == 0, "2-spsp %lld passes BPSW", (long long)spl[i]); }
  /* 1e5 random 100-bit odd n: production vs literal mpz reference vs GMP probable-prime */
  long d2 = 0, d3 = 0, nprime = 0;
  for(long it = 0; it < 100000; it++){
    u128 x = rnd128(100) | 1 | ((u128)1 << 99); mpz_set_u128(n, x);
    if(it % 4 == 0){ mpz_nextprime(n, n); } /* include many primes */
    if(it % 50 == 1){ mpz_set_u128(n, rnd128(50) | 1); mpz_mul(n, n, n); } /* perfect squares */
    int g = bpsw_mpz(&b, n), r = refz_bpsw(n), pp = mpz_probab_prime_p(n, 30) != 0;
    if(g != r) d2++; if(g != pp) d3++; nprime += pp;
  }
  CHECK(d2 == 0, "random 100-bit: production vs mpz reference: %ld", d2);
  CHECK(d3 == 0, "random 100-bit: BPSW vs probable prime: %ld", d3);
  printf("SELFTEST bpsw: 1e5 random ~100-bit odd n (%ld primes): prod-vs-ref disagreements %ld, BPSW-vs-probab_prime %ld\n", nprime, d2, d3);
  /* hit reporter / factorisation printer smoke test */
  { char buf[1024]; mpz_set_str(n, "3825123056546413051", 10); factor_string(n, buf, sizeof buf); printf("SELFTEST factor_string(3825123056546413051) = %s\n", buf);
    CHECK(strcmp(buf, "149491 * 747451 * 34233211") == 0, "factor_string"); }
  bctx_clear(&b); mpz_clear(n); free(comp);
}

static void selftest_table(void){
  long Ds[] = {5, -7, 13, -11, 17};
  for(int di = 0; di < 5; di++){
    gD = Ds[di]; gQ = (1 - gD) / 4;
    u64 X = 60000;
    gSeg = 4096; build_table(X);
    size_t n1 = NT; u64 *cp = malloc(8 * NT); u32 *c1 = malloc(4 * NT), *c2 = malloc(4 * NT); memcpy(cp, Tp, 8 * NT); memcpy(c1, Ti1, 4 * NT); memcpy(c2, Ti2, 4 * NT);
    gSeg = (u64)1 << 18; build_table(X);
    CHECK(n1 == NT && !memcmp(cp, Tp, 8 * NT) && !memcmp(c1, Ti1, 4 * NT) && !memcmp(c2, Ti2, 4 * NT), "table independent of segment size (D=%ld)", gD);
    free(cp); free(c1); free(c2);
    /* completeness of the inert list and brute-force ord / z */
    size_t k = 0; long bad = 0;
    for(u64 p = 3; p <= X; p += 2){
      int pr = 1; for(u64 q = 3; q * q <= p; q += 2) if(p % q == 0){ pr = 0; break; }
      if(!pr) continue;
      if(jacobi_s64(gD, p) != -1 || (u64)llabs(gQ) % p == 0) continue;
      if(k >= NT || Tp[k] != p){ bad++; continue; }
      u64 o = 1, x = 2 % p; while(x != 1){ x = x * 2 % p; o++; }
      u64 z = 1; { int64_t Um = 0, U = 1; while(U != 0){ int64_t Un = (U - (int64_t)(((i128)gQ * Um) % (int64_t)p)) % (int64_t)p; if(Un < 0) Un += p; Um = U; U = Un; z++; } }
      if(o != tab_ord(k) || z != tab_z(k)) bad++;
      k++;
    }
    CHECK(bad == 0 && k == NT, "table D=%ld X=%llu: bad=%ld k=%zu NT=%zu", gD, (unsigned long long)X, bad, k, NT);
    printf("SELFTEST table D=%ld X=%llu: %zu inert primes, brute-force ord/z mismatches %ld\n", gD, (unsigned long long)X, NT, bad);
  }
  /* certificate checks of single segments at large p (sieve factorisation at scale) */
  gD = 5; gQ = -1;
  u64 los[] = { 4294967296ULL - 300000, 9999700001ULL, 17179869184ULL - 200000 };
  for(int s = 0; s < 3; s++){
    u64 lo = los[s] | 1, hi = lo + 60000; init_base_primes(hi + 2);
    segws w; u64 sv = gSeg; gSeg = 1 << 18; segws_alloc(&w); seg_process(lo, hi, &w); gSeg = sv;
    long bad = 0, cnt = 0, missing = 0; size_t k = 0;
    mpz_t z; mpz_init(z);
    for(u64 p = lo; p < hi; p++){
      mpz_set_ui(z, p); if(!(p & 1) || !mpz_probab_prime_p(z, 30) || jacobi_s64(gD, p) != -1) continue;
      if(k >= w.nout || w.op[k] != p){ missing++; continue; }
      u64 o = (p - 1) / w.o1[k], zz = (p + 1) / w.o2[k]; if(!certify_entry(p, o, zz)) bad++; cnt++; k++;
    }
    CHECK(bad == 0 && missing == 0 && k == w.nout, "segment certificate lo=%llu bad=%ld missing=%ld", (unsigned long long)lo, bad, missing);
    printf("SELFTEST segment [%llu, %llu): %ld inert primes, ord/z certificates failed %ld, missing %ld\n", (unsigned long long)lo, (unsigned long long)hi, cnt, bad, missing);
    mpz_clear(z); segws_free(&w);
  }
}

static void selftest_sieve(void){
  /* every candidate killed by the sieve must fail Fermat mod m; random progressions from a real table */
  gD = 5; gQ = -1; gMode = MODE_BPSW; build_table(200000);
  mpz_t n, e, b, r, m; mpz_inits(n, e, b, r, m, NULL); mpz_set_ui(b, 2);
  u8 *dead = malloc(1 << 16); long killed = 0, wrong = 0, tot = 0;
  const char *Ss[] = {"1", "1009"};
  for(int si = 0; si < 2; si++){
    pass_setup(Ss[si], "1", "1"); sieve_pass_init();
    for(int it = 0; it < 3000; it++){
      size_t i3 = 100 + rnd() % (NT - 100); progA pa; if(!setupA(i3, &pa) || !pa.count) continue;
      u64 len = pa.count < 4096 ? pa.count : 4096; u64 p3 = Tp[i3];
      sieveA(p3, pa.m0, pa.step, len, dead);
      u128 mm = pa.m0;
      for(u64 d = 0; d < len; d++, mm += pa.step){
        tot++; mpz_set_u128(m, mm); make_n_m(n, p3, mm); mpz_sub_ui(e, n, 1); mpz_powm(r, b, e, m); int pass = mpz_cmp_ui(r, 1) == 0;
        if(dead[d]){ killed++; if(pass) wrong++; }
        u64 ev[16]; int nl = build_nm1(p3, mm, ev); if(fermat_m(mm, ev, nl) != pass) wrong++;
      }
    }
  }
  pass_setup("1", "1", "1");
  CHECK(wrong == 0, "sieve killed a Fermat-mod-m survivor or fermat_m mismatch: %ld", wrong);
  printf("SELFTEST sieve: %ld candidates, %ld killed by small-prime sieve, wrong kills %ld\n", tot, killed, wrong);
  /* progression count vs direct enumeration (small p3, with split parts) */
  long bad = 0, ntest = 0;
  const char *SAB[][3] = { {"1","1","1"}, {"1009","28","6"}, {"11","5","1"}, {"1","4","6"} };
  for(int si = 0; si < 4; si++){
    if(!pass_setup(SAB[si][0], SAB[si][1], SAB[si][2])) continue;
    for(size_t i3 = 0; i3 < 60 && i3 < NT; i3++){
      progA pa; int ok = setupA(i3, &pa); u64 p3 = Tp[i3]; u128 M3 = get_M(i3); u64 cnt = 0; u128 first = 0;
      for(u128 x = 15; x < (u128)p3 * p3; x++){
        if(!(x & 1)) continue;
        if(!w1_holds(M3, (u64)x, 1)) continue;
        mpz_set_u128(n, x); mpz_mul_ui(n, n, p3); mpz_mul(n, n, PC.S);
        if(!PC.W1){ mpz_mod(r, n, PC.W); if(mpz_cmp(r, PC.nu)) continue; }
        if(!cnt) first = x; cnt++;
      }
      ntest++;
      if(!ok){ if(cnt) bad++; continue; }
      if(cnt != pa.count || (cnt && first != pa.m0)) bad++;
    }
  }
  pass_setup("1", "1", "1");
  CHECK(bad == 0, "branch-A progression vs direct enumeration: %ld bad", bad);
  printf("SELFTEST progression: %ld (p3, S) cases checked against direct enumeration of m < p3^2, bad %ld\n", ntest, bad);
  free(dead); mpz_clears(n, e, b, r, m, NULL);
}

static int hit_cmp(const void *a, const void *b){ const hitrec *x = a, *y = b;
  if(x->p3 != y->p3) return x->p3 < y->p3 ? -1 : 1; if(x->p1 != y->p1) return x->p1 < y->p1 ? -1 : 1; if(x->p2 != y->p2) return x->p2 < y->p2 ? -1 : 1; if(x->m != y->m) return x->m < y->m ? -1 : 1; return 0; }
static void selftest_engine(void){
  /* cap independence of the control modes (and bpsw) on a small range */
  gD = 5; gQ = -1; gX = 20000; build_table(gX);
  int modes[] = { MODE_WEAK, MODE_FERMAT, MODE_LUCAS, MODE_BPSW, MODE_FERMAT, MODE_LUCAS };
  const char *msab[][3] = { {"1","1","1"}, {"1","1","1"}, {"1","1","1"}, {"1","1","1"}, {"11","10","1"}, {"11","1","10"} };
  long long caps[] = { -1, 0, -2, 40, 400 };
  int sq = gQuiet; gQuiet = 1;
  for(int mi = 0; mi < 6; mi++){
    gMode = modes[mi]; hitrec *ref = NULL; size_t nref = 0;
    for(int ci = 0; ci < 5; ci++){
      gCap = caps[ci]; NHITS = 0; pass_setup(msab[mi][0], msab[mi][1], msab[mi][2]); double wt; ctrs c = run_pass(&wt);
      qsort(HITS, NHITS, sizeof(hitrec), hit_cmp);
      if(ci == 0){ ref = malloc(sizeof(hitrec) * (NHITS + 1)); memcpy(ref, HITS, sizeof(hitrec) * NHITS); nref = NHITS; }
      else { int same = (NHITS == nref); for(size_t k = 0; same && k < NHITS; k++) if(hit_cmp(&ref[k], &HITS[k])) same = 0;
        CHECK(same, "mode %s cap %lld: %zu hits vs %zu at cap -1 (or different triples)", MODENAME[gMode], gCap, NHITS, nref); }
      printf("SELFTEST engine mode=%s S=%s A=%s B=%s X=%llu cap=%lld: hits=%zu Aprimes=%llu Cprimes=%llu Acands=%llu Cpairs=%llu\n", MODENAME[gMode], PC.Ss, PC.As, PC.Bs, (unsigned long long)gX, gCap, NHITS,
        (unsigned long long)c.Aprimes, (unsigned long long)c.Cprimes, (unsigned long long)c.Acands, (unsigned long long)c.Cpairs);
    }
    free(ref);
  }
  /* synthetic plants: overrides M(p_i) := gcd(p_i^2 - 1, n/p_i - 1); split parts chosen so that (W2) holds:
     pass 0: S = 1; pass 1: S odd with A = gcd(n-1, K), B = gcd(n+1, K) (W < 2^63, u128 CRT paths);
     pass 2: S odd with A = n - 1, B = (n+1)/2 (W > 2^127: GMP CRT in branch A, deferred (W2) in branch C) */
  gMode = MODE_PLANT; NOV = 0; NPL = 0; free(OV); free(PL); OV = malloc(sizeof(ovr) * 64); PL = malloc(sizeof(plant) * 16);
  { size_t idx[][3] = { {3, 10, 40}, {20, 70, 150}, {5, 200, 290}, {7, 30, 400}, {50, 120, 600}, {9, 60, 500} };
    const char *Sv[] = { "1", "1", "1000000007", "100000000003", "1", "1000000000000000000000007" };
    int kind[] = { 0, 0, 1, 2, 1, 2 };
    mpz_t nz, t, K, Az, Bz; mpz_inits(nz, t, K, Az, Bz, NULL);
    mpz_set_str(K, "363896953291872000", 10);   /* 2^8 3^5 5^3 7^2 11 13 17 19 23 29 31 */
    for(int k = 0; k < 6; k++){
      u64 p[3] = { Tp[idx[k][0]], Tp[idx[k][1]], Tp[idx[k][2]] };
      mpz_set_str(t, Sv[k], 10); mpz_set_ui(nz, p[0]); mpz_mul_ui(nz, nz, p[1]); mpz_mul_ui(nz, nz, p[2]); mpz_mul(nz, nz, t);
      for(int i = 0; i < 3; i++){ mpz_t c, g; mpz_inits(c, g, NULL); mpz_divexact_ui(c, nz, p[i]); mpz_sub_ui(c, c, 1); mpz_set_ui(g, p[i]); mpz_mul_ui(g, g, p[i]); mpz_sub_ui(g, g, 1); mpz_gcd(g, g, c);
        u128 M; mpz_get_u128(&M, g); OV[NOV].p = p[i]; OV[NOV].M = M; NOV++; mpz_clears(c, g, NULL); }
      if(kind[k] == 0){ mpz_set_ui(Az, 1); mpz_set_ui(Bz, 1); }
      else if(kind[k] == 1){ mpz_sub_ui(t, nz, 1); mpz_gcd(Az, t, K); mpz_add_ui(t, nz, 1); mpz_gcd(Bz, t, K); }
      else { mpz_sub_ui(Az, nz, 1); mpz_add_ui(Bz, nz, 1); mpz_tdiv_q_2exp(Bz, Bz, 1); }
      memset(&PL[NPL], 0, sizeof(plant)); PL[NPL].p1 = p[0]; PL[NPL].p2 = p[1]; PL[NPL].p3 = p[2];
      strcpy(PL[NPL].S, Sv[k]); mpz_get_str(PL[NPL].A, 10, Az); mpz_get_str(PL[NPL].B, 10, Bz); NPL++;
    }
    mpz_clears(nz, t, K, Az, Bz, NULL);
    qsort(OV, NOV, sizeof(ovr), ovr_cmp);
    for(size_t i = 1; i < NOV; i++) if(OV[i].p == OV[i - 1].p) die("selftest plant override clash");
    /* passes: plants with identical (S, A, B) share a pass */
    int np = 0; for(size_t k = 0; k < NPL; k++){ int f = -1; for(size_t j = 0; j < k; j++) if(!strcmp(PL[j].S, PL[k].S) && !strcmp(PL[j].A, PL[k].A) && !strcmp(PL[j].B, PL[k].B)) f = PL[j].pass; PL[k].pass = f >= 0 ? f : np++; }
    for(int ci = 0; ci < 5; ci++){
      gCap = caps[ci]; int allfound = 1;
      for(int q = 0; q < np; q++){
        size_t k0 = 0; while(PL[k0].pass != q) k0++;
        if(!pass_setup(PL[k0].S, PL[k0].A, PL[k0].B)) die("selftest plant pass setup");
        plants_prepare(q); NHITS = 0; double wt; run_pass(&wt);
        for(size_t k = 0; k < NPL; k++) if(PL[k].pass == q){ CHECK(PL[k].valid, "synthetic plant %zu invalid: %s", k, PL[k].why); if(PL[k].found != 1) allfound = 0;
          if(ci == 0) printf("SELFTEST plant %zu: p = %llu %llu %llu S=%s Wbits=%zu Wsmall=%d\n", k, (unsigned long long)PL[k].p1, (unsigned long long)PL[k].p2, (unsigned long long)PL[k].p3, PL[k].S, mpz_sizeinbase(PC.W, 2), PC.Wsmall); }
      }
      CHECK(allfound, "synthetic plants not all found exactly once at cap %lld", gCap);
      printf("SELFTEST plant cap=%lld: %zu plants in %d passes found exactly once: %s\n", gCap, NPL, np, allfound ? "yes" : "NO");
    }
  }
  NOV = 0; NPL = 0; gQuiet = sq; gMode = MODE_BPSW; gCap = -2;
}

static int selftest(void){
  double t0 = wall();
  init_sieve_primes(gSieveB);
  selftest_arith(); printf("SELFTEST arithmetic done (%ld failures so far) %.1fs\n", st_fail, wall() - t0); fflush(stdout);
  selftest_bpsw(); printf("SELFTEST bpsw done (%ld failures so far) %.1fs\n", st_fail, wall() - t0); fflush(stdout);
  selftest_table(); printf("SELFTEST table done (%ld failures so far) %.1fs\n", st_fail, wall() - t0); fflush(stdout);
  selftest_sieve(); printf("SELFTEST sieve/progression done (%ld failures so far) %.1fs\n", st_fail, wall() - t0); fflush(stdout);
  selftest_engine(); printf("SELFTEST engine done (%ld failures so far) %.1fs\n", st_fail, wall() - t0);
  printf("SELFTEST %s: %ld failures, %.1fs\n", st_fail ? "FAILED" : "PASSED", st_fail, wall() - t0);
  return st_fail ? 1 : 0; }



/* ================================ plan-only cost estimate ================================ */
static void plan_only(void){
  double t0 = wall();
  double sumN = 0, sumW = 0, pairs = 0; u64 nA = 0, nC = 0, nU = 0; u64 hist[24]; memset(hist, 0, sizeof hist);
  const int TOP = 16; u64 topc[16]; u64 topp[16]; for(int k = 0; k < TOP; k++){ topc[k] = 0; topp[k] = 0; }
  double heavyA = 0; u64 nheavy = 0;
  #pragma omp parallel num_threads(gThreads)
  {
    double sN = 0, sW = 0, pr = 0, hA = 0; u64 a = 0, c = 0, u = 0, nh = 0, h[24]; memset(h, 0, sizeof h); u64 tc[16], tp[16]; for(int k = 0; k < TOP; k++){ tc[k] = 0; tp[k] = 0; }
    #pragma omp for schedule(dynamic, 4096)
    for(size_t i = 0; i < NT; i++){
      progA pa; int br = decide(i, &pa);
      u64 cnt = pa.count; int dec = 0; while(dec < 23 && cnt >= ((u64)1 << (dec + 1))) dec++; if(br) h[cnt ? dec + 1 : 0]++;
      if(cnt > 1000000){ nh++; hA += cnt; }
      int k = TOP - 1; if(cnt > tc[k]){ while(k > 0 && cnt > tc[k - 1]){ tc[k] = tc[k - 1]; tp[k] = tp[k - 1]; k--; } tc[k] = cnt; tp[k] = Tp[i]; }
      if(br == 0){ u++; continue; }
      if(br == 1){ a++; double p3 = (double)Tp[i], fw = 1.0 - 9.223372036854775808e18 / (p3 * p3); if(fw < 0) fw = 0; sW += cnt * fw; sN += cnt * (1 - fw); }
      else { c++; pr += (double)i; }
    }
    #pragma omp critical
    { sumN += sN; sumW += sW; pairs += pr; nA += a; nC += c; nU += u; heavyA += hA; nheavy += nh; for(int k = 0; k < 24; k++) hist[k] += h[k];
      for(int k = 0; k < TOP; k++){ u64 cc = tc[k], pp = tp[k]; int q = TOP - 1; if(cc > topc[q]){ while(q > 0 && cc > topc[q - 1]){ topc[q] = topc[q - 1]; topp[q] = topp[q - 1]; q--; } topc[q] = cc; topp[q] = pp; } } }
  }
  double cpuA = sumN * gCtest + sumW * gCtestW, cpuC = pairs * gCpair;
  printf("PLAN D=%lld X=%llu inert=%zu Aprimes=%llu Cprimes=%llu unsolv=%llu Acands_narrow=%.4g Acands_wide=%.4g Cpairs=%.4g time=%.1fs\n", (long long)gD, (unsigned long long)gX, NT,
    (unsigned long long)nA, (unsigned long long)nC, (unsigned long long)nU, sumN, sumW, pairs, wall() - t0);
  printf("PLAN estimate (CPU, us/cand narrow %.3f wide %.3f, us/pair %.3f): branch A %.0f s, branch C %.0f s, total %.0f CPU-s = %.2f h on 4 uncontended cores\n",
    gCtest, gCtestW, gCpair, cpuA / 1e6, cpuC / 1e6, (cpuA + cpuC) / 1e6, (cpuA + cpuC) / 1e6 / 4 / 3600);
  printf("PLAN heavy tail: %llu p3 with count_A > 1e6, total count %.4g\n", (unsigned long long)nheavy, heavyA);
  printf("PLAN count_A histogram (bucket k: 2^(k-1) <= count < 2^k; bucket 0: count 0):"); for(int k = 0; k < 24; k++) if(hist[k]) printf(" %d:%llu", k, (unsigned long long)hist[k]); printf("\n");
  printf("PLAN largest count_A:"); for(int k = 0; k < TOP; k++) if(topc[k]) printf(" %llu@%llu", (unsigned long long)topc[k], (unsigned long long)topp[k]); printf("\n");
  fflush(stdout); }

/* ================================ micro-benchmarks (single thread, CPU time) ================================ */
static void bench(void){
  init_sieve_primes(gSieveB);
  gD = 5; gQ = -1; gMode = MODE_BPSW; pass_setup("1", "1", "1");
  const int N = 1000000; u128 *ms = malloc(sizeof(u128) * N); u64 *ps = malloc(8 * N); volatile u64 sink = 0;
  for(int i = 0; i < N; i++){ ps[i] = (rnd() % 100000000ULL) | 1; ms[i] = ((u128)(rnd() % ((u64)ps[i] * ps[i])) | 1); if(ms[i] < 15) ms[i] = 15; }
  double t = tcpu(); for(int i = 0; i < N; i++){ u64 e[4]; int nl = build_nm1(ps[i], ms[i], e); sink += nl; } double tb = tcpu() - t;
  t = tcpu(); for(int i = 0; i < N; i++){ u64 e[4]; int nl = build_nm1(ps[i], ms[i], e); sink += fermat_m(ms[i], e, nl); } double tf = tcpu() - t;
  printf("BENCH build_nm1 %.1f ns; build+fermat_m (m < 2^54, n ~ 2^80) %.1f ns\n", 1e9 * tb / N, 1e9 * tf / N);
  for(int i = 0; i < N; i++){ ps[i] = 5000000000ULL + (rnd() % 5000000000ULL) | 1; ms[i] = (((u128)rnd() << 3) | ((u128)1 << 66) | 1); }
  t = tcpu(); for(int i = 0; i < N; i++){ u64 e[4]; int nl = build_nm1(ps[i], ms[i], e); sink += fermat_m(ms[i], e, nl); } tf = tcpu() - t;
  printf("BENCH build+fermat_m (m ~ 2^66, m128 path, n ~ 2^99) %.1f ns\n", 1e9 * tf / N);
  { u64 nv[4], e0[4], e1[4]; u128 mv[4]; int res[4];
    t = tcpu(); for(int i = 0; i + 4 <= N; i += 4){ for(int k = 0; k < 4; k++){ mv[k] = ms[i + k]; u128 E = (u128)ps[i + k] * mv[k] - 1; e0[k] = (u64)E; e1[k] = (u64)(E >> 64); } m128_pow2_is1_x4(mv, e0, e1, res); sink += res[0] + res[3]; } tf = tcpu() - t;
    printf("BENCH 4-lane m128 fermat (m ~ 2^66, n ~ 2^99) %.1f ns per candidate\n", 1e9 * tf / N);
    for(int i = 0; i < N; i++){ ps[i] = (rnd() % 100000000ULL) | 1; ms[i] = ((u128)(rnd() % ((u64)ps[i] * ps[i])) | 1); if(ms[i] < 15) ms[i] = 15; }
    t = tcpu(); for(int i = 0; i + 4 <= N; i += 4){ for(int k = 0; k < 4; k++){ nv[k] = (u64)ms[i + k]; u128 E = (u128)ps[i + k] * ms[i + k] - 1; e0[k] = (u64)E; e1[k] = (u64)(E >> 64); } m64_pow2_is1_x4(nv, e0, e1, res); sink += res[0] + res[3]; } tf = tcpu() - t;
    printf("BENCH 4-lane m64 fermat (m < 2^54, n ~ 2^80) %.1f ns per candidate\n", 1e9 * tf / N); }
  t = tcpu(); for(int i = 0; i < N; i++){ u128 inv; sink += invmod128(ps[i] % 99999989ULL, (u128)ps[i] * 1000003ULL, &inv); } tf = tcpu() - t;
  printf("BENCH invmod128 (64-bit modulus) %.1f ns\n", 1e9 * tf / N);
  t = tcpu(); for(int i = 0; i < N; i++){ m64 c; m64_init_light(&c, (u64)ms[i] | 1); sink += c.r1; } tf = tcpu() - t;
  printf("BENCH m64_init_light %.1f ns\n", 1e9 * tf / N);
  u8 *dead = malloc(1 << 16);
  for(u64 len = 16; len <= 32768; len *= 8){
    t = tcpu(); int reps = (int)(2000000 / len); u64 k = 0;
    for(int r = 0; r < reps; r++) k += sieveA(ps[r], ms[r] % ((u128)1 << 50), (u128)(rnd() % 1000000000ULL) * 2 + 1000000, len, dead);
    tf = tcpu() - t; printf("BENCH sieveA len=%llu: %.1f ns per candidate (%.1f us per call), killed %.3f\n", (unsigned long long)len, 1e9 * tf / (reps * (double)len), 1e6 * tf / reps, (double)k / (reps * (double)len));
  }
  printf("(sink %llu)\n", (unsigned long long)sink); free(ms); free(ps); free(dead); }


/* ================================ BPSW verdict dump (cross-check against src/bpsw.py itself) ================================ */
static void dump_bpsw(const char *fn){
  FILE *f = fopen(fn, "w"); if(!f) die("cannot write dump file");
  bctx b; bctx_init(&b); mpz_t n; mpz_init(n);
  fprintf(f, "# odd n < 1000000 passing production BPSW\n");
  for(long x = 1; x < 1000000; x += 2){ mpz_set_si(n, x); if(bpsw_mpz(&b, n)) fprintf(f, "P %ld\n", x); }
  fprintf(f, "# random odd n with production verdict (R n v); includes primes, squares, Selfridge-D-heavy n\n");
  for(long it = 0; it < 20000; it++){
    u128 x = rnd128(100) | 1 | ((u128)1 << 99); mpz_set_u128(n, x);
    if(it % 4 == 0) mpz_nextprime(n, n);
    if(it % 50 == 1){ mpz_set_u128(n, rnd128(50) | 1); mpz_mul(n, n, n); }
    if(it % 10 == 2){ mpz_set_u128(n, rnd128(40) | 1); mpz_t q; mpz_init(q); mpz_nextprime(q, n); mpz_mul_ui(n, q, 3); mpz_nextprime(q, q); mpz_mul(n, n, q); mpz_clear(q); }
    char *ns = mpz_get_str(NULL, 10, n); fprintf(f, "R %s %d\n", ns, bpsw_mpz(&b, n)); free(ns);
  }
  const long known[] = {5459, 5777, 10877, 16109, 18971, 2047, 3277, 4033, 4681, 8321};
  for(int i = 0; i < 10; i++){ mpz_set_si(n, known[i]); fprintf(f, "R %ld %d\n", known[i], bpsw_mpz(&b, n)); }
  fclose(f); bctx_clear(&b); mpz_clear(n); }

/* ================================ main ================================ */
static u64 parse_u64(const char *s){ char *end; long double v = strtold(s, &end); if(*end || v < 0 || v > 1.8e19L) die("bad number"); u64 x = (u64)v; if((long double)x != v) die("non-integer number"); return x; }
int main(int argc, char **argv){
  int do_self = 0, do_plan = 0; const char *plantfile = NULL;
  char *splits[64][3]; int nsplit = 0;
  for(int i = 1; i < argc; i++){
    const char *a = argv[i];
    #define NEED(k) do { if(i + (k) >= argc) die("missing argument"); } while(0)
    if(!strcmp(a, "--D")){ NEED(1); gD = atoll(argv[++i]); }
    else if(!strcmp(a, "--X")){ NEED(1); gX = parse_u64(argv[++i]); }
    else if(!strcmp(a, "--mode")){ NEED(1); const char *m = argv[++i]; int f = 0; for(int k = 0; k < 5; k++) if(!strcmp(m, MODENAME[k])){ gMode = k; f = 1; } if(!f) die("bad mode"); }
    else if(!strcmp(a, "--cap")){ NEED(1); gCap = atoll(argv[++i]); if(gCap < -1) die("bad cap"); }
    else if(!strcmp(a, "--split")){ NEED(3); if(nsplit == 64) die("too many splits"); splits[nsplit][0] = argv[i + 1]; splits[nsplit][1] = argv[i + 2]; splits[nsplit][2] = argv[i + 3]; nsplit++; i += 3; }
    else if(!strcmp(a, "--plants")){ NEED(1); plantfile = argv[++i]; }
    else if(!strcmp(a, "--threads")){ NEED(1); gThreads = atoi(argv[++i]); if(gThreads < 1) gThreads = 1; }
    else if(!strcmp(a, "--nosieve")) gSieve = 0;
    else if(!strcmp(a, "--sieveB")){ NEED(1); gSieveB = (u32)parse_u64(argv[++i]); }
    else if(!strcmp(a, "--ctest")){ NEED(1); gCtest = atof(argv[++i]); }
    else if(!strcmp(a, "--cpair")){ NEED(1); gCpair = atof(argv[++i]); }
    else if(!strcmp(a, "--chunkA")){ NEED(1); gChunkA = parse_u64(argv[++i]); if(!gChunkA) die("chunkA"); }
    else if(!strcmp(a, "--chunkC")){ NEED(1); gChunkC = parse_u64(argv[++i]); if(!gChunkC) die("chunkC"); }
    else if(!strcmp(a, "--progress")) gProgress = 1;
    else if(!strcmp(a, "--planonly")) do_plan = 1;
    else if(!strcmp(a, "--dumpbpsw")){ NEED(1); dump_bpsw(argv[++i]); return 0; }
    else if(!strcmp(a, "--ctestW")){ NEED(1); gCtestW = atof(argv[++i]); }
    else if(!strcmp(a, "--sample")){ NEED(1); char *c = strchr(argv[++i], ':'); if(!c) die("--sample lo:hi"); *c = 0; gX0 = parse_u64(argv[i]); gX = parse_u64(c + 1); gCap = -1; }
    else if(!strcmp(a, "--quiet")) gQuiet = 1;
    else if(!strcmp(a, "--selftest")) do_self = 1;
    else if(!strcmp(a, "--bench")) do_self = 2;
    else { fprintf(stderr, "unknown option %s\n", a); return 1; }
  }
  omp_set_num_threads(gThreads);
  if(do_self == 2){ bench(); return 0; }
  if(do_self) return selftest();
  if(!gX) die("--X required");
  if(gX0 > 3 && (gMode != MODE_BPSW || gCap != -1)) die("--sample needs bpsw mode with every p3 in branch A (table holds only the window)");
  if(gX0 > 3) printf("SAMPLE window: p3 in [%llu, %llu], branch A only; Cprimes/Cpairs meaningless here\n", (unsigned long long)gX0, (unsigned long long)gX);
  if(((1 - gD) % 4) != 0) die("D must be 1 mod 4");
  { int64_t r = (int64_t)sqrtl((long double)llabs(gD)); while(r * r > llabs(gD)) r--; while((r + 1) * (r + 1) <= llabs(gD)) r++; if(gD > 0 && r * r == gD) die("D is a square"); }
  gQ = (1 - gD) / 4;
  if(gX > ((u64)1 << 36)) die("X too large for this build (p3^3 S bookkeeping assumes X < 2^36)");
  if(gMode != MODE_BPSW && gX > ((u64)1 << 31)) die("control modes need X < 2^31 (u64 factoring of m)");
  if(gMode == MODE_WEAK && gX > 10000000ULL) fprintf(stderr, "warning: weak mode above 1e7 (spec allows <= 1e7)\n");
  init_sieve_primes(gSieveB);
  printf("anchor3 D=%lld Q=%lld X=%llu mode=%s threads=%d sieve=%s sieveB=%u ctest=%.3f ctestW=%.3f cpair=%.3f chunkA=%llu chunkC=%llu\n", (long long)gD, (long long)gQ, (unsigned long long)gX,
    MODENAME[gMode], gThreads, gSieve ? "on" : "off", gSieveB, gCtest, gCtestW, gCpair, (unsigned long long)gChunkA, (unsigned long long)gChunkC);
  double t0 = wall(), tc0 = pcpu(); build_table(gX); double tt = wall() - t0, tc = pcpu() - tc0;
  printf("TABLE D=%lld X=%llu inert=%zu overflow=%zu bytes=%.1fMB time=%.2fs cpu=%.2fs\n", (long long)gD, (unsigned long long)gX, NT, NTov, (16.0 * NT) / 1e6, tt, tc); fflush(stdout);
  if(do_plan){ if(!nsplit) pass_setup("1", "1", "1"); else if(!pass_setup(splits[0][0], splits[0][1], splits[0][2])) die("W2 incompatible"); plan_only(); return 0; }
  if(gMode == MODE_PLANT){
    if(!plantfile) die("--plants FILE required in plant mode");
    load_plants(plantfile);
    int np = 0; char pS[64][512], pA[64][512], pB[64][512];
    for(size_t k = 0; k < NPL; k++){ int f = -1; for(int q = 0; q < np; q++) if(!strcmp(pS[q], PL[k].S) && !strcmp(pA[q], PL[k].A) && !strcmp(pB[q], PL[k].B)) f = q;
      if(f < 0){ if(np == 64) die("too many plant passes"); strcpy(pS[np], PL[k].S); strcpy(pA[np], PL[k].A); strcpy(pB[np], PL[k].B); f = np++; } PL[k].pass = f; }
    printf("PLANTFILE %s: %zu overrides, %zu triples, %d passes\n", plantfile, NOV, NPL, np);
    if(!np){ np = 1; strcpy(pS[0], "1"); strcpy(pA[0], "1"); strcpy(pB[0], "1"); }
    for(int q = 0; q < np; q++){
      if(!pass_setup(pS[q], pA[q], pB[q])){ printf("PASS S=%s A=%s B=%s: (W2) incompatible (gcd(A,B) does not divide 2): no solutions\n", pS[q], pA[q], pB[q]); continue; }
      plants_prepare(q); double wt; ctrs c = run_pass(&wt); print_pass(&c, wt); plants_report(q);
    }
    return 0;
  }
  if(!nsplit){ splits[0][0] = "1"; splits[0][1] = "1"; splits[0][2] = "1"; nsplit = 1; }
  for(int q = 0; q < nsplit; q++){
    if(!pass_setup(splits[q][0], splits[q][1], splits[q][2])){ printf("PASS S=%s A=%s B=%s: (W2) incompatible (gcd(A,B) does not divide 2): no solutions\n", splits[q][0], splits[q][1], splits[q][2]); continue; }
    double wt; ctrs c = run_pass(&wt); print_pass(&c, wt);
  }
  return 0; }

/* ======================================================================================================
 * FIXES (post-review change log; appended here so the tag line numbers in the header stay valid)
 * ======================================================================================================
 * F1 (2026-10-24, review 2, low severity, confirmed by reading).  Finding: with --quiet, record_hit() printed no
 *    HIT line, so a bpsw-mode hit (a counterexample) would survive only as the count hits=N in the PASS line,
 *    without n or its factors; SPEC requires every hit to be printed with its factors.  (No production driver
 *    passes --quiet; the search itself was unaffected.)  Change: record_hit() line 668, the print block now runs
 *    when !gQuiet || gMode == MODE_BPSW, i.e. --quiet silences only the control-mode hit lines (FERMAT, LUCAS,
 *    WEAK, PLANTFOUND, W1ONE), which the self-test relies on.  Usage text line 16 updated in place, line 1 points
 *    here.  No enumeration, prefilter or test code changed; line count of lines 1-1485 unchanged.
 * Not applied (not defects of this file): review 1 found no defect; its notes concern the driver run_E.sh (check
 *    rc and one PASS / 'incompatible' line per --split), the usability of --cap -2 / --sample (optional), and a
 *    branch-C observation (no valid p1 < p2 with the true two-torus M; expected arithmetic).
 * VALIDATION OF F1 (2026-10-25, fixed binary; scripts and table: src/te/val/fixF1/; 62/62 rows PASS):
 *  - --selftest: PASSED, 0 failures.
 *  - Hit-print test (test-only variant hp2: bpsw mode with M := ord and final test strong Fermat, so bpsw-mode HIT
 *    lines occur; D = 5, X = 2e4 caps auto/0/-1, X = 1e5 caps auto/0; branch A and C print paths): fixed source,
 *    with and without --quiet: HIT lines == hits=N (35, 35, 2009, 122, 122) and identical n-sets.  Pre-fix source
 *    with --quiet: 0 HIT lines for the same hits=N (the reviewed defect reproduced); verbose n-sets identical.
 *  - Calibration bpsw D = 5, X = 1e7: hits=0, all PASS/DIAG counters identical to the pre-fix run; --quiet equal.
 *  - C1 one-torus brute-force equality 72/72, W1-final-test variant 54/54, split parts S > 1 108/108; C2 bpsw/weak
 *    vs reference 72/72; C3 planted: default cap 792/792, cap 0 642/642, cap -1 217/217 (all 23 plant runs PASS).
 *  - Every control output (c1 144, c1w 54, c6 108, c3 23 files) is identical to the pre-fix output after stripping
 *    timing fields.
 */
