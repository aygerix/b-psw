/* rv1.c -- adversarial review harness for src/te/anchor3.c (review only; not part of the engine).
   build: python3 mkhook.py anchor3.c anchor3_hook.c && gcc -O2 -fopenmp -o rv1 rv1.c -lgmp -lm
   usage: rv1 table  D X WIDTH [LO2 HI2]          certify the PRODUCTION segments (lo = 3 + s*2^18) covering
                                                 [X-WIDTH, X] (last one clipped at X+1) and [LO2, HI2]
          rv1 setupA D LO HI S A B [S A B]...      setupA vs independent GMP progression, every p3 of window [LO,HI]
          rv1 runA   D LO HI CAP S A B [...]       branch-A survivors (hooked run_pass) vs GMP; sieve marks vs
                                                 trial-division reference; Montgomery residues vs mpz_powm
          rv1 runC   D LO HI NSEL P3F S A B [...]  branch-C p1 sets (hooked runC) vs GMP CRT, NSEL smallest-M p3
                                                 of the window's top half (+ P3F if nonzero and in the table) */
#define main anchor3_main
#include "anchor3_hook.c"
#undef main

typedef struct { u64 p2, p3, p1; } crec;
static crec *CR = NULL; static size_t NCR = 0, CapCR = 0;
static void hookC(u64 p2, u64 p3, u64 p1){
  #pragma omp critical(hookc)
  { if(NCR == CapCR){ CapCR = CapCR ? 2 * CapCR : (1 << 16); CR = realloc(CR, sizeof(crec) * CapCR); if(!CR) die("CR"); }
    CR[NCR].p2 = p2; CR[NCR].p3 = p3; CR[NCR].p1 = p1; NCR++; } }
static void setD(long D){ gD = D; gQ = (1 - D) / 4; }
static void mpz_from_M(mpz_t z, size_t i){ u64 o = tab_ord(i), zz = tab_z(i); mpz_t a; mpz_init_set_ui(a, zz); mpz_set_ui(z, o); mpz_lcm(z, z, a); mpz_clear(a); }
static void mpz_from_getM(mpz_t z, size_t i){ mpz_set_u128(z, get_M(i)); }
/* overrides for T4 power: M(p) := a divisor of p^2 - 1 of about 2^bits (greedy over its prime factors) */
static u64 tfac(u64 x, u64 *f){ int n = 0; for(u64 q = 2; q * q <= x; q += (q == 2 ? 1 : 2)) while(x % q == 0){ f[n++] = q; x /= q; } if(x > 1) f[n++] = x; return n; }
static void make_overrides(size_t i0, size_t i1, int bits3, int bits2, size_t *sel, int ns){
  free(OV); OV = malloc(sizeof(ovr) * (i1 - i0 + 1)); NOV = 0; u64 rs0 = 0x9876543210ULL;
  for(size_t i = i0; i < i1; i++){
    u64 p = Tp[i], f[160]; int n = (int)tfac(p - 1, f); n += (int)tfac(p + 1, f + n);
    int bits = bits2; for(int k = 0; k < ns; k++) if(sel[k] == i) bits = bits3;
    for(int a = n - 1; a > 0; a--){ rs0 = rs0 * 6364136223846793005ULL + 1442695040888963407ULL; int b = (int)((rs0 >> 33) % (u64)(a + 1)); u64 t = f[a]; f[a] = f[b]; f[b] = t; }
    u128 M = 1, lim = (u128)1 << bits; for(int a = 0; a < n; a++) if(M * f[a] <= lim) M *= f[a];
    OV[NOV].p = p; OV[NOV].M = M; NOV++; }
}

/* ------------------------------------------------------------------ T1 table */
static long test_table(long D, u64 X, u64 width, u64 lo2, u64 hi2){
  setD(D); u64 seg = (u64)1 << 18; gSeg = seg; init_base_primes((hi2 > X ? hi2 : X) + 2);
  u64 los[64], his[64]; int ns = 0;
  for(u64 s = (X - width - 3) / seg; s <= (X - 3) / seg; s++){ u64 lo = 3 + s * seg, hi = lo + seg; if(hi > X + 1) hi = X + 1; los[ns] = lo; his[ns] = hi; ns++; }
  if(hi2 > lo2) for(u64 s = (lo2 - 3) / seg; s <= (hi2 - 3) / seg && ns < 64; s++){ u64 lo = 3 + s * seg; los[ns] = lo; his[ns] = lo + seg; ns++; }
  long tbad = 0, tmiss = 0, textra = 0, tcnt = 0, tov = 0; u64 aQ = (u64)labs(gQ);
  #pragma omp parallel for schedule(dynamic, 1) num_threads(gThreads) reduction(+:tbad,tmiss,textra,tcnt,tov)
  for(int k = 0; k < ns; k++){
    segws w; segws_alloc(&w); seg_process(los[k], his[k], &w);
    mpz_t z; mpz_init(z); size_t j = 0; long bad = 0, miss = 0, cnt = 0, ov = 0;
    for(u64 p = los[k]; p < his[k]; p++){
      if(!(p & 1)) continue; mpz_set_ui(z, p); if(!mpz_probab_prime_p(z, 25)) continue;
      if(mpz_si_kronecker(D, z) != -1 || aQ % p == 0) continue;
      if(j >= w.nout || w.op[j] != p){ miss++; continue; }
      if(!w.o1[j] || !w.o2[j]){ ov++; j++; continue; }
      u64 o = (p - 1) / w.o1[j], zz = (p + 1) / w.o2[j];
      if((p - 1) % w.o1[j] || (p + 1) % w.o2[j] || !certify_entry(p, o, zz)){ bad++; if(bad < 5) printf("  BAD entry D=%ld p=%llu o=%llu z=%llu\n", D, (unsigned long long)p, (unsigned long long)o, (unsigned long long)zz); }
      cnt++; j++;
    }
    long extra = (long)w.nout - (long)j;
    tbad += bad; tmiss += miss; textra += extra; tcnt += cnt; tov += ov; mpz_clear(z); segws_free(&w);
  }
  printf("T1 table D=%ld X=%llu width=%llu extra=[%llu,%llu): %d production segments, %ld inert primes certified, bad=%ld missing=%ld extra=%ld overflow=%ld\n",
    D, (unsigned long long)X, (unsigned long long)width, (unsigned long long)lo2, (unsigned long long)hi2, ns, tcnt, tbad, tmiss, textra, tov);
  fflush(stdout); return tbad + tmiss + textra + tov; }

/* ------------------------------------------------------------------ reference progression (GMP) */
/* class of m: m S == 1 (M3), m p3 S == nu (W); odd m in [15, p3^2). returns count (mpz), m0, step */
static void ref_prog(size_t i3, mpz_t cnt, mpz_t m0, mpz_t step){
  u64 p3 = Tp[i3]; mpz_t M3, R, L, t, inv, r2, g, Lg, Wg, bnd; mpz_inits(M3, R, L, t, inv, r2, g, Lg, Wg, bnd, NULL);
  mpz_set_ui(cnt, 0); mpz_set_ui(m0, 0); mpz_set_ui(step, 0);
  mpz_from_M(M3, i3);
  mpz_mod(t, PC.S, M3); if(!mpz_invert(R, t, M3)) goto out;
  mpz_set(L, M3);
  if(mpz_cmp_ui(PC.W, 1) > 0){
    mpz_mul_ui(t, PC.S, p3); mpz_mod(t, t, PC.W); if(!mpz_invert(inv, t, PC.W)) goto out;
    mpz_mul(r2, inv, PC.nu); mpz_mod(r2, r2, PC.W);
    mpz_gcd(g, L, PC.W); mpz_sub(t, r2, R); if(!mpz_divisible_p(t, g)) goto out;
    mpz_divexact(Lg, L, g); mpz_divexact(Wg, PC.W, g); mpz_divexact(t, t, g);
    if(mpz_cmp_ui(Wg, 1) > 0){ if(!mpz_invert(inv, Lg, Wg)) die("ref crt"); mpz_mul(t, t, inv); mpz_mod(t, t, Wg); } else mpz_set_ui(t, 0);
    mpz_mul(t, t, L); mpz_add(R, R, t); mpz_mul(L, L, Wg); mpz_mod(R, R, L);
    /* verify */
    mpz_mod(t, R, M3); mpz_t q; mpz_init(q); mpz_mod(q, PC.S, M3); mpz_mul(t, t, q); mpz_mod(t, t, M3); if(mpz_cmp_ui(t, 1) && mpz_cmp_ui(M3, 1)) die("ref verify W1");
    mpz_mul_ui(t, R, p3); mpz_mul(t, t, PC.S); mpz_mod(t, t, PC.W); if(mpz_cmp(t, PC.nu)) die("ref verify W2"); mpz_clear(q);
  }
  if(mpz_even_p(L)){ if(mpz_even_p(R)) goto out; mpz_set(step, L); }
  else { mpz_mul_2exp(step, L, 1); if(mpz_even_p(R)) mpz_add(R, R, L); }
  if(mpz_cmp_ui(R, 15) < 0){ mpz_ui_sub(t, 15, R); mpz_cdiv_q(t, t, step); mpz_addmul(R, t, step); }
  mpz_set(m0, R); mpz_set_ui(bnd, p3); mpz_mul_ui(bnd, bnd, p3);
  if(mpz_cmp(R, bnd) < 0){ mpz_sub(t, bnd, R); mpz_sub_ui(t, t, 1); mpz_fdiv_q(cnt, t, step); mpz_add_ui(cnt, cnt, 1); }
out: mpz_clears(M3, R, L, t, inv, r2, g, Lg, Wg, bnd, NULL); }

static int setup_passes(int np, char **v, int q){ return pass_setup(v[3 * q], v[3 * q + 1], v[3 * q + 2]); }
static int check_nu(void){ /* nu == 1 mod A, == -1 mod B, W = lcm(A,B) (independent check of pass_setup) */
  mpz_t t, l; mpz_inits(t, l, NULL); int ok = 1;
  mpz_mod(t, PC.nu, PC.A); if(mpz_cmp_ui(PC.A, 1) && mpz_cmp_ui(t, 1)) ok = 0;
  mpz_add_ui(t, PC.nu, 1); mpz_mod(t, t, PC.B); if(mpz_sgn(t)) ok = 0;
  mpz_lcm(l, PC.A, PC.B); if(mpz_cmp(l, PC.W)) ok = 0; mpz_clears(t, l, NULL); return ok; }

/* ------------------------------------------------------------------ T2 setupA */
static long test_setupA(long D, u64 lo, u64 hi, int np, char **v){
  setD(D); gMode = MODE_BPSW; gX0 = lo; gX = hi; build_table(hi); long tot = 0;
  for(int q = 0; q < np; q++){
    if(!setup_passes(np, v, q)){ printf("T2 pass S=%s A=%s B=%s incompatible\n", v[3*q], v[3*q+1], v[3*q+2]); continue; }
    if(!check_nu()){ printf("T2 nu/W WRONG for S=%s\n", PC.Ss); tot++; }
    long bad = 0, nz = 0, gmp = 0, one = 0; double sumc = 0;
    #pragma omp parallel for schedule(dynamic, 256) num_threads(gThreads) reduction(+:bad,nz,gmp,one,sumc)
    for(size_t i = 0; i < NT; i++){
      progA pa; int ok = setupA(i, &pa); mpz_t c, m0, st; mpz_inits(c, m0, st, NULL); ref_prog(i, c, m0, st);
      u64 ec = ok ? pa.count : 0; int b = 0;
      if(mpz_cmp_ui(c, ec)) b = 1;
      if(!b && ec >= 1){ mpz_t x; mpz_init(x); mpz_set_u128(x, pa.m0); if(mpz_cmp(x, m0)) b = 2; if(ec >= 2){ mpz_set_u128(x, pa.step); if(mpz_cmp(x, st)) b = 3; } mpz_clear(x); }
      if(ok && pa.step == 0) gmp++; if(ec == 1) one++;
      if(b){ bad++; if(bad < 6){ char b1[48], b2[48]; char *s1 = mpz_get_str(NULL, 10, c), *s2 = mpz_get_str(NULL, 10, m0);
        printf("  T2 MISMATCH(%d) p3=%llu ok=%d count=%llu m0=%s step=%s | ref count=%s m0=%s\n", b, (unsigned long long)Tp[i], ok, (unsigned long long)ec, u128str(pa.m0, b1), u128str(pa.step, b2), s1, s2); free(s1); free(s2); } }
      if(ec) nz++; sumc += (double)ec;
      mpz_clears(c, m0, st, NULL);
    }
    printf("T2 setupA D=%ld [%llu,%llu] S=%s A=%s B=%s Wbits=%zu Wsmall=%d: %zu p3, nonempty %ld, count==1 %ld, gmp-bigL %ld, sum count %.0f, MISMATCH %ld\n", D, (unsigned long long)lo, (unsigned long long)hi,
      PC.Ss, PC.As, PC.Bs, mpz_sizeinbase(PC.W, 2), PC.Wsmall, NT, nz, one, gmp, sumc, bad);
    fflush(stdout); tot += bad;
  }
  return tot; }

/* ------------------------------------------------------------------ T3 runA */
static u64 rv_res64(const u64 *e, int nl, u64 n){
  m64 c; m64_init_light(&c, n); int w = nl - 1; while(w > 0 && !e[w]) w--;
  u64 word = e[w]; if(!word) return 1 % n;
  int b = 63 - __builtin_clzll(word); u64 x = m64_add(c.r1, c.r1, &c);
  for(b--;;){ for(; b >= 0; b--){ x = m64_mul(x, x, &c); x = m64_add(x, x & ((u64)0 - ((word >> b) & 1)), &c); } if(--w < 0) break; word = e[w]; b = 63; }
  return m64_redc((u128)x, &c); }
static u128 rv_res128(const u64 *e, int nl, u128 n){
  m128 c; m128_init(&c, n); int w = nl - 1; while(w > 0 && !e[w]) w--;
  u64 word = e[w]; if(!word) return 1;
  int b = 63 - __builtin_clzll(word); u128 x = m128_add(c.r1, c.r1, &c);
  for(b--;;){ for(; b >= 0; b--){ x = m128_mul(x, x, &c); x = m128_add(x, x & ((u128)0 - ((word >> b) & 1)), &c); } if(--w < 0) break; word = e[w]; b = 63; }
  return m128_mul(x, 1, &c); }
typedef struct { u64 p3; u128 m; } sv;
static int sv_cmp(const void *a, const void *b){ const sv *x = a, *y = b; if(x->p3 != y->p3) return x->p3 < y->p3 ? -1 : 1; if(x->m != y->m) return x->m < y->m ? -1 : 1; return 0; }

static long test_runA(long D, u64 lo, u64 hi, long long cap, int np, char **v){
  setD(D); gMode = MODE_BPSW; gX0 = lo; gX = hi; gCap = cap; build_table(hi); init_sieve_primes(4096);
  long tot = 0;
  for(int q = 0; q < np; q++){
    if(!setup_passes(np, v, q)){ printf("T3 pass incompatible\n"); continue; }
    /* engine: hooked run_pass, survivors recorded as hits */
    gQuiet = 1; gHookSurv = 1; NHITS = 0; double wt; ctrs c = run_pass(&wt); gHookSurv = 0;
    sv *E = malloc(sizeof(sv) * (NHITS + 1)); size_t nE = 0;
    for(size_t k = 0; k < NHITS; k++) if(!HITS[k].p1){ E[nE].p3 = HITS[k].p3; E[nE].m = HITS[k].m; nE++; }
    qsort(E, nE, sizeof(sv), sv_cmp);
    /* reference */
    size_t capT = 1 << 12; sv *T = malloc(sizeof(sv) * capT); size_t nT = 0;
    long sbad = 0, wrongkill = 0, resbad = 0, e_bad = 0, lanebad = 0, forced = 0, forcedbad = 0, wide = 0; double cands = 0, sdead = 0;
    #pragma omp parallel num_threads(gThreads) reduction(+:sbad,wrongkill,resbad,e_bad,lanebad,forced,forcedbad,wide,cands,sdead)
    {
      u8 *dead = malloc(gChunkA + 8), *rdead = malloc(gChunkA + 8); mpz_t n, e, m, r, two, Sz; mpz_inits(n, e, m, r, two, Sz, NULL); mpz_set_ui(two, 2);
      u64 ev[16];
      u128 qm[2][4]; u64 q0[2][4], q1[2][4]; int qt[2][4], qn[2] = {0, 0}, res[4];
      u64 rs_local = 0x12345 + omp_get_thread_num();
      #pragma omp for schedule(dynamic, 64)
      for(size_t i3 = 0; i3 < NT; i3++){
        progA pa; if(!setupA(i3, &pa)) continue; if(!(pa.count < (u64)cap)) continue;
        u64 p3 = Tp[i3];
        for(u64 k0 = 0; k0 < pa.count; k0 += gChunkA){
          u64 len = pa.count - k0 < gChunkA ? pa.count - k0 : gChunkA; u128 ms = pa.m0 + (u128)k0 * pa.step, st = pa.step;
          int us = (len >= 4 && st);
          if(us) sdead += (double)sieveA(p3, ms, st, len, dead); else memset(dead, 0, len);
          /* reference sieve marks (same prime range rule, l | step skipped as in the engine) */
          memset(rdead, 0, len);
          if(us && !((ms >> 90) || (st >> 90))){
            int wd = ((ms + (u128)(len - 1) * st) >> 63) != 0; u64 lim = (wd ? 30 : 6) * len + 64; if(lim > gSieveB) lim = gSieveB;
            for(int s = 0; s < NSV; s++){ u32 l = SV[s].l, o = SV[s].o; if(l > lim) break; if(st % l == 0) continue;
              u64 Smo = mpz_fdiv_ui(PC.S, o);
              for(u64 d = 0; d < len; d++){ u128 md = ms + (u128)d * st; if(md % l) continue;
                u64 x = (u64)(((u128)(p3 % o) * (u64)(md % o)) % o); x = (x * Smo) % o; if(x != 1 % o) rdead[d] = 1; } }
          }
          for(u64 d = 0; d < len; d++){
            if(dead[d] != rdead[d]) sbad++;
            u128 md = ms + (u128)d * st; cands += 1; if(md >> 63) wide++;
            mpz_set_u128(m, md); make_n_m(n, p3, md); mpz_sub_ui(e, n, 1); mpz_powm(r, two, e, m); int truth = (mpz_cmp_ui(r, 1) == 0);
            if(truth && dead[d]) wrongkill++;
            if(truth){
              #pragma omp critical(tlist)
              { if(nT == capT){ capT *= 2; T = realloc(T, sizeof(sv) * capT); } T[nT].p3 = p3; T[nT].m = md; nT++; }
            }
            /* engine e = n - 1 limbs and residue via engine Montgomery primitives */
            int nl = build_nm1(p3, md, ev); { mpz_t ee; mpz_init(ee); mpz_import(ee, nl, -1, 8, 0, 0, ev); if(mpz_cmp(ee, e)) e_bad++; mpz_clear(ee); }
            if(md >> 63){ u128 x = rv_res128(ev, nl, md), w; mpz_get_u128(&w, r); if(x != w) resbad++; }
            else { u64 x = rv_res64(ev, nl, (u64)md); if(x != mpz_get_ui(r)) resbad++; }
            if(fermat_m(md, ev, nl) != truth) lanebad++;
            /* the engine's 4-lane routines on the same (m, n-1) (S = 1 only: E < 2^128) */
            if(PC.S1){ int wq = (md >> 63) ? 1 : 0; u128 E1 = (u128)p3 * md - 1;
              #define RV_PUSH(MM, EE, TT) do { qm[wq][qn[wq]] = (MM); q0[wq][qn[wq]] = (u64)(EE); q1[wq][qn[wq]] = (u64)((EE) >> 64); qt[wq][qn[wq]] = (TT); qn[wq]++; \
                if(qn[wq] == 4){ if(wq == 0){ u64 nv[4]; for(int z = 0; z < 4; z++) nv[z] = (u64)qm[0][z]; m64_pow2_is1_x4(nv, q0[0], q1[0], res); } else m128_pow2_is1_x4(qm[1], q0[1], q1[1], res); \
                  for(int z = 0; z < 4; z++){ if(qt[wq][z] == 2){ if(res[z] != 1) forcedbad++; } else if(res[z] != qt[wq][z]) lanebad++; } qn[wq] = 0; } } while(0)
              RV_PUSH(md, E1, truth);
              /* forced-1 lane: prime m, E = (m-1)*k */
              rs_local = rs_local * 6364136223846793005ULL + 1442695040888963407ULL;
              if((rs_local >> 60) == 0 && mpz_probab_prime_p(m, 15)){ u64 kk = 1 + ((rs_local >> 20) & 0x3fffffff); u128 F = (md - 1) * kk;
                if((F / kk) == md - 1){ forced++; RV_PUSH(md, F, 2); } }
            }
          }
        }
      }
      for(int w = 0; w < 2; w++) if(qn[w]){ int k = qn[w]; for(int z = k; z < 4; z++){ qm[w][z] = qm[w][0]; q0[w][z] = q0[w][0]; q1[w][z] = q1[w][0]; qt[w][z] = qt[w][0]; }
        if(w == 0){ u64 nv[4]; for(int z = 0; z < 4; z++) nv[z] = (u64)qm[0][z]; m64_pow2_is1_x4(nv, q0[0], q1[0], res); } else m128_pow2_is1_x4(qm[1], q0[1], q1[1], res);
        for(int z = 0; z < k; z++){ if(qt[w][z] == 2){ if(res[z] != 1) forcedbad++; } else if(res[z] != qt[w][z]) lanebad++; } }
      free(dead); free(rdead); mpz_clears(n, e, m, r, two, Sz, NULL);
    }
    qsort(T, nT, sizeof(sv), sv_cmp);
    long setbad = (nT != nE); for(size_t k = 0; !setbad && k < nT; k++) if(sv_cmp(&T[k], &E[k])) setbad = 1;
    printf("T3 runA D=%ld [%llu,%llu] cap=%lld S=%s A=%s B=%s: engine Aprimes=%llu Acands=%llu (ref %.0f, wide %ld) Asieved=%llu (ref %.0f) prefilter=%llu | "
           "survivors engine %zu ref %zu SETMISMATCH=%ld | sieve-mark mismatches %ld, wrong kills %ld, n-1 limb errors %ld, residue errors %ld, lane/fermat_m errors %ld, forced-1 lanes %ld bad %ld\n",
      D, (unsigned long long)lo, (unsigned long long)hi, cap, PC.Ss, PC.As, PC.Bs, (unsigned long long)c.Aprimes, (unsigned long long)c.Acands, cands, wide, (unsigned long long)c.Asieved, sdead,
      (unsigned long long)c.Apre, nE, nT, setbad, sbad, wrongkill, e_bad, resbad, lanebad, forced, forcedbad);
    for(size_t k = 0; k < nT && k < 5; k++){ char b[48]; printf("  ref survivor p3=%llu m=%s\n", (unsigned long long)T[k].p3, u128str(T[k].m, b)); }
    fflush(stdout);
    tot += setbad + sbad + wrongkill + e_bad + resbad + lanebad + forcedbad + ((double)c.Acands != cands) + ((double)c.Asieved != sdead);
    free(E); free(T);
  }
  return tot; }

/* ------------------------------------------------------------------ T4 runC */
static int crec_cmp(const void *a, const void *b){ const crec *x = a, *y = b; if(x->p3 != y->p3) return x->p3 < y->p3 ? -1 : 1; if(x->p2 != y->p2) return x->p2 < y->p2 ? -1 : 1; if(x->p1 != y->p1) return x->p1 < y->p1 ? -1 : 1; return 0; }
static long test_runC(long D, u64 lo, u64 hi, int nsel, u64 p3f, int np, char **v){
  setD(D); gMode = MODE_BPSW; gX0 = lo; gX = hi; build_table(hi); init_sieve_primes(4096); pass_setup("1", "1", "1");
  size_t sel[64]; u128 selM[64]; int ns = 0;
  for(int k = 0; k < nsel && k < 60; k++){ selM[k] = ~(u128)0; sel[k] = 0; }
  for(size_t i = NT / 2; i < NT; i++){ u128 M = get_M(i); int k = nsel - 1; if(M >= selM[k]) continue; while(k > 0 && M < selM[k - 1]){ selM[k] = selM[k - 1]; sel[k] = sel[k - 1]; k--; } selM[k] = M; sel[k] = i; }
  ns = nsel; if(p3f){ int64_t f = tab_find(p3f); if(f >= 0){ sel[ns] = (size_t)f; selM[ns] = get_M(f); ns++; } else printf("T4 forced p3 %llu not in table\n", (unsigned long long)p3f); }
  char *ovs = getenv("RV_OVR");
  if(ovs){ int b3 = atoi(ovs), b2 = atoi(strchr(ovs, ':') + 1); size_t imin = NT; for(int k = 0; k < ns; k++) if(sel[k] < imin) imin = sel[k];
    size_t imax = 0; for(int k = 0; k < ns; k++) if(sel[k] > imax) imax = sel[k];
    make_overrides(0, imax + 1, b3, b2, sel, ns); gMode = MODE_PLANT; printf("T4 overrides: %zu primes, M(p3) ~ 2^%d, M(p2) ~ 2^%d (plant-mode get_M)\n", NOV, b3, b2); }
  long tot = 0;
  for(int q = 0; q < np; q++){
    if(!setup_passes(np, v, q)){ printf("T4 pass incompatible\n"); continue; }
    sieve_pass_init();
    tstate ts; memset(&ts, 0, sizeof ts); bctx_init(&ts.b); mpz_inits(ts.n, ts.t1, NULL); ts.dead = malloc(gChunkA + 8); ts.e = malloc(8 * (PC.nSl + 8));
    for(int k = 0; k < ns; k++){
      size_t i3 = sel[k]; u64 p3 = Tp[i3];
      NCR = 0; gHookC = 1; item it; memset(&it, 0, sizeof it); it.type = 2; it.i3 = i3; it.k0 = 0; it.k1 = i3; ctrs before = ts.c; runC(&ts, &it); gHookC = 0;
      qsort(CR, NCR, sizeof(crec), crec_cmp);
      /* reference */
      mpz_t M3, M2, R, L, t, inv, r2, g, Lg, Wg, P1; mpz_inits(M3, M2, R, L, t, inv, r2, g, Lg, Wg, P1, NULL);
      mpz_from_getM(M3, i3); size_t ptr = 0; long bad = 0, nref = 0, longp = 0; u64 maxinner = 0;
      for(size_t j = 0; j < i3; j++){
        u64 p2 = Tp[j]; int none = 0;
        /* p1 p2 S == 1 (M3) */
        mpz_mul_ui(t, PC.S, p2); mpz_mod(t, t, M3); if(!mpz_invert(R, t, M3)) none = 1; mpz_set(L, M3);
        /* p1 p3 S == 1 (M2) */
        if(!none){ mpz_from_getM(M2, j); mpz_mul_ui(t, PC.S, p3); mpz_mod(t, t, M2); if(!mpz_invert(r2, t, M2)) none = 1;
          else { mpz_gcd(g, L, M2); mpz_sub(t, r2, R); if(!mpz_divisible_p(t, g)) none = 1;
            else { mpz_divexact(Lg, L, g); mpz_divexact(Wg, M2, g); mpz_divexact(t, t, g); if(mpz_cmp_ui(Wg, 1) > 0){ mpz_invert(inv, Lg, Wg); mpz_mul(t, t, inv); mpz_mod(t, t, Wg); } else mpz_set_ui(t, 0);
              mpz_mul(t, t, L); mpz_add(R, R, t); mpz_mul(L, L, Wg); mpz_mod(R, R, L); } } }
        /* p1 p2 p3 S == nu (W) */
        if(!none && mpz_cmp_ui(PC.W, 1) > 0){ mpz_mul_ui(t, PC.S, p2); mpz_mul_ui(t, t, p3); mpz_mod(t, t, PC.W); if(!mpz_invert(inv, t, PC.W)) none = 1;
          else { mpz_mul(r2, inv, PC.nu); mpz_mod(r2, r2, PC.W); mpz_gcd(g, L, PC.W); mpz_sub(t, r2, R); if(!mpz_divisible_p(t, g)) none = 1;
            else { mpz_divexact(Lg, L, g); mpz_divexact(Wg, PC.W, g); mpz_divexact(t, t, g); if(mpz_cmp_ui(Wg, 1) > 0){ mpz_invert(inv, Lg, Wg); mpz_mul(t, t, inv); mpz_mod(t, t, Wg); } else mpz_set_ui(t, 0);
              mpz_mul(t, t, L); mpz_add(R, R, t); mpz_mul(L, L, Wg); mpz_mod(R, R, L); } } }
        /* enumerate p1 in [3, p2) */
        u64 inner = 0;
        if(!none){ mpz_set(P1, R); if(mpz_cmp_ui(P1, 3) < 0){ mpz_ui_sub(t, 3, P1); mpz_cdiv_q(t, t, L); mpz_addmul(P1, t, L); }
          for(; mpz_cmp_ui(P1, p2) < 0; mpz_add(P1, P1, L)){
            u64 p1 = mpz_get_ui(P1); nref++; inner++;
            if(ptr < NCR && CR[ptr].p2 == p2 && CR[ptr].p1 == p1) ptr++;
            else { bad++; if(bad < 6) printf("  T4 MISSING p3=%llu p2=%llu p1=%llu (engine rec at ptr: p2=%llu p1=%llu)\n", (unsigned long long)p3, (unsigned long long)p2, (unsigned long long)p1,
                       ptr < NCR ? (unsigned long long)CR[ptr].p2 : 0ULL, ptr < NCR ? (unsigned long long)CR[ptr].p1 : 0ULL); }
          } }
        if(inner > maxinner) maxinner = inner; if(inner > 1000000) longp++;
        while(ptr < NCR && CR[ptr].p2 == p2){ bad++; if(bad < 6) printf("  T4 EXTRA p3=%llu p2=%llu p1=%llu\n", (unsigned long long)p3, (unsigned long long)p2, (unsigned long long)CR[ptr].p1); ptr++; }
      }
      if(ptr != NCR) bad += (long)(NCR - ptr);
      char b[48];
      printf("T4 runC D=%ld p3=%llu M3=%s S=%s A=%s B=%s Wsmall=%d: pairs engine %llu (ref %zu) records engine %zu ref %ld, max inner %llu, long pairs ref %ld engine %llu, MISMATCH %ld\n", D, (unsigned long long)p3, u128str(get_M(i3), b), PC.Ss, PC.As, PC.Bs, PC.Wsmall,
        (unsigned long long)(ts.c.Cpairs - before.Cpairs), i3, NCR, nref, (unsigned long long)maxinner, longp, (unsigned long long)(ts.c.longpairs - before.longpairs), bad);
      fflush(stdout); tot += bad + ((ts.c.Cpairs - before.Cpairs) != i3);
      mpz_clears(M3, M2, R, L, t, inv, r2, g, Lg, Wg, P1, NULL);
    }
    bctx_clear(&ts.b); mpz_clears(ts.n, ts.t1, NULL); free(ts.dead); free(ts.e);
  }
  return tot; }


/* ------------------------------------------------------------------ T5 runC CRT chain (no filtering): engine u128 chain vs mpz, every pair */
static long test_crtC(long D, u64 lo, u64 hi, int nsel, int np, char **v){
  setD(D); gMode = MODE_BPSW; gX0 = lo; gX = hi; build_table(hi); pass_setup("1", "1", "1");
  size_t sel[64]; u128 selM[64]; for(int k = 0; k < nsel; k++){ selM[k] = ~(u128)0; sel[k] = 0; }
  for(size_t i = NT / 2; i < NT; i++){ u128 M = get_M(i); int k = nsel - 1; if(M >= selM[k]) continue; while(k > 0 && M < selM[k - 1]){ selM[k] = selM[k - 1]; sel[k] = sel[k - 1]; k--; } selM[k] = M; sel[k] = i; }
  long tot = 0;
  for(int q = 0; q < np; q++){
    if(!setup_passes(np, v, q)) continue;
    long bad = 0, chk = 0, big = 0;
    #pragma omp parallel for schedule(dynamic, 1) num_threads(gThreads) reduction(+:bad,chk,big)
    for(int k = 0; k < nsel; k++){
      size_t i3 = sel[k]; u64 p3 = Tp[i3]; u128 M3 = get_M(i3), s3, c3 = 0;
      mpz_t A, B, R, L, t, inv, g, Lg, Wg, x; mpz_inits(A, B, R, L, t, inv, g, Lg, Wg, x, NULL);
      if(!invmod128(S_mod(M3), M3, &s3)) goto done;
      if(!PC.W1 && PC.Wsmall){ u128 W = PC.Ws, tt = mulmod128(p3 % W, S_mod(W), W), ti; if(!invmod128(tt, W, &ti)) goto done; c3 = mulmod128(ti, PC.nus, W); }
      for(size_t j = 0; j < i3; j++){
        u64 p2 = Tp[j]; u128 r, Lu = M3, iv, M2 = get_M(j), a2;
        if(!invmod128(p2, M3, &iv)) continue;
        r = PC.S1 ? iv : mulmod128(iv, s3, M3);
        u128 tt = mulmod128(p3 % M2, S_mod(M2), M2);
        int ok2 = invmod128(tt, M2, &a2);
        /* mpz reference for the (W1 p3) x (W1 p2) CRT */
        mpz_set_u128(L, M3); mpz_mul_ui(t, PC.S, p2); mpz_mod(t, t, L); int okr = mpz_invert(R, t, L);
        mpz_set_u128(A, M2); mpz_mul_ui(t, PC.S, p3); mpz_mod(t, t, A); int oka = mpz_invert(B, t, A);
        if(!okr){ bad++; continue; }
        mpz_set_u128(x, r); if(mpz_cmp(x, R)) bad++;
        if(oka != ok2){ bad++; continue; } if(!ok2) continue;
        mpz_set_u128(x, a2); if(mpz_cmp(x, B)) bad++;
        u128 rr, LL; int st = crt128(r, Lu, a2, M2, &rr, &LL);
        mpz_gcd(g, L, A); mpz_sub(t, B, R); int comp = mpz_divisible_p(t, g);
        if((st != 0) != comp){ bad++; continue; } if(!comp) continue;
        mpz_divexact(Lg, L, g); mpz_divexact(Wg, A, g); mpz_divexact(t, t, g); if(mpz_cmp_ui(Wg, 1) > 0){ mpz_invert(inv, Lg, Wg); mpz_mul(t, t, inv); mpz_mod(t, t, Wg); } else mpz_set_ui(t, 0);
        mpz_mul(t, t, L); mpz_add(R, R, t); mpz_mul(L, L, Wg); mpz_mod(R, R, L);
        if(st == 2){ big++; if(mpz_sizeinbase(L, 2) < 127) bad++; continue; }
        chk++; mpz_set_u128(x, rr); if(mpz_cmp(x, R)) bad++; mpz_set_u128(x, LL); if(mpz_cmp(x, L)) bad++;
        if(!PC.W1 && PC.Wsmall && !(LL >> 100)){ u128 W = PC.Ws, ip; if(!invmod128(p2 % W, W, &ip)) continue; u128 r3, L3; int s3t = crt128(rr, LL, mulmod128(c3, ip, W), W, &r3, &L3);
          if(s3t == 1){ mpz_set_u128(x, r3); mpz_mul_ui(t, x, p2); mpz_mul_ui(t, t, p3); mpz_mul(t, t, PC.S); mpz_mod(t, t, PC.W); if(mpz_cmp(t, PC.nu)) bad++; mpz_set_u128(t, L3); mpz_mod(x, x, L); if(mpz_cmp(x, R)) bad++; } }
      }
    done: mpz_clears(A, B, R, L, t, inv, g, Lg, Wg, x, NULL);
    }
    printf("T5 crtC D=%ld [%llu,%llu] nsel=%d S=%s A=%s B=%s: CRT results checked %ld, lcm>=2^127 %ld, MISMATCH %ld\n", D, (unsigned long long)lo, (unsigned long long)hi, nsel, PC.Ss, PC.As, PC.Bs, chk, big, bad);
    fflush(stdout); tot += bad;
  }
  return tot; }

int main(int argc, char **argv){
  if(argc < 3){ fprintf(stderr, "see header\n"); return 1; }
  const char *t = argv[1]; long D = atol(argv[2]); gThreads = 3; long bad = 0; char *th = getenv("RV_THREADS"); if(th) gThreads = atoi(th);
  if(!strcmp(t, "table")){ u64 X = parse_u64(argv[3]), W = parse_u64(argv[4]), lo2 = argc > 6 ? parse_u64(argv[5]) : 0, hi2 = argc > 6 ? parse_u64(argv[6]) : 0; bad = test_table(D, X, W, lo2, hi2); }
  else if(!strcmp(t, "setupA")){ u64 lo = parse_u64(argv[3]), hi = parse_u64(argv[4]); bad = test_setupA(D, lo, hi, (argc - 5) / 3, argv + 5); }
  else if(!strcmp(t, "runA")){ u64 lo = parse_u64(argv[3]), hi = parse_u64(argv[4]); long long cap = atoll(argv[5]); bad = test_runA(D, lo, hi, cap, (argc - 6) / 3, argv + 6); }
  else if(!strcmp(t, "runC")){ u64 lo = parse_u64(argv[3]), hi = parse_u64(argv[4]); int ns = atoi(argv[5]); u64 pf = parse_u64(argv[6]); bad = test_runC(D, lo, hi, ns, pf, (argc - 7) / 3, argv + 7); }
  else if(!strcmp(t, "crtC")){ u64 lo = parse_u64(argv[3]), hi = parse_u64(argv[4]); int ns = atoi(argv[5]); bad = test_crtC(D, lo, hi, ns, (argc - 6) / 3, argv + 6); }
  else { fprintf(stderr, "bad test\n"); return 1; }
  printf("RESULT %s D=%ld: %s (%ld problems)\n", t, D, bad ? "FAIL" : "ok", bad); return bad ? 1 : 0; }
