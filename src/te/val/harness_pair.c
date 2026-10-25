/* validation harness (not part of the engine): branch-C pair cost at large p2.
   Table restricted to a window [lo, hi] (as --sample does); the nsel p3 of the window with the largest
   branch-A count (smallest M(p3), i.e. typical branch-C primes) are run through runC over every p2 of the
   window below them; reports pairs, inner iterations and thread CPU per pair.
   build: gcc -O2 -fopenmp -o harness_pair harness_pair.c -lgmp -lm   (anchor3.c in the same directory) */
#define main anchor3_main
#include "anchor3.c"
#undef main
int main(int argc, char **argv){
  if(argc < 4){ fprintf(stderr, "usage: harness_pair LO HI NSEL\n"); return 1; }
  u64 lo = parse_u64(argv[1]), hi = parse_u64(argv[2]); int nsel = atoi(argv[3]);
  gD = 5; gQ = -1; gMode = MODE_BPSW; gThreads = 1; gX0 = lo; gX = hi; gCap = -2;
  init_sieve_primes(4096); build_table(hi); pass_setup("1", "1", "1"); sieve_pass_init();
  printf("window [%llu, %llu]: %zu inert primes\n", (unsigned long long)lo, (unsigned long long)hi, NT);
  size_t *sel = calloc(nsel, sizeof(size_t)); u64 *sc = calloc(nsel, sizeof(u64));
  for(size_t i = NT / 2; i < NT; i++){ progA pa; if(!setupA(i, &pa)) continue; u64 c = pa.count;
    int k = nsel - 1; if(c <= sc[k]) continue; while(k > 0 && c > sc[k - 1]){ sc[k] = sc[k - 1]; sel[k] = sel[k - 1]; k--; } sc[k] = c; sel[k] = i; }
  tstate ts; memset(&ts, 0, sizeof ts); bctx_init(&ts.b); mpz_inits(ts.n, ts.t1, NULL); ts.dead = malloc(gChunkA + 8); ts.e = malloc(64);
  double totc = 0; u64 totp = 0;
  for(int k = 0; k < nsel; k++){ if(!sc[k]) continue; size_t i3 = sel[k];
    item it; memset(&it, 0, sizeof it); it.type = 2; it.i3 = i3; it.k0 = 0; it.k1 = i3;
    ctrs before = ts.c; double t = tcpu(); runC(&ts, &it); double dt = tcpu() - t;
    u64 pr = ts.c.Cpairs - before.Cpairs, in = ts.c.Cinner - before.Cinner; totc += dt; totp += pr;
    u128 M3 = get_M(i3); char b[48];
    printf("p3=%llu M3=%s (%d bits) countA=%llu pairs=%llu inner=%llu noinv=%llu incompat=%llu cpu=%.3fs ns/pair=%.1f\n", (unsigned long long)Tp[i3], u128str(M3, b),
      M3 ? 128 - (M3 >> 64 ? __builtin_clzll((u64)(M3 >> 64)) : 64 + __builtin_clzll((u64)M3)) : 0,
      (unsigned long long)sc[k], (unsigned long long)pr, (unsigned long long)in, (unsigned long long)(ts.c.Cnoinv - before.Cnoinv), (unsigned long long)(ts.c.Cincompat - before.Cincompat), dt, 1e9 * dt / (pr ? pr : 1));
  }
  /* share of window primes with M(p) >= 2^64 */
  size_t big = 0; for(size_t i = 0; i < NT; i++) if(get_M(i) >> 64) big++;
  printf("TOTAL pairs=%llu cpu=%.3fs ns/pair=%.1f ; share of window primes with M >= 2^64: %.3f\n", (unsigned long long)totp, totc, 1e9 * totc / (totp ? totp : 1), (double)big / NT);
  return 0; }
