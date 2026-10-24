// V2: two-sided-smooth pool for D = 5, class (1,2), WITHOUT a side rule.
// p = 4B - 1 prime with B and A = (p-1)/2 = 2B - 1 both composed of odd primes q <= T, q != 5.
// Then p == 3 (mod 8); we keep (5/p) = -1 (p == 2, 3 mod 5), which gives class (1,2) (Lemma 18).
// Output per candidate: p, then "A:" and "B:" factorizations as q^e lists (sides are assigned later).
// Usage: ./freegen T log10Xmax out.txt
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
typedef unsigned __int128 u128;
static int Q[3000], nQ, expB[3000];
static uint64_t Xmax, Bmax;
static long long cntB, cntA, cntP;
static FILE *lf;
static uint64_t mulmod(uint64_t a, uint64_t b, uint64_t m) { return (uint64_t)((u128)a * b % m); }
static uint64_t powmod(uint64_t a, uint64_t e, uint64_t m) { uint64_t r = 1; a %= m; while (e) { if (e & 1) r = mulmod(r, a, m); a = mulmod(a, a, m); e >>= 1; } return r; }
static int isprime64(uint64_t n) {
  if (n < 2) return 0;
  static const uint64_t sp[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
  for (int i = 0; i < 12; i++) { if (n % sp[i] == 0) return n == sp[i]; }
  uint64_t d = n - 1; int s = 0; while (!(d & 1)) { d >>= 1; s++; }
  static const uint64_t bases[] = {2, 325, 9375, 28178, 450775, 9780504, 1795265022};
  for (int i = 0; i < 7; i++) {
    uint64_t a = bases[i] % n; if (a == 0) continue;
    uint64_t x = powmod(a, d, n);
    if (x == 1 || x == n - 1) continue;
    int ok = 0;
    for (int r = 1; r < s; r++) { x = mulmod(x, x, n); if (x == n - 1) { ok = 1; break; } }
    if (!ok) return 0;
  }
  return 1;
}
static void test(uint64_t B) {
  cntB++;
  uint64_t A = 2 * B - 1, a = A;
  int eA[3000];
  for (int i = 0; i < nQ; i++) { int e = 0; while (a % Q[i] == 0) { a /= Q[i]; e++; } eA[i] = e; }
  if (a != 1) return;
  cntA++;
  uint64_t p = 4 * B - 1;
  if (p % 5 != 2 && p % 5 != 3) return;
  if (!isprime64(p)) return;
  cntP++;
  fprintf(lf, "%llu A:", (unsigned long long)p);
  for (int i = 0; i < nQ; i++) if (eA[i]) fprintf(lf, " %d^%d", Q[i], eA[i]);
  fprintf(lf, " B:");
  for (int i = 0; i < nQ; i++) if (expB[i]) fprintf(lf, " %d^%d", Q[i], expB[i]);
  fprintf(lf, "\n");
}
static void rec(int i, uint64_t B) {
  test(B);
  for (int j = i; j < nQ; j++) {
    uint64_t q = Q[j];
    if (B > Bmax / q) break;
    uint64_t b = B; int e = 0;
    while (b <= Bmax / q) { b *= q; e++; expB[j] = e; rec(j + 1, b); }
    expB[j] = 0;
  }
}
int main(int argc, char **argv) {
  int T = atoi(argv[1]); double lx = atof(argv[2]);
  Xmax = (uint64_t)pow(10.0, lx); Bmax = (Xmax + 1) / 4;
  lf = fopen(argv[3], "w");
  for (int q = 3; q <= T; q++) { int pr = 1; for (int d = 2; d * d <= q; d++) if (q % d == 0) { pr = 0; break; }
    if (pr && q != 5) Q[nQ++] = q; }
  rec(0, 1);
  printf("T=%d Xmax=1e%.1f: odd primes used %d; B candidates %lld, A smooth %lld, candidates p %lld\n", T, lx, nQ, cntB, cntA, cntP);
  fclose(lf);
  return 0;
}
