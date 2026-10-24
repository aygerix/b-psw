// Direct generation of Pomerance's BPSW pool for D = 5 (class (1,2), mod-4 side rule), beyond CHM.
// p = 4B - 1 prime with B composed of primes q == 3 mod 4 (q <= T)          [L side: z(p) | p+1 = 4B]
//                and A = (p-1)/2 = 2B - 1 composed of primes q == 1 mod 4, q != 5, q <= T   [F side]
// Then p == 3 mod 8, (5/p) = -1, class (1,2) automatically, and every such p is side-compatible with the
// Pomerance assignment. Prints counts by size and the exponent profile needed for log2|G|.
// Usage: ./pomgen T log10Xmax [listfile]
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
typedef unsigned __int128 u128;
static int Lq[2000], nL, Fq[2000], nF;
static uint64_t Xmax, Bmax;
static long long cntB, cntA, cntP, bucket[80];
static int emaxF[20001], emaxL[20001];
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
static int expB[2000];
static void test(uint64_t B) {
  cntB++;
  uint64_t A = 2 * B - 1;
  int eA[2000]; int nz = 0;
  for (int i = 0; i < nF && A > 1; i++) { int e = 0; while (A % Fq[i] == 0) { A /= Fq[i]; e++; } eA[i] = e; if (e) nz = i + 1; }
  if (A != 1) return;
  for (int i = nz; i < nF; i++) eA[i] = 0;
  cntA++;
  uint64_t p = 4 * B - 1;
  if (!isprime64(p)) return;
  cntP++;
  bucket[(int)floor(log10((double)p))]++;
  for (int i = 0; i < nF; i++) if (eA[i] > emaxF[Fq[i]]) emaxF[Fq[i]] = eA[i];
  for (int i = 0; i < nL; i++) if (expB[i] > emaxL[Lq[i]]) emaxL[Lq[i]] = expB[i];
  if (lf) fprintf(lf, "%llu\n", (unsigned long long)p);
}
static void rec(int i, uint64_t B) {
  // B is composed of Lq[0..i-1] (with exponents expB); extend with primes Lq[j], j >= i
  test(B);
  for (int j = i; j < nL; j++) {
    uint64_t q = Lq[j];
    if (B > Bmax / q) break;
    uint64_t b = B; int e = 0;
    while (b <= Bmax / q) { b *= q; e++; expB[j] = e; rec(j + 1, b); }
    expB[j] = 0;
  }
}
int main(int argc, char **argv) {
  int T = atoi(argv[1]); double lx = atof(argv[2]);
  Xmax = (uint64_t)pow(10.0, lx); Bmax = (Xmax + 1) / 4;
  lf = argc > 3 ? fopen(argv[3], "w") : NULL;
  for (int q = 3; q <= T; q++) { int pr = 1; for (int d = 2; d * d <= q; d++) if (q % d == 0) { pr = 0; break; }
    if (!pr) continue; if (q % 4 == 3) Lq[nL++] = q; else if (q != 5) Fq[nF++] = q; }
  rec(0, 1);
  double lg = 1.0; int used = 0;
  for (int q = 3; q <= T; q++) {
    int e = emaxF[q] > emaxL[q] ? emaxF[q] : emaxL[q];
    if (e) { lg += log2((double)(q - 1)) + (e - 1) * log2((double)q); used++; }
  }
  printf("T=%d Xmax=1e%.1f: L-side primes %d, F-side primes %d; B candidates %lld, A smooth %lld, pool primes N=%lld\n",
         T, lx, nL, nF, cntB, cntA, cntP);
  printf("  moduli used %d, log2|G| (p-1, p+1 exponents) = %.1f, N/log2|G| = %.3f\n", used, lg, cntP / lg);
  printf("  by decade:"); for (int d = 0; d < 80; d++) if (bucket[d]) printf(" 1e%d:%lld", d, bucket[d]); printf("\n");
  if (lf) fclose(lf);
  return 0;
}
