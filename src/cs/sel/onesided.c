// One-sided index stragglers for the Pomerance pool (D = 5, class (1,2), mod-4 side rule), beyond the census range.
// mode 1: B composed of L-primes (q == 3 mod 4, q <= T), p = 4B - 1, A = (p-1)/2 = 2B - 1 = s1 r1 (s1 = F-part).
//         p is USABLE iff p prime, 5 does not divide A, and 2^(2 s1) == 1 (mod p)   [odd(ord_p 2) | s1]
// mode 2: A composed of F-primes (q == 1 mod 4, q != 5, q <= T), p = 2A + 1, B = (p+1)/4 = (A+1)/2 = s2 r2.
//         p is USABLE iff p prime, 5 does not divide B, and F_(4 s2) == 0 (mod p)        [odd(z(p)) | s2]
// r = 1 is the Pomerance pool itself (both modes reproduce pomgen3). Kummer/Chebotarev prediction for a prime
// candidate: P(usable) = 1/r. We accumulate the prediction over prime candidates with r <= RCUT (primality is tested
// only there; the neglected tail is printed as an upper bound using the candidate count).
// Output: per decade of p: pool primes, observed stragglers, predicted stragglers (sum 1/r over prime candidates, r > 1).
// Records "p o r" (o = odd order of the smooth side's element) for every observed straggler and pool prime -> recfile.
// Usage: ./onesided mode T log10X recfile
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mont.h"

static dv_t Sd[400]; static int nS;            // primes of the OTHER side (to factor the tested side)
static uint64_t Gq[400]; static int nG;        // generating primes (smooth side)
static uint64_t Xmax, Ymax; static int mode;
static long long ncand, nprimecand, decPool[24], decObs[24]; static double decPred[24], tailBound;
static FILE *rf;
static const uint64_t RCUT = 10000000ULL;

static void test(uint64_t Y) {
  // mode 1: Y = B, p = 4B - 1, other = A = 2B - 1.  mode 2: Y = A, p = 2A + 1, other = B = (A+1)/2
  uint64_t p, other;
  if (mode == 1) { p = 4 * Y - 1; other = 2 * Y - 1; } else { p = 2 * Y + 1; other = (Y + 1) / 2; }
  if (p < 7) return;
  if (other % 5 == 0) return;                 // p == +-1 mod 5: split
  ncand++;
  uint64_t a = other, s = 1; uint8_t e[400];
  for (int t = 0; t < nS; t++) { e[t] = 0; while (dvdiv(&Sd[t], &a)) { s *= Sd[t].q; e[t]++; } }
  uint64_t r = a;                             // absorbed part (rough part of the other side)
  int dk = (int)floor(log10((double)p));
  if (r > RCUT) { tailBound += 1.0 / (double)r; }
  int isP = -1;
  if (r <= RCUT) { isP = isprime64(p); if (isP) { nprimecand++; if (r > 1) decPred[dk] += 1.0 / (double)r; } }
  mont_t M; minit(&M, p);
  int ok;
  if (mode == 1) { uint64_t two = mto(&M, 2); ok = (r == 1) || (mpow(&M, two, 2 * s) == M.r1); }
  else ok = (r == 1) || (mfib(&M, 4 * s) == 0);
  if (!ok) return;
  if (isP < 0) isP = isprime64(p);
  if (!isP) return;
  // exact odd order on the tested side
  uint64_t o = s;
  for (int t = 0; t < nS; t++) for (int k = 0; k < e[t]; k++) {
    uint64_t q = Sd[t].q; int red;
    if (mode == 1) red = (mpow(&M, mto(&M, 2), 2 * (o / q)) == M.r1); else red = (mfib(&M, 4 * (o / q)) == 0);
    if (red) o /= q; else break;
  }
  if (r == 1) decPool[dk]++; else decObs[dk]++;
  if (rf) fprintf(rf, "%llu %llu %llu\n", (unsigned long long)p, (unsigned long long)o, (unsigned long long)r);
}
static void rec(int i, uint64_t Y) {
  test(Y);
  for (int j = i; j < nG; j++) {
    uint64_t q = Gq[j];
    if (Y > Ymax / q) break;
    uint64_t y = Y;
    while (y <= Ymax / q) { y *= q; rec(j + 1, y); }
  }
}
int main(int argc, char **argv) {
  if (argc < 5) { fprintf(stderr, "usage: onesided mode T log10X recfile\n"); return 1; }
  mode = atoi(argv[1]); int T = atoi(argv[2]); double lx = atof(argv[3]); rf = fopen(argv[4], "w");
  Xmax = (uint64_t)pow(10.0, lx);
  Ymax = mode == 1 ? (Xmax + 1) / 4 : (Xmax - 1) / 2;
  for (int q = 3; q <= T; q++) {
    int pr = 1; for (int d = 2; d * d <= q; d++) if (q % d == 0) { pr = 0; break; }
    if (!pr) continue;
    int isL = (q % 4 == 3), isF = (q % 4 == 1 && q != 5);
    if ((mode == 1 && isL) || (mode == 2 && isF)) Gq[nG++] = q;
    if ((mode == 1 && isF) || (mode == 2 && isL)) dvinit(&Sd[nS++], q);
  }
  rec(0, 1);
  fclose(rf);
  printf("mode %d, T = %d, X = 1e%.1f: smooth-side candidates %lld, prime candidates with r <= %llu: %lld\n",
         mode, T, lx, ncand, (unsigned long long)RCUT, nprimecand);
  long long P = 0, O = 0; double E = 0;
  printf("decade: pool, observed stragglers, predicted stragglers (r > 1), ratio (pool+obs)/pool\n");
  for (int k = 0; k < 24; k++) if (decPool[k] || decObs[k] || decPred[k] > 0) {
    P += decPool[k]; O += decObs[k]; E += decPred[k];
    printf("  1e%-2d: %6lld %6lld %9.2f  %.3f\n", k, decPool[k], decObs[k], decPred[k],
           decPool[k] ? (double)(decPool[k] + decObs[k]) / decPool[k] : 0.0);
  }
  printf("total: pool %lld, observed stragglers %lld, predicted %.2f (+ tail <= %.3g over r > RCUT, all candidates)\n",
         P, O, E, tailBound);
  return 0;
}
