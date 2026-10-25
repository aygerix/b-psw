// Census of index-aware usable primes vs the Pomerance pool, D = 5, class (1,2), mod-4 side rule.
//
// For every prime p in [X0, X1) with p == 3 mod 8 and p == +-2 mod 5 (inert, class (1,2) automatic):
//   A = (p-1)/2, B = (p+1)/4 (both odd). F = primes q == 1 mod 4, q != 5; L = primes q == 3 mod 4.
//   s1(T) = F(T)-part of A, r1(T) = A/s1(T);  s2(T) = L(T)-part of B, r2(T) = B/s2(T).
//   POOL(T):   r1(T) = r2(T) = 1           (p - 1, p + 1 supported on the sides: Pomerance pool, pomgen3)
//   USABLE(T): odd(ord_p 2) | s1(T) and odd(z(p)) | s2(T)   <=>  2^(2 s1) == 1 and F_(4 s2) == 0 (mod p)
//              (the index absorbs r1 and r2: 'index stragglers' when r1 r2 > 1)
// Chebotarev/Kummer (Theorem D, exact on fibres) predicts P(USABLE | factorisation) = 1/(r1 r2).
// Output: summary (observed vs predicted, by category and by size of the absorbed part) and one record per
// USABLE(Tmax) prime: "p o1 o2 r1 r2" with o1 = odd(ord_p 2), o2 = odd(z(p)).
// Usage: ./census X0 X1 Tmax recfile
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mont.h"

#define NT 8
static const int Tlist[NT] = {50, 100, 200, 300, 500, 700, 1000, 1500};
static dv_t Fd[400], Ld[400]; static int nF, nL;
static int Fidx[NT], Lidx[NT];          // number of F (L) primes <= Tlist[j]

int main(int argc, char **argv) {
  if (argc < 5) { fprintf(stderr, "usage: census X0 X1 Tmax recfile\n"); return 1; }
  uint64_t X0 = strtoull(argv[1], 0, 10), X1 = strtoull(argv[2], 0, 10); int Tmax = atoi(argv[3]);
  FILE *rf = fopen(argv[4], "w");
  for (int q = 3; q <= Tmax; q++) {
    int pr = 1; for (int d = 2; d * d <= q; d++) if (q % d == 0) { pr = 0; break; }
    if (!pr) continue;
    if (q % 4 == 3) dvinit(&Ld[nL++], q); else if (q != 5) dvinit(&Fd[nF++], q);
  }
  for (int j = 0; j < NT; j++) {
    Fidx[j] = 0; while (Fidx[j] < nF && (int)Fd[Fidx[j]].q <= Tlist[j]) Fidx[j]++;
    Lidx[j] = 0; while (Lidx[j] < nL && (int)Ld[Lidx[j]].q <= Tlist[j]) Lidx[j]++;
  }
  int jmax = 0; while (jmax < NT && Tlist[jmax] <= Tmax) jmax++;   // thresholds used: 0..jmax-1
  // base primes for the sieve
  uint64_t sq = (uint64_t)sqrt((double)X1) + 2;
  char *bs = calloc(sq + 1, 1); uint32_t *bp = malloc(sizeof(uint32_t) * (sq / 2 + 10)); int nbp = 0;
  for (uint64_t i = 2; i <= sq; i++) if (!bs[i]) { if (i > 2) bp[nbp++] = (uint32_t)i; for (uint64_t k = i * i; k <= sq; k += i) bs[k] = 1; }
  const uint64_t SEG = 1ULL << 23;      // odd numbers per segment
  char *seg = malloc(SEG);
  long long ncand = 0, nU = 0, cat[4] = {0};      // cat: 0 pool, 1 side-1 straggler, 2 side-2, 3 both (at Tmax)
  double Ew[NT][4]; memset(Ew, 0, sizeof Ew);     // predicted counts per threshold and category
  long long Nc[NT][4]; memset(Nc, 0, sizeof Nc);  // candidates per threshold and category (pool: r1=r2=1)
  // size buckets for the absorbed part R = r1 r2 at Tmax (log10 R)
  long long bN[24] = {0}, bO[24] = {0}; double bE[24] = {0};
  long long decU[24] = {0}, decP[24] = {0};
  uint64_t lo = X0 | 1;
  for (uint64_t base = lo; base < X1; base += 2 * SEG) {
    uint64_t hi = base + 2 * SEG; if (hi > X1) hi = X1;
    uint64_t cnt = (hi - base + 1) / 2;
    memset(seg, 0, cnt);
    for (int i = 0; i < nbp; i++) {
      uint64_t q = bp[i]; if (q * q >= hi) break;
      uint64_t st = q * q; if (st < base) { st = ((base + q - 1) / q) * q; if (!(st & 1)) st += q; }
      for (uint64_t k = st; k < hi; k += 2 * q) seg[(k - base) >> 1] = 1;
    }
    for (uint64_t i = 0; i < cnt; i++) {
      if (seg[i]) continue;
      uint64_t p = base + 2 * i; if (p < 7 || p >= X1) continue;
      if ((p & 7) != 3) continue;
      int m5 = p % 5; if (m5 != 2 && m5 != 3) continue;
      ncand++;
      uint64_t A = (p - 1) >> 1, B = (p + 1) >> 2;
      // factor A over F, B over L, recording partial products at thresholds
      uint64_t a = A, s1T[NT], s2T[NT]; uint64_t s1 = 1; int j = 0;
      uint8_t eA[400]; memset(eA, 0, nF);
      for (int t = 0; t < nF; t++) {
        while (j < jmax && t == Fidx[j]) { s1T[j++] = s1; }
        while (dvdiv(&Fd[t], &a)) { s1 *= Fd[t].q; eA[t]++; }
      }
      while (j < jmax) s1T[j++] = s1;
      uint64_t b = B, s2 = 1; j = 0; uint8_t eB[400]; memset(eB, 0, nL);
      for (int t = 0; t < nL; t++) {
        while (j < jmax && t == Lidx[j]) { s2T[j++] = s2; }
        while (dvdiv(&Ld[t], &b)) { s2 *= Ld[t].q; eB[t]++; }
      }
      while (j < jmax) s2T[j++] = s2;
      for (j = 0; j < jmax; j++) {
        double r1 = (double)(A / s1T[j]), r2 = (double)(B / s2T[j]);
        int c = (r1 > 1.5) + 2 * (r2 > 1.5);
        Ew[j][c] += 1.0 / (r1 * r2); Nc[j][c]++;
      }
      uint64_t r1 = A / s1, r2 = B / s2;
      double lR = log10((double)r1 * (double)r2); int bk = (int)floor(lR + 1e-9); if (bk > 23) bk = 23;
      bN[bk]++; bE[bk] += 1.0 / ((double)r1 * (double)r2);
      mont_t M; minit(&M, p);
      uint64_t two = mto(&M, 2);
      int ok1 = (r1 == 1) || (mpow(&M, two, 2 * s1) == M.r1);
      if (!ok1) continue;
      int ok2 = (r2 == 1) || (mfib(&M, 4 * s2) == 0);
      if (!ok2) continue;
      // usable at Tmax: exact odd orders
      uint64_t o1 = s1;
      for (int t = 0; t < nF; t++) for (int e = 0; e < eA[t]; e++) {
        uint64_t q = Fd[t].q; if (mpow(&M, two, 2 * (o1 / q)) == M.r1) o1 /= q; else break;
      }
      uint64_t o2 = s2;
      for (int t = 0; t < nL; t++) for (int e = 0; e < eB[t]; e++) {
        uint64_t q = Ld[t].q; if (mfib(&M, 4 * (o2 / q)) == 0) o2 /= q; else break;
      }
      int c = (r1 > 1) + 2 * (r2 > 1);
      nU++; cat[c]++; bO[bk]++;
      int dk = (int)floor(log10((double)p)); decU[dk]++; if (c == 0) decP[dk]++;
      fprintf(rf, "%llu %llu %llu %llu %llu\n", (unsigned long long)p, (unsigned long long)o1,
              (unsigned long long)o2, (unsigned long long)r1, (unsigned long long)r2);
    }
  }
  fclose(rf);
  printf("census [%llu, %llu), Tmax = %d: candidates (inert, p == 3 mod 8) %lld\n",
         (unsigned long long)X0, (unsigned long long)X1, Tmax, ncand);
  printf("USABLE(Tmax) = %lld: pool %lld, side-1 stragglers %lld, side-2 %lld, both %lld\n", nU, cat[0], cat[1], cat[2], cat[3]);
  printf("predicted (sum 1/(r1 r2)) by threshold T: [pool, s1, s2, both] and candidates per category\n");
  for (int j = 0; j < jmax; j++)
    printf("  T=%4d: pred %.1f %.2f %.2f %.3f  | cand %lld %lld %lld %lld\n", Tlist[j], Ew[j][0], Ew[j][1], Ew[j][2], Ew[j][3],
           Nc[j][0], Nc[j][1], Nc[j][2], Nc[j][3]);
  printf("by log10 of the absorbed part R = r1 r2 at Tmax: [bucket: candidates, predicted, observed]\n");
  for (int k = 0; k < 24; k++) if (bN[k]) printf("  R in [1e%d,1e%d): %lld %.2f %lld\n", k, k + 1, bN[k], bE[k], bO[k]);
  printf("by decade of p: usable / pool:");
  for (int k = 0; k < 24; k++) if (decU[k]) printf(" 1e%d:%lld/%lld", k, decU[k], decP[k]);
  printf("\n");
  return 0;
}
