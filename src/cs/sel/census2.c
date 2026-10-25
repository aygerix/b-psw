// Side-rule-free census: every inert prime p < X (D = 5) whose odd orders odd(ord_p 2) and odd(z(p)) are both
// supported on odd primes <= T0 (5 excluded). Such p is usable for EVERY side partition (F, L) of the primes <= T0
// with supp(o1) in F and supp(o2) in L (o1 | p-1, o2 | p+1, so the supports are disjoint automatically).
// Records "p j1 j2 o1 o2" with (j1, j2) = (v2 ord_p 2, v2 z(p)) the 2-adic class.
// Test: s1 = T0-part of odd(p-1), s2 = T0-part of odd(p+1); usable iff 2^(2^v2(p-1) s1) == 1 and F_(2^v2(p+1) s2) == 0.
// Usage: ./census2 X T0 recfile
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mont.h"

static dv_t Pd[200]; static int nP;

int main(int argc, char **argv) {
  if (argc < 4) { fprintf(stderr, "usage: census2 X T0 recfile\n"); return 1; }
  uint64_t X = strtoull(argv[1], 0, 10); int T0 = atoi(argv[2]); FILE *rf = fopen(argv[3], "w");
  for (int q = 3; q <= T0; q++) {
    int pr = 1; for (int d = 2; d * d <= q; d++) if (q % d == 0) { pr = 0; break; }
    if (pr && q != 5) dvinit(&Pd[nP++], q);
  }
  uint64_t sq = (uint64_t)sqrt((double)X) + 2;
  char *bs = calloc(sq + 1, 1); uint32_t *bp = malloc(sizeof(uint32_t) * (sq / 2 + 10)); int nbp = 0;
  for (uint64_t i = 2; i <= sq; i++) if (!bs[i]) { if (i > 2) bp[nbp++] = (uint32_t)i; for (uint64_t k = i * i; k <= sq; k += i) bs[k] = 1; }
  const uint64_t SEG = 1ULL << 23; char *seg = malloc(SEG);
  long long ninert = 0, nrec = 0;
  for (uint64_t base = 1; base < X; base += 2 * SEG) {
    uint64_t hi = base + 2 * SEG; if (hi > X) hi = X;
    uint64_t cnt = (hi - base + 1) / 2;
    memset(seg, 0, cnt);
    for (int i = 0; i < nbp; i++) {
      uint64_t q = bp[i]; if (q * q >= hi) break;
      uint64_t st = q * q; if (st < base) { st = ((base + q - 1) / q) * q; if (!(st & 1)) st += q; }
      for (uint64_t k = st; k < hi; k += 2 * q) seg[(k - base) >> 1] = 1;
    }
    for (uint64_t i = 0; i < cnt; i++) {
      if (seg[i]) continue;
      uint64_t p = base + 2 * i; if (p < 7 || p >= X) continue;
      int m5 = p % 5; if (m5 != 2 && m5 != 3) continue;
      ninert++;
      uint64_t A = p - 1, B = p + 1; int a2 = 0, b2 = 0;
      while (!(A & 1)) { A >>= 1; a2++; } while (!(B & 1)) { B >>= 1; b2++; }
      uint64_t a = A, s1 = 1; uint8_t eA[200];
      for (int t = 0; t < nP; t++) { eA[t] = 0; while (dvdiv(&Pd[t], &a)) { s1 *= Pd[t].q; eA[t]++; } }
      mont_t M; minit(&M, p); uint64_t two = mto(&M, 2);
      if (a != 1 && mpow(&M, two, s1 << a2) != M.r1) continue;
      uint64_t b = B, s2 = 1; uint8_t eB[200];
      for (int t = 0; t < nP; t++) { eB[t] = 0; while (dvdiv(&Pd[t], &b)) { s2 *= Pd[t].q; eB[t]++; } }
      if (b != 1 && mfib(&M, s2 << b2) != 0) continue;
      // exact odd orders and 2-adic classes
      uint64_t o1 = s1;
      for (int t = 0; t < nP; t++) for (int e = 0; e < eA[t]; e++) {
        uint64_t q = Pd[t].q; if (mpow(&M, two, (o1 / q) << a2) == M.r1) o1 /= q; else break;
      }
      int j1 = a2; while (j1 > 0 && mpow(&M, two, o1 << (j1 - 1)) == M.r1) j1--;
      uint64_t o2 = s2;
      for (int t = 0; t < nP; t++) for (int e = 0; e < eB[t]; e++) {
        uint64_t q = Pd[t].q; if (mfib(&M, (o2 / q) << b2) == 0) o2 /= q; else break;
      }
      int j2 = b2; while (j2 > 0 && mfib(&M, o2 << (j2 - 1)) == 0) j2--;
      nrec++;
      fprintf(rf, "%llu %d %d %llu %llu\n", (unsigned long long)p, j1, j2, (unsigned long long)o1, (unsigned long long)o2);
    }
  }
  fclose(rf);
  printf("census2 X = %llu, T0 = %d: inert primes %lld, records (both odd orders T0-supported) %lld\n",
         (unsigned long long)X, T0, ninert, nrec);
  return 0;
}
