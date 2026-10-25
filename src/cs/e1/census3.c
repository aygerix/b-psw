// E1 pool census, all ten Selfridge D at once: every prime 3 <= p < X (p not dividing D*Q) whose two odd orders
//   o1 = odd(ord_p 2)            (ord_p 2 | p - 1)
//   o2 = odd(z_D(p))             (rank of apparition in U(P=1, Q=(1-D)/4); z | p - (D/p))
// are both supported on the odd primes <= T0 (5 INCLUDED as a row; for D = 5 it never occurs on inert p).
// Inert (D/p) = -1 and split (D/p) = +1 primes are both recorded; a split prime is kept only when gcd(o1, o2) = 1
// (otherwise a row q | gcd would need n == +1 and n == -1 mod q).
// The Fermat test is D-independent and is done once per p; the Lucas test only on survivors.
// Record line: "D p t j1 j2 o1 o2"   t = 0 inert, 1 split;  (j1, j2) = (v2 ord_p 2, v2 z_D(p)).
// Usage: ./census3 X T0 recfile [D1 D2 ...]   (default D list: 5 -7 -11 13 -15 17 -19 -23 29 -31)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../sel/mont.h"

static dv_t Pd[200]; static int nP;

static int jacobi(int64_t a, uint64_t n) {        // n odd > 0
  int64_t am = a % (int64_t)n; if (am < 0) am += (int64_t)n;
  uint64_t x = (uint64_t)am, y = n; int r = 1;
  while (x) {
    while (!(x & 1)) { x >>= 1; uint64_t m8 = y & 7; if (m8 == 3 || m8 == 5) r = -r; }
    uint64_t t = x; x = y; y = t;
    if ((x & 3) == 3 && (y & 3) == 3) r = -r;
    x %= y;
  }
  return y == 1 ? r : 0;
}

// U_k(P = 1, Q) mod n, Montgomery form; Qm = Q in Montgomery form. (a, b) = (U_j, U_{j+1}).
// U_2j = U_j (2 U_{j+1} - U_j);  U_2j+1 = U_{j+1}^2 - Q U_j^2;  U_2j+2 = U_2j+1 - Q U_2j.
static uint64_t mlucU(const mont_t *M, uint64_t Qm, uint64_t k) {
  uint64_t a = 0, b = M->r1;
  int top = 63; while (top >= 0 && !((k >> top) & 1)) top--;
  for (int i = top; i >= 0; i--) {
    uint64_t c = mmul(M, a, msub(M, madd(M, b, b), a));
    uint64_t d = msub(M, mmul(M, b, b), mmul(M, Qm, mmul(M, a, a)));
    if ((k >> i) & 1) { a = d; b = msub(M, d, mmul(M, Qm, c)); } else { a = c; b = d; }
  }
  return a;
}

int main(int argc, char **argv) {
  if (argc < 4) { fprintf(stderr, "usage: census3 X T0 recfile [D ...]\n"); return 1; }
  uint64_t X = (uint64_t)strtod(argv[1], 0); int T0 = atoi(argv[2]); FILE *rf = fopen(argv[3], "w");
  int64_t Dl[32]; int nD = 0;
  if (argc > 4) { for (int i = 4; i < argc && nD < 32; i++) Dl[nD++] = atoll(argv[i]); }
  else { int64_t d0[10] = {5, -7, -11, 13, -15, 17, -19, -23, 29, -31}; for (int i = 0; i < 10; i++) Dl[nD++] = d0[i]; }
  int64_t Ql[32]; for (int i = 0; i < nD; i++) Ql[i] = (1 - Dl[i]) / 4;
  for (int q = 3; q <= T0; q++) {
    int pr = 1; for (int d = 2; d * d <= q; d++) if (q % d == 0) { pr = 0; break; }
    if (pr) dvinit(&Pd[nP++], q);
  }
  uint64_t sq = (uint64_t)sqrt((double)X) + 2;
  char *bs = calloc(sq + 1, 1); uint32_t *bp = malloc(sizeof(uint32_t) * (sq / 2 + 10)); int nbp = 0;
  for (uint64_t i = 2; i <= sq; i++) if (!bs[i]) { if (i > 2) bp[nbp++] = (uint32_t)i; for (uint64_t k = i * i; k <= sq; k += i) bs[k] = 1; }
  const uint64_t SEG = 1ULL << 23; char *seg = malloc(SEG);
  long long nprime = 0, nferm = 0, nrec[32] = {0}, nsplitrec[32] = {0};
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
      uint64_t p = base + 2 * i; if (p < 3 || p >= X) continue;
      nprime++;
      // Fermat side (D-independent)
      uint64_t A = p - 1; int a2 = 0; while (!(A & 1)) { A >>= 1; a2++; }
      uint64_t a = A, s1 = 1; uint8_t eA[200];
      for (int t = 0; t < nP; t++) { eA[t] = 0; while (dvdiv(&Pd[t], &a)) { s1 *= Pd[t].q; eA[t]++; } }
      mont_t M; minit(&M, p); uint64_t two = mto(&M, 2);
      if (a != 1 && mpow(&M, two, s1 << a2) != M.r1) continue;
      nferm++;
      uint64_t o1 = s1;
      for (int t = 0; t < nP; t++) for (int e = 0; e < eA[t]; e++) {
        uint64_t q = Pd[t].q; if (mpow(&M, two, (o1 / q) << a2) == M.r1) o1 /= q; else break;
      }
      int j1 = a2; while (j1 > 0 && mpow(&M, two, o1 << (j1 - 1)) == M.r1) j1--;
      // Lucas side, per D
      for (int di = 0; di < nD; di++) {
        int64_t D = Dl[di], Q = Ql[di];
        int64_t Da = D < 0 ? -D : D, Qa = Q < 0 ? -Q : Q;
        if ((uint64_t)Da % p == 0 || (Qa != 0 && (uint64_t)Qa % p == 0)) continue;
        int js = jacobi(D, p); if (js == 0) continue;
        uint64_t m = js == -1 ? p + 1 : p - 1;
        uint64_t B = m; int b2 = 0; while (!(B & 1)) { B >>= 1; b2++; }
        uint64_t b = B, s2 = 1; uint8_t eB[200];
        for (int t = 0; t < nP; t++) { eB[t] = 0; while (dvdiv(&Pd[t], &b)) { s2 *= Pd[t].q; eB[t]++; } }
        int64_t Qr = Q % (int64_t)p; if (Qr < 0) Qr += (int64_t)p;
        uint64_t Qm = mto(&M, (uint64_t)Qr);
        if (b != 1 && mlucU(&M, Qm, s2 << b2) != 0) continue;
        uint64_t o2 = s2;
        for (int t = 0; t < nP; t++) for (int e = 0; e < eB[t]; e++) {
          uint64_t q = Pd[t].q; if (mlucU(&M, Qm, (o2 / q) << b2) == 0) o2 /= q; else break;
        }
        int j2 = b2; while (j2 > 0 && mlucU(&M, Qm, o2 << (j2 - 1)) == 0) j2--;
        if (js == 1) {          // split: both orders divide p - 1; need coprime odd parts
          uint64_t g = o1, h = o2; while (h) { uint64_t t = g % h; g = h; h = t; }
          if (g != 1) continue;
          nsplitrec[di]++;
        }
        nrec[di]++;
        fprintf(rf, "%lld %llu %d %d %d %llu %llu\n", (long long)D, (unsigned long long)p, js == 1, j1, j2,
                (unsigned long long)o1, (unsigned long long)o2);
      }
    }
  }
  fclose(rf);
  printf("census3 X = %llu, T0 = %d: primes %lld, Fermat-smooth %lld\n", (unsigned long long)X, T0, nprime, nferm);
  for (int di = 0; di < nD; di++) printf("  D = %lld: records %lld (split %lld)\n", (long long)Dl[di], nrec[di], nsplitrec[di]);
  return 0;
}
