// Independent check of the Selfridge-reachability theorem (adversarial FX check).
// For every odd n in [3, X], not a perfect square: run Selfridge Method A (D = 5, -7, 9, -11, ...; first D with
// Jacobi(D, n) = -1), classify the selected D, check the stated conditions, and for D in {21, -27} run an independent
// strong Lucas test (P = 1, Q = (1 - D)/4) and require failure.  Own Jacobi, own Lucas (x^k mod x^2 - x + Q).
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
typedef unsigned __int128 u128;
static int jac(int64_t a, uint64_t n) { // n odd positive
  int64_t r0 = a % (int64_t)n; if (r0 < 0) r0 += n; uint64_t x = r0, y = n; int s = 1;
  while (x) { while (!(x & 1)) { x >>= 1; if ((y & 7) == 3 || (y & 7) == 5) s = -s; }
    uint64_t t = x; x = y; y = t; if ((x & 3) == 3 && (y & 3) == 3) s = -s; x %= y; }
  return y == 1 ? s : 0;
}
static int isprime_small(int k) { if (k < 2) return 0; for (int d = 2; d * d <= k; d++) if (k % d == 0) return 0; return 1; }
static uint64_t mulm(uint64_t a, uint64_t b, uint64_t m) { return (uint64_t)((u128)a * b % m); }
// x^k mod (x^2 - x + Q) over Z/n: returns (a, b) with x^k = a x + b; a = U_k, and V_k = 2b + a (P = 1)
static void xpow(uint64_t k, uint64_t Qm, uint64_t n, uint64_t *A, uint64_t *B) {
  uint64_t a = 0, b = 1; // x^0
  int top = 63; while (top >= 0 && !((k >> top) & 1)) top--;
  for (int i = top; i >= 0; i--) {
    // square: (a x + b)^2 = a^2 x^2 + 2ab x + b^2, x^2 = x - Q
    uint64_t a2 = mulm(a, a, n), ab = mulm(a, b, n), b2 = mulm(b, b, n);
    uint64_t na = (a2 + 2 * ab % n) % n; uint64_t nb = (b2 + n - mulm(a2, Qm, n)) % n;
    a = na; b = nb;
    if ((k >> i) & 1) { // multiply by x: (a x + b) x = a x^2 + b x = a (x - Q) + b x = (a + b) x - a Q
      uint64_t na2 = (a + b) % n, nb2 = (n - mulm(a, Qm, n)) % n; a = na2; b = nb2; }
  }
  *A = a; *B = b;
}
static int slucas(uint64_t n, int64_t D) {
  int64_t Q = (1 - D) / 4; int64_t qr = Q % (int64_t)n; if (qr < 0) qr += n; uint64_t Qm = qr;
  uint64_t d = n + 1; int s = 0; while (!(d & 1)) { d >>= 1; s++; }
  uint64_t a, b; xpow(d, Qm, n, &a, &b);
  if (a == 0) return 1;
  for (int r = 0; r < s; r++) { uint64_t k = d << r; xpow(k, Qm, n, &a, &b); uint64_t V = (2 * b + a) % n; if (V == 0) return 1; }
  return 0;
}
int main(int argc, char **argv) {
  uint64_t X = (uint64_t)atof(argv[1]);
  long long cnt[400] = {0}; long long viol = 0, bad21 = 0, bad27 = 0, pass_bad = 0, unexpected = 0, n21 = 0, n27 = 0, nm15 = 0;
  for (uint64_t n = 3; n <= X; n += 2) {
    uint64_t r = (uint64_t)sqrtl((long double)n); while (r * r > n) r--; while ((r + 1) * (r + 1) <= n) r++;
    if (r * r == n) continue;
    int64_t D = 5; int k;
    for (;;) { if (jac(D, n) == -1) break; D = D > 0 ? -D - 2 : -D + 2; }
    k = (int)(D < 0 ? -D : D);
    if (k < 400) cnt[k]++;
    int j3 = jac(n, 3);
    if (isprime_small(k) && k >= 5) { if (j3 == -1 && k > 27) viol++; continue; }
    if (D == -15) { nm15++; if (!(j3 == -1 && jac(n, 5) == 1 && jac(n, 7) != -1 && jac(n, 11) != -1 && jac(n, 13) != -1)) viol++; }
    else if (D == 21) { n21++; if (!(j3 == -1 && n % 5 == 0 && n % 7 != 0 && jac(n, 7) == 1)) bad21++; if (slucas(n, D)) pass_bad++; }
    else if (D == -27) { n27++; if (!(j3 == -1 && n % 35 == 0)) bad27++; if (slucas(n, D)) pass_bad++; }
    else { unexpected++; if (unexpected < 10) printf("unexpected D %lld at n %llu\n", (long long)D, (unsigned long long)n); }
  }
  printf("X=%llu viol(prime/-15)=%lld bad21=%lld bad27=%lld unexpectedD=%lld strongLucasPass(21,-27)=%lld  n(-15)=%lld n(21)=%lld n(-27)=%lld\n",
    (unsigned long long)X, viol, bad21, bad27, unexpected, pass_bad, nm15, n21, n27);
  printf("counts by |D|:"); for (int k = 5; k < 400; k++) if (cnt[k]) printf(" %d:%lld", k, cnt[k]); printf("\n");
  // sanity: the Lucas implementation on primes must pass
  int okp = 0, tot = 0; for (uint64_t p = 1000003; p < 1100000; p += 2) { int pr = 1; for (uint64_t d = 3; d * d <= p; d += 2) if (p % d == 0) { pr = 0; break; }
    if (!pr) continue; int64_t D = 5; while (jac(D, p) != -1) D = D > 0 ? -D - 2 : -D + 2; tot++; okp += slucas(p, D); }
  printf("sanity: primes in [1000003,1.1e6) passing own strong Lucas: %d / %d\n", okp, tot);
  return 0;
}
// (appended) known strong Lucas pseudoprimes sanity is done in slpsp.c
