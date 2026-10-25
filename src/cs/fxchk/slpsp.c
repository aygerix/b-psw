#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>
typedef unsigned __int128 u128;
static int jac(int64_t a, uint64_t n) { int64_t r0 = a % (int64_t)n; if (r0 < 0) r0 += n; uint64_t x = r0, y = n; int s = 1;
  while (x) { while (!(x & 1)) { x >>= 1; if ((y & 7) == 3 || (y & 7) == 5) s = -s; } uint64_t t = x; x = y; y = t; if ((x & 3) == 3 && (y & 3) == 3) s = -s; x %= y; }
  return y == 1 ? s : 0; }
static uint64_t mulm(uint64_t a, uint64_t b, uint64_t m) { return (uint64_t)((u128)a * b % m); }
static void xpow(uint64_t k, uint64_t Qm, uint64_t n, uint64_t *A, uint64_t *B) { uint64_t a = 0, b = 1; int top = 63; while (top >= 0 && !((k >> top) & 1)) top--;
  for (int i = top; i >= 0; i--) { uint64_t a2 = mulm(a, a, n), ab = mulm(a, b, n), b2 = mulm(b, b, n); uint64_t na = (a2 + 2 * ab % n) % n; uint64_t nb = (b2 + n - mulm(a2, Qm, n)) % n; a = na; b = nb;
    if ((k >> i) & 1) { uint64_t na2 = (a + b) % n, nb2 = (n - mulm(a, Qm, n)) % n; a = na2; b = nb2; } } *A = a; *B = b; }
static int slucas(uint64_t n, int64_t D) { int64_t Q = (1 - D) / 4; int64_t qr = Q % (int64_t)n; if (qr < 0) qr += n; uint64_t Qm = qr;
  uint64_t d = n + 1; int s = 0; while (!(d & 1)) { d >>= 1; s++; } uint64_t a, b; xpow(d, Qm, n, &a, &b); if (a == 0) return 1;
  for (int r = 0; r < s; r++) { xpow(d << r, Qm, n, &a, &b); if ((2 * b + a) % n == 0) return 1; } return 0; }
int main(void) { // composite strong Lucas psp (Selfridge) below 1e5: OEIS A217255 begins 5459, 5777, 10877, 16109, 18971, 22499, 24569, 25199, 40309, 58519, 75077, 97439
  int c = 0; for (uint64_t n = 9; n < 100000; n += 2) { uint64_t r = (uint64_t)sqrt((double)n); while (r*r > n) r--; while ((r+1)*(r+1) <= n) r++; if (r*r == n) continue;
    int pr = 1; for (uint64_t d = 3; d * d <= n; d += 2) if (n % d == 0) { pr = 0; break; } if (pr) continue;
    int64_t D = 5; while (jac(D, n) != -1) D = D > 0 ? -D - 2 : -D + 2; if (slucas(n, D)) { printf("%llu ", (unsigned long long)n); c++; } }
  printf("\ncount %d\n", c); return 0; }
