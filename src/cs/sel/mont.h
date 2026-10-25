// 64-bit Montgomery arithmetic for odd moduli n < 2^63, plus helpers used by the selector-line programs.
#include <stdint.h>
typedef unsigned __int128 u128;
typedef struct { uint64_t n, ninv, r1, r2; } mont_t;   // r1 = R mod n (Montgomery 1), r2 = R^2 mod n
static inline uint64_t mredc(const mont_t *M, u128 T) {
  uint64_t m = (uint64_t)T * M->ninv;
  u128 t = (T + (u128)m * M->n) >> 64;
  uint64_t r = (uint64_t)t;
  return r >= M->n ? r - M->n : r;
}
static inline uint64_t mmul(const mont_t *M, uint64_t a, uint64_t b) { return mredc(M, (u128)a * b); }
static inline uint64_t madd(const mont_t *M, uint64_t a, uint64_t b) { uint64_t s = a + b; return s >= M->n ? s - M->n : s; }
static inline uint64_t msub(const mont_t *M, uint64_t a, uint64_t b) { return a >= b ? a - b : a + M->n - b; }
static void minit(mont_t *M, uint64_t n) {
  M->n = n; uint64_t x = n;               // Newton: x = n^-1 mod 2^64
  for (int i = 0; i < 6; i++) x *= 2 - n * x;
  M->ninv = (uint64_t)0 - x;
  M->r1 = (uint64_t)(((u128)1 << 64) % n);
  M->r2 = (uint64_t)((u128)M->r1 * M->r1 % n);
}
static inline uint64_t mto(const mont_t *M, uint64_t a) { return mmul(M, a % M->n, M->r2); }
static inline uint64_t mfrom(const mont_t *M, uint64_t a) { return mredc(M, a); }
// a^e in Montgomery form (a given in Montgomery form)
static uint64_t mpow(const mont_t *M, uint64_t a, uint64_t e) {
  uint64_t r = M->r1;
  while (e) { if (e & 1) r = mmul(M, r, a); a = mmul(M, a, a); e >>= 1; }
  return r;
}
// Fibonacci F_k mod n (returned in Montgomery form): fast doubling from the top bit
static uint64_t mfib(const mont_t *M, uint64_t k) {
  uint64_t a = 0, b = M->r1;               // (F_0, F_1)
  int top = 63; while (top >= 0 && !((k >> top) & 1)) top--;
  for (int i = top; i >= 0; i--) {
    uint64_t c = mmul(M, a, msub(M, madd(M, b, b), a));      // F_2j = F_j (2F_{j+1} - F_j)
    uint64_t d = madd(M, mmul(M, a, a), mmul(M, b, b));      // F_2j+1 = F_j^2 + F_{j+1}^2
    if ((k >> i) & 1) { a = d; b = madd(M, c, d); } else { a = c; b = d; }
  }
  return a;
}
static int isprime64(uint64_t n) {
  if (n < 2) return 0;
  static const uint64_t sp[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
  for (int i = 0; i < 12; i++) if (n % sp[i] == 0) return n == sp[i];
  mont_t M; minit(&M, n);
  uint64_t d = n - 1; int s = 0; while (!(d & 1)) { d >>= 1; s++; }
  static const uint64_t bases[] = {2, 325, 9375, 28178, 450775, 9780504, 1795265022};
  uint64_t one = M.r1, mone = M.n - M.r1;
  for (int i = 0; i < 7; i++) {
    uint64_t a = bases[i] % n; if (a == 0) continue;
    uint64_t x = mpow(&M, mto(&M, a), d);
    if (x == one || x == mone) continue;
    int ok = 0;
    for (int r = 1; r < s; r++) { x = mmul(&M, x, x); if (x == mone) { ok = 1; break; } }
    if (!ok) return 0;
  }
  return 1;
}
// divisibility by odd q via multiplicative inverse (Granlund-Montgomery)
typedef struct { uint64_t q, inv, lim; } dv_t;
static void dvinit(dv_t *D, uint64_t q) { uint64_t x = q; for (int i = 0; i < 6; i++) x *= 2 - q * x; D->q = q; D->inv = x; D->lim = UINT64_MAX / q; }
static inline int dvdiv(const dv_t *D, uint64_t *a) { uint64_t t = *a * D->inv; if (t <= D->lim) { *a = t; return 1; } return 0; }
