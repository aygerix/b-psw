// Independent check of FX rules R1, R2 (D = -15 and Q = -c^2), FX3, and inert class shares.
// j1 = v2(ord_p 2); z = rank of apparition of U(1, Q) mod p computed from x^k mod (x^2 - x + Q) (coefficient of x = U_k);
// j2 = v2(z), v = v2(p - (D/p)).
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
typedef unsigned __int128 u128;
static int jac(int64_t a, uint64_t n) { int64_t r0 = a % (int64_t)n; if (r0 < 0) r0 += n; uint64_t x = r0, y = n; int s = 1;
  while (x) { while (!(x & 1)) { x >>= 1; if ((y & 7) == 3 || (y & 7) == 5) s = -s; } uint64_t t = x; x = y; y = t; if ((x & 3) == 3 && (y & 3) == 3) s = -s; x %= y; }
  return y == 1 ? s : 0; }
static uint64_t mulm(uint64_t a, uint64_t b, uint64_t m) { return (uint64_t)((u128)a * b % m); }
static uint64_t powm(uint64_t b, uint64_t e, uint64_t m) { uint64_t r = 1 % m; b %= m; while (e) { if (e & 1) r = mulm(r, b, m); b = mulm(b, b, m); e >>= 1; } return r; }
static uint64_t Uk(uint64_t k, uint64_t Qm, uint64_t n) { uint64_t a = 0, b = 1; int top = 63; while (top >= 0 && !((k >> top) & 1)) top--;
  for (int i = top; i >= 0; i--) { uint64_t a2 = mulm(a, a, n), ab = mulm(a, b, n), b2 = mulm(b, b, n); uint64_t na = (a2 + 2 * ab % n) % n; uint64_t nb = (b2 + n - mulm(a2, Qm, n)) % n; a = na; b = nb;
    if ((k >> i) & 1) { uint64_t na2 = (a + b) % n, nb2 = (n - mulm(a, Qm, n)) % n; a = na2; b = nb2; } } return a; }
static int factor(uint64_t m, uint64_t *qs) { int k = 0; for (uint64_t d = 2; d * d <= m; d++) if (m % d == 0) { qs[k++] = d; while (m % d == 0) m /= d; } if (m > 1) qs[k++] = m; return k; }
static int v2(uint64_t x) { int v = 0; while (!(x & 1)) { x >>= 1; v++; } return v; }
#define ND 24
int main(int argc, char **argv) {
  uint64_t X = (uint64_t)atof(argv[1]);
  int64_t DL[ND] = {5, -7, -11, 13, -15, 17, -19, -23, 29, -31, 37, 41, -43, -47, 53, -59, 61, -67, -71, 73, -79, -83, 89, 97};
  char *sv = calloc(X + 1, 1); for (uint64_t i = 2; i * i <= X; i++) if (!sv[i]) for (uint64_t j = i * i; j <= X; j += i) sv[j] = 1;
  long long pairs = 0, vR1 = 0, vR2m15 = 0, vR2neg = 0, vFX3 = 0, vQ2 = 0, nR2neg = 0, nQ2 = 0;
  static long long tab[ND][3][12][12]; // [D][eps+1][j1][j2] capped at 11
  long long inertJ1two[ND][3]; memset(inertJ1two, 0, sizeof inertJ1two);
  for (uint64_t p = 7; p <= X; p += 2) {
    if (sv[p]) continue;
    uint64_t qs[64];
    // j1
    uint64_t o = p - 1; int nq = factor(o, qs);
    for (int i = 0; i < nq; i++) while (o % qs[i] == 0 && powm(2, o / qs[i], p) == 1) o /= qs[i];
    int j1 = v2(o);
    for (int di = 0; di < ND; di++) {
      int64_t D = DL[di], Q = (1 - D) / 4;
      if ((D % (int64_t)p) == 0 || (Q % (int64_t)p) == 0) continue;
      int e = jac(D, p); uint64_t m = p - e; int v = v2(m);
      int64_t qr = Q % (int64_t)p; if (qr < 0) qr += p; uint64_t Qm = qr;
      if (Uk(m, Qm, p) != 0) { printf("ERROR z does not divide m p=%llu D=%lld\n", (unsigned long long)p, (long long)D); return 1; }
      uint64_t z = m; nq = factor(m, qs);
      for (int i = 0; i < nq; i++) while (z % qs[i] == 0 && Uk(z / qs[i], Qm, p) == 0) z /= qs[i];
      int j2 = v2(z); pairs++;
      int jq = jac(Q, p);
      if ((j2 == v) != (jq == -1)) vR1++;
      if (D == -15) { if (jac(10, p) == -1) { if (j2 != v - 1) vR2m15++; } else { int mx = v - 2 > 0 ? v - 2 : 0; if (j2 > mx) vR2m15++; } }
      if (D == 5 || D == 17 || D == 37) { // Q = -1, -4, -9
        if (p % 4 == 3) { nR2neg++; if (j2 != v) vR2neg++; }
        if (p % 4 == 1 && e == -1) { nR2neg++; if (j2 != 0) vR2neg++; } }
      if (e == -1 && j1 == 2) { if (p % 4 != 1 || v != 1 || ((j2 == 0) != (jq == 1))) vFX3++; inertJ1two[di][j2 < 2 ? j2 : 2]++; }
      if ((D == -7 || D == -31 || D == -71 || D == 73) && e == -1 && j1 == 2) { nQ2++; int exp = (p % 8 == 5) ? 1 : 0; if (j2 != exp) vQ2++; }
      tab[di][e + 1][j1 > 11 ? 11 : j1][j2 > 11 ? 11 : j2]++;
    }
  }
  printf("primes 7..%llu, (p,D) pairs over 24 D: %lld\n", (unsigned long long)X, pairs);
  printf("violations: R1 %lld | R2(D=-15) %lld | R2(Q=-c^2: D=5,17,37; %lld tests) %lld | FX3 %lld | Q=+-2c^2 inert j1=2 p==5 mod 8 -> (2,1), p==1 mod 8 -> (2,0) (%lld tests) %lld\n",
         vR1, vR2m15, nR2neg, vR2neg, vFX3, nQ2, vQ2);
  printf("inert j1=2 primes per D: j2=0 / j2=1 / other\n");
  for (int di = 0; di < ND; di++) printf("  D %4lld: %lld / %lld / %lld\n", (long long)DL[di], inertJ1two[di][0], inertJ1two[di][1], inertJ1two[di][2]);
  int show[2] = {0, 4};
  for (int s = 0; s < 2; s++) { int di = show[s]; long long tot = 0;
    for (int a = 0; a < 12; a++) for (int b = 0; b < 12; b++) tot += tab[di][0][a][b];
    printf("D = %lld inert (%lld primes) class shares >= 0.002:\n", (long long)DL[di], tot);
    for (int a = 0; a < 12; a++) for (int b = 0; b < 12; b++) if (tab[di][0][a][b] * 500 >= tot) printf("   (%d,%d) %lld %.4f\n", a, b, tab[di][0][a][b], (double)tab[di][0][a][b] / tot);
    printf("   (1,2) %.4f (1,0) %.4f (1,1) %.4f (0,0) %.4f (0,1) %.4f (0,2) %.4f (2,0) %.4f\n", (double)tab[di][0][1][2]/tot, (double)tab[di][0][1][0]/tot, (double)tab[di][0][1][1]/tot,
      (double)tab[di][0][0][0]/tot, (double)tab[di][0][0][1]/tot, (double)tab[di][0][0][2]/tot, (double)tab[di][0][2][0]/tot);
  }
  return 0;
}
