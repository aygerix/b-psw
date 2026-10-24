// Exhaustive enumeration of side-consistent subsets inside each 2-adic class (input from export.py).
// A subset S is an independent set of the conflict graph. For |S| >= 2 with an odd number of inert primes
// we test the exact Lemma 1 conditions: n = prod S == 1 mod q^eF(q) for every F-side prime power,
// n == -1 mod q^eL(q) for every L-side prime power (eF, eL = max exponent over S), and
// v2(n-1) >= j1, v2(n+1) >= j2.  Hits are printed as index lists for exact re-verification.
// Usage: ./enum export_file [node_cap]
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned __int128 u128;
#define MAXQ 96
#define MAXDEPTH 512
static int k, nq, j1, j2, D;
static uint64_t mods[MAXQ];
static uint64_t pw[MAXQ][64];            // pw[q][e] = q^e (prime powers), q index < nq-1
static int qbase[MAXQ];                   // the prime q for index (0 for the 2-power slot)
static uint64_t *res;                     // res[i*nq + t] = p_i mod mods[t]
static int *eps;
static int *nF, *nL, (*F)[2*16], (*L)[2*16];
static uint64_t *conf;                    // conf[i*W + w]
static char **pstr;
static int W;
static uint64_t forb[MAXDEPTH][80];
static uint64_t cur[MAXDEPTH][MAXQ];
static int eF[MAXQ], eL[MAXQ];
static int chosen[MAXDEPTH];
static long long nodes, checks, hits, cap;
static int maxdepth, maxsize = 1 << 30;

static int noparity_g;
static void check(int depth, int par) {
  if (depth < 2 || (par == 0 && !noparity_g)) return;
  checks++;
  uint64_t *r = cur[depth];
  for (int t = 0; t < nq - 1; t++) {
    if (eF[t]) { if (r[t] % pw[t][eF[t]] != 1) return; }
    else if (eL[t]) { uint64_t m = pw[t][eL[t]]; if (r[t] % m != m - 1) return; }
  }
  uint64_t r2 = r[nq - 1];
  if (j1 > 0 && ((r2 - 1) & ((1ULL << j1) - 1))) return;
  if (j2 > 0 && ((r2 + 1) & ((1ULL << j2) - 1))) return;
  hits++;
  printf("HIT D=%d class (%d,%d) size %d :", D, j1, j2, depth);
  for (int i = 0; i < depth; i++) printf(" %s", pstr[chosen[i]]);
  printf("\n"); fflush(stdout);
}

static void rec(int start, int depth, int par) {
  nodes++;
  if (depth > maxdepth) maxdepth = depth;
  check(depth, par);
  if (cap && nodes >= cap) return;
  if (depth >= maxsize) return;
  if (depth + 1 >= MAXDEPTH) { fprintf(stderr, "depth overflow\n"); exit(1); }
  for (int i = start; i < k; i++) {
    if (forb[depth][i >> 6] >> (i & 63) & 1) continue;
    // new forbidden set: only words >= (i+1)/64 matter
    int w0 = (i + 1) >> 6;
    for (int w = w0; w < W; w++) forb[depth + 1][w] = forb[depth][w] | conf[(size_t)i * W + w];
    uint64_t *ri = res + (size_t)i * nq;
    for (int t = 0; t < nq; t++) cur[depth + 1][t] = (uint64_t)(((u128)cur[depth][t] * ri[t]) % mods[t]);
    int saveF[16], saveL[16];
    for (int a = 0; a < nF[i]; a++) { int t = F[i][2*a], e = F[i][2*a+1]; saveF[a] = eF[t]; if (e > eF[t]) eF[t] = e; }
    for (int a = 0; a < nL[i]; a++) { int t = L[i][2*a], e = L[i][2*a+1]; saveL[a] = eL[t]; if (e > eL[t]) eL[t] = e; }
    chosen[depth] = i;
    rec(i + 1, depth + 1, par ^ (eps[i] == -1));
    for (int a = nL[i] - 1; a >= 0; a--) eL[L[i][2*a]] = saveL[a];
    for (int a = nF[i] - 1; a >= 0; a--) eF[F[i][2*a]] = saveF[a];
    if (cap && nodes >= cap) return;
  }
}

int main(int argc, char **argv) {
  FILE *f = fopen(argv[1], "r"); if (!f) { perror("open"); return 1; }
  cap = argc > 2 ? atoll(argv[2]) : 0;
  if (argc > 3) maxsize = atoi(argv[3]);
  noparity_g = argc > 4 && atoi(argv[4]);
  char tag[16];
  long long totnodes = 0, totchecks = 0, tothits = 0; int incomplete = 0;
  while (fscanf(f, "%15s", tag) == 1) {
    if (strcmp(tag, "CLASS") != 0) { fprintf(stderr, "parse error\n"); return 1; }
    fscanf(f, "%d %d %d %d %d", &D, &j1, &j2, &k, &nq);
    if (nq > MAXQ) { fprintf(stderr, "nq too large %d\n", nq); return 1; }
    fscanf(f, "%15s", tag);
    for (int t = 0; t < nq; t++) fscanf(f, "%llu", (unsigned long long *)&mods[t]);
    // prime powers: recover q from mods[t] (q^emax); pw[t][e]
    for (int t = 0; t < nq - 1; t++) {
      uint64_t m = mods[t], q = 0;
      for (uint64_t c = 3; c <= m; c += 2) if (m % c == 0) { q = c; break; }
      qbase[t] = (int)q; pw[t][0] = 1;
      for (int e = 1; e < 64; e++) { pw[t][e] = pw[t][e-1] <= m / q ? pw[t][e-1] * q : pw[t][e-1]; }
    }
    W = (k + 63) / 64;
    if (W > 80) { fprintf(stderr, "class too large %d\n", k); return 1; }
    res = malloc(sizeof(uint64_t) * (size_t)k * nq);
    eps = malloc(sizeof(int) * k); nF = malloc(sizeof(int) * k); nL = malloc(sizeof(int) * k);
    F = malloc(sizeof(*F) * k); L = malloc(sizeof(*L) * k);
    conf = calloc((size_t)k * W, sizeof(uint64_t));
    pstr = malloc(sizeof(char *) * k);
    for (int i = 0; i < k; i++) {
      char buf[64];
      fscanf(f, "%15s %63s %d", tag, buf, &eps[i]);
      pstr[i] = strdup(buf);
      for (int t = 0; t < nq; t++) fscanf(f, "%llu", (unsigned long long *)&res[(size_t)i * nq + t]);
      fscanf(f, "%d", &nF[i]); if (nF[i] > 16) { fprintf(stderr, "nF\n"); return 1; }
      for (int a = 0; a < nF[i]; a++) fscanf(f, "%d %d", &F[i][2*a], &F[i][2*a+1]);
      fscanf(f, "%d", &nL[i]); if (nL[i] > 16) { fprintf(stderr, "nL\n"); return 1; }
      for (int a = 0; a < nL[i]; a++) fscanf(f, "%d %d", &L[i][2*a], &L[i][2*a+1]);
      int nc; fscanf(f, "%d", &nc);
      for (int a = 0; a < nc; a++) { int j; fscanf(f, "%d", &j); conf[(size_t)i * W + (j >> 6)] |= 1ULL << (j & 63); }
    }
    memset(forb[0], 0, sizeof(forb[0]));
    for (int t = 0; t < nq; t++) cur[0][t] = 1 % mods[t];
    memset(eF, 0, sizeof(eF)); memset(eL, 0, sizeof(eL));
    nodes = checks = hits = 0; maxdepth = 0;
    rec(0, 0, 0);
    int done = !(cap && nodes >= cap);
    printf("D=%d class (%d,%d) k=%d nq=%d: independent sets %lld (max size %d), parity-odd checks %lld, hits %lld%s\n",
           D, j1, j2, k, nq, nodes - 1, maxdepth, checks, hits, done ? "" : " INCOMPLETE");
    fflush(stdout);
    totnodes += nodes - 1; totchecks += checks; tothits += hits; incomplete |= !done;
    free(res); free(eps); free(nF); free(nL); free(F); free(L); free(conf);
    for (int i = 0; i < k; i++) free(pstr[i]);
    free(pstr);
  }
  printf("TOTAL independent sets %lld, checks %lld, hits %lld, %s\n", totnodes, totchecks, tothits,
         incomplete ? "INCOMPLETE" : "complete");
  return 0;
}
