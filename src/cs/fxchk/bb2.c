// bb2: independent re-implementation (adversarial FX check) of an exact upper bound for
//   M(lambda) = max over side-consistent configurations C of  N(C) - lambda * bits(C)
// on a bbc instance file.  Written from scratch: own parser, own Dinic, own branching rule (product of side weights),
// decision mode with a fixed threshold T (prune when the relaxation bound <= T).
// Relaxation at a node: maximum-weight closure over alive primes with every free row offered on both sides
// (prime weight 1, item weight -lambda*cost, costs rounded DOWN to multiples of 1/KS: bound >= true relaxed value).
// Deletion: primes outside the MAXIMAL optimal closure are dropped in the subtree.  Justification (my own): the closure
// value f is supermodular on prime sets, so for any feasible C in the subtree f(C n Smax) >= f(C) + f(Smax) - f(C u Smax)
// >= f(C), and C n Smax is feasible.  Leaf: minimal optimal closure side-consistent -> exact.
// Output: M_ub = max(T, max leaf bound) and the best exact configuration found (N, bits, value).
// Usage: ./bb2 INST lam T [nodecap]
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#define KS (1LL << 24)
#define BIG (1LL << 61)
static int nP, nI, nR, *Irow, *Iside, *Ilev, *Ipar; static double *Icost;
static long long *P; static int *Pk, **Pit;
static int **needS[2], *needC[2];
// graph
static int NV, NE, *hd, *nx, *tt, *lv, *it, *qq; static long long *cp;
static void adde(int u, int v, long long c) { tt[NE] = v; cp[NE] = c; nx[NE] = hd[u]; hd[u] = NE++; tt[NE] = u; cp[NE] = 0; nx[NE] = hd[v]; hd[v] = NE++; }
static int bfs(void) { for (int i = 0; i < NV; i++) lv[i] = -1; int h = 0, t = 0; qq[t++] = 0; lv[0] = 0;
  while (h < t) { int u = qq[h++]; for (int e = hd[u]; e >= 0; e = nx[e]) if (cp[e] > 0 && lv[tt[e]] < 0) { lv[tt[e]] = lv[u] + 1; qq[t++] = tt[e]; } }
  return lv[1] >= 0; }
static long long aug(int u, long long f) { if (u == 1) return f;
  for (; it[u] >= 0; it[u] = nx[it[u]]) { int e = it[u], v = tt[e];
    if (cp[e] > 0 && lv[v] == lv[u] + 1) { long long d = aug(v, f < cp[e] ? f : cp[e]); if (d > 0) { cp[e] -= d; cp[e ^ 1] += d; return d; } } }
  return 0; }
static int *pnode, *inode, *ilist; static uint8_t *rs, *rt, *smin, *smax, *imn, *imx;
static double relax(const uint8_t *alive, double lam) {
  int na = 0; for (int i = 0; i < nP; i++) pnode[i] = alive[i] ? 2 + na++ : -1;
  int ni = 0; for (int j = 0; j < nI; j++) inode[j] = -1;
  for (int i = 0; i < nP; i++) if (alive[i]) for (int r = 0; r < Pk[i]; r++) for (int j = Pit[i][r]; j >= 0 && inode[j] == -1; j = Ipar[j]) inode[j] = -2;
  for (int j = 0; j < nI; j++) if (inode[j] == -2) { inode[j] = 2 + na + ni; ilist[ni++] = j; }
  NV = 2 + na + ni; NE = 0; for (int v = 0; v < NV; v++) hd[v] = -1;
  for (int x = 0; x < ni; x++) { int j = ilist[x]; long long c = (long long)floor(lam * Icost[j] * (double)KS); if (c > 0) adde(inode[j], 1, c); if (Ipar[j] >= 0) adde(inode[j], inode[Ipar[j]], BIG); }
  for (int i = 0; i < nP; i++) if (alive[i]) { adde(0, pnode[i], KS); for (int r = 0; r < Pk[i]; r++) adde(pnode[i], inode[Pit[i][r]], BIG); }
  long long flow = 0; while (bfs()) { for (int v = 0; v < NV; v++) it[v] = hd[v]; long long f; while ((f = aug(0, BIG)) > 0) flow += f; }
  // residual reachability
  memset(rs, 0, NV); memset(rt, 0, NV); int h = 0, t = 0; qq[t++] = 0; rs[0] = 1;
  while (h < t) { int u = qq[h++]; for (int e = hd[u]; e >= 0; e = nx[e]) if (cp[e] > 0 && !rs[tt[e]]) { rs[tt[e]] = 1; qq[t++] = tt[e]; } }
  h = t = 0; qq[t++] = 1; rt[1] = 1;   // nodes u with a residual path u -> T: reverse edges; edge e: x -> tt[e]; residual u->v on edge e from u
  while (h < t) { int v = qq[h++]; for (int e = hd[v]; e >= 0; e = nx[e]) { int u = tt[e]; /* edge e^1 is u -> v */ if (cp[e ^ 1] > 0 && !rt[u]) { rt[u] = 1; qq[t++] = u; } } }
  for (int i = 0; i < nP; i++) { smin[i] = alive[i] && rs[pnode[i]]; smax[i] = alive[i] && !rt[pnode[i]]; }
  for (int j = 0; j < nI; j++) { imn[j] = inode[j] >= 0 && rs[inode[j]]; imx[j] = inode[j] >= 0 && !rt[inode[j]]; }
  return ((double)na * KS - (double)flow) / (double)KS;
}
static double value_of(const uint8_t *sel, double lam, int *N, double *bits, int *consistent) {
  static uint8_t *u = 0; if (!u) u = malloc(nI); memset(u, 0, nI); int n = 0;
  for (int i = 0; i < nP; i++) if (sel[i]) { n++; for (int r = 0; r < Pk[i]; r++) for (int j = Pit[i][r]; j >= 0 && !u[j]; j = Ipar[j]) u[j] = 1; }
  double b = 0; int *sd = malloc(nR * sizeof(int)); for (int r = 0; r < nR; r++) sd[r] = -1; int ok = 1;
  for (int j = 0; j < nI; j++) if (u[j]) { b += Icost[j]; if (Ilev[j] == 1) { if (sd[Irow[j]] >= 0 && sd[Irow[j]] != Iside[j]) ok = 0; sd[Irow[j]] = Iside[j]; } }
  free(sd); *N = n; *bits = b; *consistent = ok; return n - lam * b;
}
typedef struct { uint8_t *alive; } nd_t;
int main(int argc, char **argv) {
  FILE *f = fopen(argv[1], "r"); double lam = atof(argv[2]), T = atof(argv[3]); long long cap = argc > 4 ? atoll(argv[4]) : 2000000000LL;
  if (fscanf(f, "%d %d %d", &nP, &nI, &nR) != 3) return 1;
  Irow = malloc(nI * 4); Iside = malloc(nI * 4); Ilev = malloc(nI * 4); Ipar = malloc(nI * 4); Icost = malloc(nI * 8);
  for (int j = 0; j < nI; j++) if (fscanf(f, "%d %d %d %d %lf", &Irow[j], &Iside[j], &Ilev[j], &Ipar[j], &Icost[j]) != 5) return 1;
  P = malloc(nP * 8); Pk = malloc(nP * 4); Pit = malloc(nP * sizeof(int *)); long long ne = 0;
  for (int i = 0; i < nP; i++) { if (fscanf(f, "%lld %d", &P[i], &Pk[i]) != 2) return 1; Pit[i] = malloc((Pk[i] + 1) * 4); ne += Pk[i] + 1;
    for (int r = 0; r < Pk[i]; r++) if (fscanf(f, "%d", &Pit[i][r]) != 1) return 1; }
  fclose(f);
  for (int s = 0; s < 2; s++) { needC[s] = calloc(nR, 4); needS[s] = malloc(nR * sizeof(int *)); }
  for (int i = 0; i < nP; i++) for (int r = 0; r < Pk[i]; r++) { int j = Pit[i][r]; needC[Iside[j]][Irow[j]]++; }
  for (int s = 0; s < 2; s++) for (int r = 0; r < nR; r++) { needS[s][r] = malloc((needC[s][r] + 1) * 4); needC[s][r] = 0; }
  for (int i = 0; i < nP; i++) for (int r = 0; r < Pk[i]; r++) { int j = Pit[i][r]; needS[Iside[j]][Irow[j]][needC[Iside[j]][Irow[j]]++] = i; }
  int nvm = 2 + nP + nI; long long nem = 2 * (ne + 2 * nI + nP + 8);
  hd = malloc(nvm * 4); nx = malloc(nem * 4); tt = malloc(nem * 4); cp = malloc(nem * 8); lv = malloc(nvm * 4); it = malloc(nvm * 4); qq = malloc(nvm * 4);
  pnode = malloc(nP * 4); inode = malloc(nI * 4); ilist = malloc(nI * 4); rs = malloc(nvm); rt = malloc(nvm); smin = malloc(nP); smax = malloc(nP); imn = malloc(nI); imx = malloc(nI);
  uint8_t *uF = malloc(nR), *uL = malloc(nR), *fixedr = malloc(nR);
  int capst = 1024, sp = 0; nd_t *st = malloc(capst * sizeof(nd_t));
  st[sp].alive = malloc(nP); memset(st[sp].alive, 1, nP); sp++;
  double Mub = T, best = -1e300, bestbits = 0; int bestN = 0, bestcons = 1; long long nodes = 0, leaves_exact = 0, pruned = 0; int status = 0;
  clock_t t0 = clock();
  while (sp > 0) {
    nd_t nd = st[--sp];
    if (++nodes > cap) { status = 1; free(nd.alive); break; }
    double ub = relax(nd.alive, lam);
    if (ub <= T + 1e-9) { pruned++; free(nd.alive); continue; }
    // rows used on both sides by the minimal closure (level-1 items)
    memset(uF, 0, nR); memset(uL, 0, nR);
    for (int j = 0; j < nI; j++) if (imn[j] && Ilev[j] == 1) { if (Iside[j]) uL[Irow[j]] = 1; else uF[Irow[j]] = 1; }
    int conf = 0; for (int r = 0; r < nR; r++) if (uF[r] && uL[r]) { conf = 1; break; }
    if (!conf) { int N, cons; double bits; double v = value_of(smin, lam, &N, &bits, &cons); leaves_exact++;
      if (ub > Mub) Mub = ub; if (v > best) { best = v; bestN = N; bestbits = bits; bestcons = cons; } free(nd.alive); continue; }
    memset(uF, 0, nR); memset(uL, 0, nR);
    for (int j = 0; j < nI; j++) if (imx[j] && Ilev[j] == 1) { if (Iside[j]) uL[Irow[j]] = 1; else uF[Irow[j]] = 1; }
    int br = -1; double bw = -1;
    for (int r = 0; r < nR; r++) if (uF[r] && uL[r]) { double wF = 0, wL = 0;
      for (int x = 0; x < needC[0][r]; x++) wF += smax[needS[0][r][x]];
      for (int x = 0; x < needC[1][r]; x++) wL += smax[needS[1][r][x]];
      double w = wF * wL; if (w > bw) { bw = w; br = r; } }
    if (br < 0) { fprintf(stderr, "no branch row\n"); return 2; }
    for (int s = 0; s < 2; s++) {   // child s: row br on side s; delete primes needing the other side
      if (sp >= capst) { capst *= 2; st = realloc(st, capst * sizeof(nd_t)); }
      uint8_t *a = malloc(nP); memcpy(a, smax, nP);
      for (int x = 0; x < needC[1 - s][br]; x++) a[needS[1 - s][br][x]] = 0;
      st[sp++].alive = a;
    }
    free(nd.alive);
  }
  printf("bb2 %s lam %.6f T %.6f status %s M_ub %.6f best_exact %.6f N %d bits %.4f cons %d nodes %lld exact %lld pruned %lld time %.1f\n",
         argv[1], lam, T, status ? "NODECAP" : "DONE", Mub, best, bestN, bestbits, bestcons, nodes, leaves_exact, pruned, (double)(clock() - t0) / CLOCKS_PER_SEC);
  return 0;
}
