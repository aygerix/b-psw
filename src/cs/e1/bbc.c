// E1 Lagrangian branch-and-bound over row sides, C port of lagr_bb.py (same model, same certificate).
// For lambda > 0, M(lambda) = max over configurations C of N(C) - lambda * bits(C).
// Node bound: maximum-weight closure (min cut, Dinic) with every FREE row offered on both sides (relaxation).
// Node closed: (a) bound <= incumbent (pruned) or (b) the minimal optimal closure uses no free row on both sides
// (feasible: exact).  Primes outside the MAXIMAL optimal closure are deleted for the whole subtree (monotonicity of
// maximal min-cut source sets when source capacities decrease).  Integer capacities: prime weight K, item cost
// floor(lambda c K) (rounded down => every bound is >= the true relaxed value).  Incumbent values are exact doubles.
// Instance file (from mkinst.py):
//   line 1: n_primes n_items n_rows
//   n_items lines: row_id side level parent_item(-1 if level 1) cost(double)
//   n_primes lines: p k item_1 ... item_k        (item_j = top item (row, side, a) of each requirement)
// Usage: ./bbc instfile lam1[,lam2,...] [nodecap]
// Output per lambda: "lam status M_ub best_value N bits nodes exact pruned seconds" and the hull-point primes.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

#define K (1LL << 20)
#define INFC (1LL << 60)
static int nP, nI, nR; static int SB = 0; static long long sbcl = 0;
static int *irow, *iside, *ilev, *ipar; static double *icost;
static uint64_t *pval; static int *pk, **pit;
static int *needF_cnt, *needL_cnt, **needF, **needL;    // per row: primes needing it on F / L

// ---- Dinic ----
static int NV, NE, *head, *nxt, *to; static int64_t *cap; static int *lvl, *itp, *q;
static void ginit(int nv, int ne) { NV = nv; NE = 0; for (int i = 0; i < nv; i++) head[i] = -1; (void)ne; }
static void addE(int u, int v, int64_t c) {
  to[NE] = v; cap[NE] = c; nxt[NE] = head[u]; head[u] = NE++;
  to[NE] = u; cap[NE] = 0; nxt[NE] = head[v]; head[v] = NE++;
}
static int bfs(int s, int t) {
  for (int i = 0; i < NV; i++) lvl[i] = -1;
  int qh = 0, qt = 0; q[qt++] = s; lvl[s] = 0;
  while (qh < qt) { int u = q[qh++]; for (int e = head[u]; e >= 0; e = nxt[e]) if (cap[e] > 0 && lvl[to[e]] < 0) { lvl[to[e]] = lvl[u] + 1; q[qt++] = to[e]; } }
  return lvl[t] >= 0;
}
static int64_t dfs(int u, int t, int64_t f) {
  if (u == t) return f;
  for (int *e = &itp[u]; *e >= 0; *e = nxt[*e]) {
    int v = to[*e];
    if (cap[*e] > 0 && lvl[v] == lvl[u] + 1) {
      int64_t d = dfs(v, t, f < cap[*e] ? f : cap[*e]);
      if (d > 0) { cap[*e] -= d; cap[*e ^ 1] += d; return d; }
    }
  }
  return 0;
}
static int64_t maxflow(int s, int t) {
  int64_t fl = 0;
  while (bfs(s, t)) { for (int i = 0; i < NV; i++) itp[i] = head[i]; int64_t f; while ((f = dfs(s, t, INFC)) > 0) fl += f; }
  return fl;
}

// ---- B&B ----
typedef struct { int8_t *fixed; uint8_t *alive; } node_t;
static int *amap, *imap, *ilist; static uint8_t *selmin, *selmax, *imin, *imax, *reachS, *reachT;
#define FB 601
static int frontN[FB];                 // frontN[b] = max N over feasible configurations seen with bits <= b (integer b)
static uint8_t *rsel;
static void record(int N, double bits) { int b0 = (int)ceil(bits - 1e-9); if (b0 < 0) b0 = 0; for (int b = b0; b < FB; b++) if (N > frontN[b]) frontN[b] = N; else if (b > b0 + 0 && frontN[b] >= N) break; }

// returns bound; fills selmin/selmax (primes), imin/imax (items) for the alive set
static double closure(const uint8_t *alive, double lam, int *nsel) {
  int na = 0; for (int i = 0; i < nP; i++) if (alive[i]) amap[i] = na++; else amap[i] = -1;
  // items used by alive primes (+ chains)
  int ni = 0; for (int j = 0; j < nI; j++) imap[j] = -1;
  for (int i = 0; i < nP; i++) if (alive[i]) for (int r = 0; r < pk[i]; r++) {
    int j = pit[i][r]; while (j >= 0 && imap[j] < 0) { imap[j] = -2; j = ipar[j]; }
  }
  for (int j = 0; j < nI; j++) if (imap[j] == -2) { imap[j] = ni; ilist[ni++] = j; }
  int S = 0, T = 1, base = 2 + na;
  ginit(base + ni, 0);
  for (int i = 0; i < nP; i++) if (alive[i]) {
    addE(S, 2 + amap[i], K);
    for (int r = 0; r < pk[i]; r++) addE(2 + amap[i], base + imap[pit[i][r]], INFC);
  }
  for (int x = 0; x < ni; x++) {
    int j = ilist[x];
    if (ipar[j] >= 0) addE(base + x, base + imap[ipar[j]], INFC);
    int64_t c = (int64_t)floor(lam * icost[j] * (double)K * (1 - 1e-12));
    addE(base + x, T, c);
  }
  int64_t fl = maxflow(S, T);
  double ub = ((double)na * K - (double)fl) / (double)K;
  // minimal closure: reachable from S in residual
  for (int v = 0; v < NV; v++) reachS[v] = reachT[v] = 0;
  int qh = 0, qt = 0; q[qt++] = S; reachS[S] = 1;
  while (qh < qt) { int u = q[qh++]; for (int e = head[u]; e >= 0; e = nxt[e]) if (cap[e] > 0 && !reachS[to[e]]) { reachS[to[e]] = 1; q[qt++] = to[e]; } }
  // nodes that can reach T: reverse BFS; edge e: u -> v (e even or odd), residual cap[e]; partner e^1 is v -> u
  qh = qt = 0; q[qt++] = T; reachT[T] = 1;
  while (qh < qt) { int v = q[qh++]; for (int e = head[v]; e >= 0; e = nxt[e]) { int u = to[e]; if (cap[e ^ 1] > 0 && !reachT[u]) { reachT[u] = 1; q[qt++] = u; } } }
  int ns = 0;
  for (int i = 0; i < nP; i++) {
    selmin[i] = selmax[i] = 0;
    if (alive[i]) { int v = 2 + amap[i]; selmin[i] = reachS[v]; selmax[i] = !reachT[v]; ns += selmin[i]; }
  }
  for (int j = 0; j < nI; j++) {
    imin[j] = imax[j] = 0;
    if (imap[j] >= 0) { int v = base + imap[j]; imin[j] = reachS[v]; imax[j] = !reachT[v]; }
  }
  *nsel = ns;
  return ub;
}

static double exact_value(const uint8_t *sel, double lam, int *N, double *bits) {
  static uint8_t *u = 0; if (!u) u = malloc(nI);
  memset(u, 0, nI); int n = 0;
  for (int i = 0; i < nP; i++) if (sel[i]) { n++; for (int r = 0; r < pk[i]; r++) { int j = pit[i][r]; while (j >= 0 && !u[j]) { u[j] = 1; j = ipar[j]; } } }
  double b = 0; for (int j = 0; j < nI; j++) if (u[j]) b += icost[j];
  *N = n; *bits = b; return n - lam * b;
}

int main(int argc, char **argv) {
  if (argc < 3) { fprintf(stderr, "usage: bbc inst lam1[,lam2] [nodecap]\n"); return 1; }
  FILE *f = fopen(argv[1], "r"); if (!f) { perror("inst"); return 1; }
  if (fscanf(f, "%d %d %d", &nP, &nI, &nR) != 3) return 1;
  irow = malloc(nI * sizeof(int)); iside = malloc(nI * sizeof(int)); ilev = malloc(nI * sizeof(int));
  ipar = malloc(nI * sizeof(int)); icost = malloc(nI * sizeof(double));
  for (int j = 0; j < nI; j++) if (fscanf(f, "%d %d %d %d %lf", &irow[j], &iside[j], &ilev[j], &ipar[j], &icost[j]) != 5) return 1;
  pval = malloc(nP * sizeof(uint64_t)); pk = malloc(nP * sizeof(int)); pit = malloc(nP * sizeof(int *));
  int maxE = 0;
  for (int i = 0; i < nP; i++) {
    unsigned long long pv; if (fscanf(f, "%llu %d", &pv, &pk[i]) != 2) return 1; pval[i] = pv;
    pit[i] = malloc((pk[i] + 1) * sizeof(int)); maxE += pk[i] + 1;
    for (int r = 0; r < pk[i]; r++) if (fscanf(f, "%d", &pit[i][r]) != 1) return 1;
  }
  fclose(f);
  needF_cnt = calloc(nR, sizeof(int)); needL_cnt = calloc(nR, sizeof(int));
  needF = malloc(nR * sizeof(int *)); needL = malloc(nR * sizeof(int *));
  for (int i = 0; i < nP; i++) for (int r = 0; r < pk[i]; r++) { int j = pit[i][r]; if (iside[j] == 0) needF_cnt[irow[j]]++; else needL_cnt[irow[j]]++; }
  for (int r = 0; r < nR; r++) { needF[r] = malloc((needF_cnt[r] + 1) * sizeof(int)); needL[r] = malloc((needL_cnt[r] + 1) * sizeof(int)); needF_cnt[r] = needL_cnt[r] = 0; }
  for (int i = 0; i < nP; i++) for (int r = 0; r < pk[i]; r++) { int j = pit[i][r]; if (iside[j] == 0) needF[irow[j]][needF_cnt[irow[j]]++] = i; else needL[irow[j]][needL_cnt[irow[j]]++] = i; }
  int nvmax = 2 + nP + nI, nemax = 2 * (maxE + 2 * nI + 4);
  head = malloc(nvmax * sizeof(int)); nxt = malloc(nemax * sizeof(int)); to = malloc(nemax * sizeof(int));
  cap = malloc(nemax * sizeof(int64_t)); lvl = malloc(nvmax * sizeof(int)); itp = malloc(nvmax * sizeof(int)); q = malloc(nvmax * sizeof(int));
  amap = malloc(nP * sizeof(int)); imap = malloc(nI * sizeof(int)); ilist = malloc(nI * sizeof(int));
  rsel = malloc(nP); selmin = malloc(nP); selmax = malloc(nP); imin = malloc(nI); imax = malloc(nI); reachS = malloc(nvmax); reachT = malloc(nvmax);
  long long nodecap = argc > 3 ? atoll(argv[3]) : 100000000LL;
  if (getenv("SB")) SB = atoi(getenv("SB"));
  char *lp = argv[2];
  while (lp && *lp) {
    double lam = strtod(lp, &lp); if (*lp == ',') lp++;
    clock_t t0 = clock();
    // DFS stack
    int cap_st = 4 * (nR + 2);
    node_t *st = malloc(cap_st * sizeof(node_t)); int sp = 0;
    st[sp].fixed = malloc(nR); memset(st[sp].fixed, -1, nR); st[sp].alive = malloc(nP); memset(st[sp].alive, 1, nP); sp++;
    memset(frontN, 0, sizeof frontN);
    double bestv = 0.0, bestbits = 0.0, Mub = -1e300; int bestN = 0; uint8_t *bestsel = calloc(nP, 1);
    long long nodes = 0, nexact = 0, npruned = 0; int status = 0;
    uint8_t *usedF = malloc(nR), *usedL = malloc(nR);
    while (sp > 0) {
      node_t nd = st[--sp];
      if (++nodes > nodecap) { status = 1; free(nd.fixed); free(nd.alive); break; }
      int nsel; double ub = closure(nd.alive, lam, &nsel);
      if (ub <= bestv + 1e-9) { npruned++; if (ub > Mub) Mub = ub; free(nd.fixed); free(nd.alive); continue; }
      // conflicts in minimal closure
      memset(usedF, 0, nR); memset(usedL, 0, nR);
      for (int j = 0; j < nI; j++) if (imin[j] && ilev[j] == 1 && nd.fixed[irow[j]] < 0) { if (iside[j] == 0) usedF[irow[j]] = 1; else usedL[irow[j]] = 1; }
      int conf = 0; for (int r = 0; r < nR; r++) if (usedF[r] && usedL[r]) { conf = 1; break; }
      if (!conf) {
        int N; double bits; double v = exact_value(selmin, lam, &N, &bits);
        nexact++; if (ub > Mub) Mub = ub; record(N, bits);
        if (v > bestv) { bestv = v; bestN = N; bestbits = bits; memcpy(bestsel, selmin, nP); }
        free(nd.fixed); free(nd.alive); continue;
      }
      // branch on a row used on both sides by the maximal closure; weight = selected primes per side
      memset(usedF, 0, nR); memset(usedL, 0, nR);
      for (int j = 0; j < nI; j++) if (imax[j] && ilev[j] == 1 && nd.fixed[irow[j]] < 0) { if (iside[j] == 0) usedF[irow[j]] = 1; else usedL[irow[j]] = 1; }
      int br = -1, bw = -1, bF = 0, bL = 0;
      for (int r = 0; r < nR; r++) if (usedF[r] && usedL[r]) {
        int wF = 0, wL = 0;
        for (int x = 0; x < needF_cnt[r]; x++) wF += selmax[needF[r][x]];
        for (int x = 0; x < needL_cnt[r]; x++) wL += selmax[needL[r][x]];
        int w = wF < wL ? wF : wL;
        if (w >= bw) { bw = w; br = r; bF = wF; bL = wL; }
      }
      if (br < 0) {   // minimal closure conflicted but maximal not (cannot happen: min subset of max); be safe
        fprintf(stderr, "internal: no branch row\n"); return 2;
      }
      {   // primal rounding: every free conflicted row to its majority side (by selmax), drop the other side
        memcpy(rsel, selmax, nP);
        for (int r = 0; r < nR; r++) if (usedF[r] && usedL[r]) {
          int wF = 0, wL = 0;
          for (int x = 0; x < needF_cnt[r]; x++) wF += rsel[needF[r][x]];
          for (int x = 0; x < needL_cnt[r]; x++) wL += rsel[needL[r][x]];
          int *lst = wF >= wL ? needL[r] : needF[r]; int cnt = wF >= wL ? needL_cnt[r] : needF_cnt[r];
          for (int x = 0; x < cnt; x++) rsel[lst[x]] = 0;
        }
        int N2; double b2; double v2 = exact_value(rsel, lam, &N2, &b2); record(N2, b2);
        if (v2 > bestv) { bestv = v2; bestN = N2; bestbits = b2; memcpy(bestsel, rsel, nP); }
      }
      if (SB > 0) {   // strong branching: among the SB heaviest conflicted rows pick the one minimising max(child bounds);
                      // max over the two children of the chosen row is a valid bound for this node (prune if <= incumbent)
        static uint8_t *savemax = 0, *ca = 0; if (!savemax) { savemax = malloc(nP); ca = malloc(nP); }
        memcpy(savemax, selmax, nP);
        int cand[64], nc = 0, cw[64];
        for (int r = 0; r < nR; r++) if (usedF[r] && usedL[r]) {
          int wF = 0, wL = 0;
          for (int x = 0; x < needF_cnt[r]; x++) wF += savemax[needF[r][x]];
          for (int x = 0; x < needL_cnt[r]; x++) wL += savemax[needL[r][x]];
          int w = wF < wL ? wF : wL;
          if (nc < 64) { cand[nc] = r; cw[nc] = w; nc++; }
        }
        for (int a = 0; a < nc; a++) for (int b = a + 1; b < nc; b++) if (cw[b] > cw[a]) { int t = cw[a]; cw[a] = cw[b]; cw[b] = t; t = cand[a]; cand[a] = cand[b]; cand[b] = t; }
        double bestscore = 1e300; int bestr = -1;
        for (int a = 0; a < nc && a < SB; a++) {
          int r = cand[a]; double sc = -1e300;
          for (int s2 = 0; s2 < 2; s2++) {
            memcpy(ca, savemax, nP);
            int *lst = s2 == 0 ? needL[r] : needF[r]; int cnt = s2 == 0 ? needL_cnt[r] : needF_cnt[r];
            for (int x = 0; x < cnt; x++) ca[lst[x]] = 0;
            int ns2; double u2 = closure(ca, lam, &ns2); sbcl++;
            if (u2 > sc) sc = u2;
          }
          if (sc < bestscore) { bestscore = sc; bestr = r; }
        }
        memcpy(selmax, savemax, nP);
        if (bestscore <= bestv + 1e-9) { npruned++; if (bestscore > Mub) Mub = bestscore; free(nd.fixed); free(nd.alive); continue; }
        br = bestr; bF = bL = 0;
        for (int x = 0; x < needF_cnt[br]; x++) bF += selmax[needF[br][x]];
        for (int x = 0; x < needL_cnt[br]; x++) bL += selmax[needL[br][x]];
      }
      int first = bF >= bL ? 0 : 1;
      for (int k = 0; k < 2; k++) {
        int s = k == 0 ? 1 - first : first;
        if (sp >= cap_st) { cap_st *= 2; st = realloc(st, cap_st * sizeof(node_t)); }
        node_t ch; ch.fixed = malloc(nR); memcpy(ch.fixed, nd.fixed, nR); ch.fixed[br] = s;
        ch.alive = malloc(nP); memcpy(ch.alive, selmax, nP);
        int *lst = s == 0 ? needL[br] : needF[br]; int cnt = s == 0 ? needL_cnt[br] : needF_cnt[br];
        for (int x = 0; x < cnt; x++) ch.alive[lst[x]] = 0;
        st[sp++] = ch;
      }
      free(nd.fixed); free(nd.alive);
    }
    while (sp > 0) { sp--; free(st[sp].fixed); free(st[sp].alive); }
    free(st);
    if (bestv > Mub) Mub = bestv;
    double secs = (double)(clock() - t0) / CLOCKS_PER_SEC;
    printf("lam %.6f %s M_ub %.6f best %.6f N %d bits %.3f nodes %lld exact %lld pruned %lld time %.2f sb %d sbclosures %lld\n", lam,
           status ? "NODECAP" : "DONE", status ? NAN : Mub, bestv, bestN, bestbits, nodes, nexact, npruned, secs, SB, sbcl);
    sbcl = 0;
    printf("  hull primes:");
    for (int i = 0; i < nP; i++) if (bestsel[i]) printf(" %llu", (unsigned long long)pval[i]);
    printf("\n  primal frontier (b:N):");
    for (int b = 20; b < FB; b += 20) printf(" %d:%d", b, frontN[b]);
    printf(" 214:%d 275:%d\n", frontN[214], frontN[275]); fflush(stdout);
    free(bestsel); free(usedF); free(usedL);
    if (status) break;
  }
  return 0;
}
