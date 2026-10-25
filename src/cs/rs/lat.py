"""Experiment B: lattice (LLL / BKZ / enumeration) on the relation lattice of sub-instances, real columns vs controls.

Sub-instance: whole rows (q, e) of the real instance, drawn at random until ~Bp bits; cyclic component Z/phi(q^e)
per row; column of p = log_g(p mod q^e); target 0 (F row) or phi/2 (L row). This is the instance restricted to
those rows, so a relation of the sub-instance is a congruence prod p^(x+) == prod p^(x-) mod prod q^e.
Columns: 'real' (pool primes), controls 'unif' (uniform residues), 'supp' (real entry where p == +-1 mod q,
uniform elsewhere: keeps the forced support pattern), 'hgt' (random integers of the same bit sizes as the chosen
primes, coprime to the rows: keeps the height, destroys primality and the support pattern).
Prime choice: 'rnd' (random pool primes) or 'low' (the d smallest pool primes: residues are unreduced, so the
height constraint is as strong as the pool allows).
Lattice (CJLOSS-type embedding, dimension d + 1 + R): rows (2 e_i, 0, K c_i), (1..1, 1, K t), (0, 0, K phi_j e_j).
A solution x in {0,1}^d gives (2x - 1, -1, 0), squared norm d + 1. After LLL the first d + 1 rows span the K-free
sublattice; BKZ-beta on it; success if a reduced basis vector is (+-1, ..., +-1, +-1). Optionally enumerate every
vector of squared norm <= d + 1 and count solutions (exact for small d).
Also reports the shortest relation (v[d] = 0 vectors are 2x with x in the relation lattice) against the Gaussian
heuristic 2 * GH(Lambda), GH = det^(1/d) sqrt(d / (2 pi e)), det = |<columns>| (= |G'| when they generate).
Usage: python3 lat.py instance.json Bp dspec(+c: d = bits + c; f: d = f * bits) trials seed mode{rnd,low} ctl{real,unif,supp,hgt} beta [enum]
"""
import sys, json, math, time
import numpy as np
from fpylll import IntegerMatrix, LLL, BKZ, GSO, Enumeration, EnumerationError

def factor(m):
    f, q = {}, 2
    while q * q <= m:
        while m % q == 0: f[q] = f.get(q, 0) + 1; m //= q
        q += 1
    if m > 1: f[m] = f.get(m, 0) + 1
    return f

def prim_root(q, e):
    m, ph = q ** e, (q - 1) * q ** (e - 1)
    for g in range(2, m):
        if g % q and all(pow(g, ph // r, m) != 1 for r in factor(ph)): return g

def table(q, e):
    m, ph = q ** e, (q - 1) * q ** (e - 1); g = prim_root(q, e)
    lg = np.full(m, -1, dtype=np.int64); x = 1
    for k in range(ph): lg[x] = k; x = x * g % m
    return lg

def main():
    a = sys.argv
    inst = json.load(open(a[1])); Bp = float(a[2]); dspec = a[3]; trials = int(a[4]); seed = int(a[5])
    mode, ctl, beta = a[6], a[7], int(a[8]); do_enum = len(a) > 9
    rng = np.random.default_rng(seed)
    ps_all = sorted(int(t[0]) for t in inst["primes"])
    rows_all = [(int(q), e, 'F') for q, e in inst["L1_factors"].items()] + [(int(q), e, 'L') for q, e in inst["L2_factors"].items()]
    tabs = {}
    stat = []
    for tr in range(trials):
        order = rng.permutation(len(rows_all)); rows, bits = [], 0.0
        for k in order:
            q, e, s = rows_all[k]; ph = (q - 1) * q ** (e - 1); lb = math.log2(ph)
            if q ** e > 200000 or bits + lb > Bp + 1: continue
            rows.append((q, e, s)); bits += lb
            if bits >= Bp - 1: break
        d = int(round(bits)) + int(dspec[1:]) if dspec[0] == '+' else int(round(float(dspec) * bits))
        ps = ps_all[:d] if mode == 'low' else [ps_all[i] for i in sorted(rng.choice(len(ps_all), d, replace=False))]
        if ctl == 'hgt':
            M = math.prod(q for q, _, _ in rows); new = []
            for p in ps:
                nb = p.bit_length()
                while True:
                    x = int(rng.integers(2 ** (nb - 1), 2 ** nb - 1)) | 1
                    if math.gcd(x, M) == 1 and x not in new: break
                new.append(x)
            ps = new
        C, mods, tg = [], [], []
        for q, e, s in rows:
            if (q, e) not in tabs: tabs[(q, e)] = table(q, e)
            lg = tabs[(q, e)]; m = q ** e; ph = (q - 1) * q ** (e - 1)
            col = np.array([lg[p % m] for p in ps], dtype=np.int64); assert np.all(col >= 0)
            if ctl in ('unif', 'supp'):
                rnd = rng.integers(0, ph, size=d)
                if ctl == 'supp':
                    frc = np.array([(p % q) in (1, q - 1) for p in ps]); col = np.where(frc, col, rnd)
                else: col = rnd
            C.append(col); mods.append(ph); tg.append(0 if s == 'F' else ph // 2)
        R = len(rows); K = 2 ** 40
        D = d + 1 + R
        B = [[0] * D for _ in range(D)]
        for i in range(d):
            B[i][i] = 2
            for j in range(R): B[i][d + 1 + j] = K * int(C[j][i])
        for i in range(d): B[d][i] = 1
        B[d][d] = 1
        for j in range(R): B[d][d + 1 + j] = K * tg[j]
        for j in range(R): B[d + 1 + j][d + 1 + j] = K * mods[j]
        A = IntegerMatrix.from_matrix(B)
        t0 = time.time()
        LLL.reduction(A)
        sub = [[A[i, j] for j in range(d + 1)] for i in range(D) if all(A[i, j] == 0 for j in range(d + 1, D))]
        assert len(sub) == d + 1, len(sub)
        S = IntegerMatrix.from_matrix(sub)
        if beta > 2: BKZ.reduction(S, BKZ.Param(block_size=min(beta, d + 1), max_loops=6,
                                                 flags=BKZ.AUTO_ABORT | BKZ.MAX_LOOPS))
        tb = time.time() - t0
        def is_sol(v): return abs(v[d]) == 1 and all(abs(x) == 1 for x in v[:d])
        found = 0; shortest_rel = None
        for i in range(d + 1):
            v = [S[i, j] for j in range(d + 1)]
            if is_sol(v): found += 1
            if v[d] == 0:
                n2 = sum(x * x for x in v)
                if shortest_rel is None or n2 < shortest_rel: shortest_rel = n2
        # Gaussian heuristic for the relation lattice: det = |<columns>| <= prod phi; we use log2 prod phi = bits
        gh = 2 ** (bits / d) * math.sqrt(d / (2 * math.pi * math.e))
        ne, npts, te = None, None, 0.0
        if do_enum:
            t1 = time.time()
            Mg = GSO.Mat(S); Mg.update_gso()
            try:
                sols = Enumeration(Mg, nr_solutions=20000).enumerate(0, d + 1, d + 1 + 0.01, 0)
            except EnumerationError:
                sols = []
            npts = len(sols); ne = 0
            for _, co in sols:
                v = [sum(int(round(co[i])) * S[i, j] for i in range(d + 1)) for j in range(d + 1)]
                if is_sol(v): ne += 1
            te = time.time() - t1
        stat.append((found, shortest_rel, gh, ne, npts, tb, te, d, bits))
        print(f"trial {tr}: rows {R} bits {bits:.1f} d {d} | BKZ-{beta} basis solutions {found}, shortest relation "
              f"|x| = {math.sqrt(shortest_rel)/2 if shortest_rel else float('nan'):.3f} vs GH {gh:.3f}, {tb:.1f}s"
              + (f" | enum: {ne} solutions among {npts} points, {te:.1f}s" if do_enum else ""), flush=True)
    f = [s[0] > 0 for s in stat]
    rel = [math.sqrt(s[1]) / 2 / s[2] for s in stat if s[1]]
    line = (f"SUMMARY mode={mode} ctl={ctl} Bp={Bp} d={dspec} beta={beta}: BKZ success {np.mean(f):.2f} ({sum(f)}/{len(f)}), "
            f"shortest-relation/GH mean {np.mean(rel):.4f} sd {np.std(rel):.4f}, time {np.mean([s[5] for s in stat]):.1f}s")
    if do_enum:
        line += (f"; enum success {np.mean([s[3] > 0 for s in stat]):.2f}, mean #solutions {np.mean([s[3] for s in stat]):.2f}, "
                 f"mean #points {np.mean([s[4] for s in stat]):.1f}, enum time {np.mean([s[6] for s in stat]):.1f}s")
    print(line)

if __name__ == "__main__":
    main()
