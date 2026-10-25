"""Experiment A: Sylow-split solver versus generic k-tree on sub-instances of the real BPSW instance.

Sub-instance: G' = H + E, H = odd Sylow components (l >= 5) of total ~b_H bits (target 0, as for every odd
component of the real instance), E = r components Z/2 (L-side rows with v2(q - 1) = 1, target 1).
Columns: the real dlog columns (comp_T1500.npz from sylow.py) or a control with uniform random columns on the
same components ('rand').

GEN  : Wagner k-tree on H + E: 2^h leaf lists, each = all 2^b subsets of b fresh primes; levels kill about b bits
       each, the root kills the rest; h is the least with (h+1) b >= bits + 2 (same rule for the SPLIT blocks). Success = a subset with image (0, t_E).
SPLIT: k = r + extra disjoint blocks, each a k-tree on H alone (target 0) returning NONEMPTY zero-sum subsets T_i;
       the empty set is the other option of block i. Then GF(2) elimination on E: find c in {0,1}^k with
       sum c_i E(T_i) = t_E. Success = solvable with c != 0.
Cost measures: primes consumed, list elements generated (time proxy), wall time.
Usage: python3 split_vs_gen.py comp.npz b bH r trials seed [rand]
"""
import sys, time, math
import numpy as np

def load(fn):
    z = np.load(fn, allow_pickle=True)
    return z['meta'], z['mat'], z['tgt']

def pick_components(meta, bH, r, rng):
    """H: odd l >= 5 components with l^v <= 2^10 (so a level key fits), drawn at random until ~bH bits.
    E: r components with l = 2, v = 1, side L (target 1)."""
    odd = [i for i, (l, v, q, e, s) in enumerate(meta) if l >= 5 and l ** v <= 1024]
    rng.shuffle(odd)
    H, bits = [], 0.0
    for i in odd:
        lb = math.log2(meta[i][0] ** meta[i][1])
        if bits + lb > bH + 0.5: continue
        H.append(i); bits += lb
        if bits >= bH - 0.5: break
    E2 = [i for i, (l, v, q, e, s) in enumerate(meta) if l == 2 and v == 1 and s == 'L']
    rng.shuffle(E2)
    return H, E2[:r], bits

def levels(mods, b, h):
    """Assign component indices (positions in mods) to h merge levels: levels 1..h-1 about b bits, root the rest."""
    order = sorted(range(len(mods)), key=lambda j: -mods[j])
    lv, cur, acc = [[] for _ in range(h)], 0, 0.0
    for j in order:
        lb = math.log2(mods[j])
        if cur < h - 1 and acc + lb > b + 0.3:
            cur += 1; acc = 0.0
        lv[cur].append(j); acc += lb
    return lv

def enc(img, cols, mods):
    k = np.zeros(img.shape[0], dtype=np.int64)
    for j in cols: k = k * int(mods[j]) + img[:, j]
    return k

def leaf(cols_block, mods):
    """All 2^b subset images of a block (b x R array of residues). Returns (imgs, masks)."""
    imgs = np.zeros((1, len(mods)), dtype=np.int64)
    for c in cols_block: imgs = np.vstack([imgs, (imgs + c) % mods])
    masks = np.arange(imgs.shape[0], dtype=np.int64)
    return imgs, masks

def merge(A, B, cols, mods, cap, rng, stats):
    """A, B: (imgs, prov). Keep pairs with imgs_a + imgs_b == 0 on cols."""
    ia, pa = A; ib, pb = B
    ka = enc(ia, cols, mods); kb = enc((-ib) % mods, cols, mods)
    ob = np.argsort(kb, kind='stable'); kbs = kb[ob]
    lo = np.searchsorted(kbs, ka, 'left'); hi = np.searchsorted(kbs, ka, 'right')
    cnt = hi - lo; tot = int(cnt.sum())
    stats['ops'] += ia.shape[0] + ib.shape[0] + tot
    if tot == 0: return np.zeros((0, len(mods)), dtype=np.int64), (np.zeros(0, np.int64), np.zeros(0, np.int64))
    ai = np.repeat(np.arange(len(ka)), cnt)
    off = np.arange(tot) - np.repeat(np.cumsum(cnt) - cnt, cnt)
    bi = ob[np.repeat(lo, cnt) + off]
    if tot > cap:
        sel = rng.choice(tot, cap, replace=False); ai, bi = ai[sel], bi[sel]
    return (ia[ai] + ib[bi]) % mods, (ai, bi)

def ktree(colmat, mods, tgt, b, h, cap, rng, stats, want_nonempty=False):
    """colmat: (2^h * b) x R residues of the primes used. Returns list of solution subsets (lists of row indices)."""
    R = len(mods); nl = 2 ** h
    lv = levels(mods, b, h)
    lists = []
    for i in range(nl):
        imgs, masks = leaf(colmat[i * b:(i + 1) * b], mods)
        if i == 0: imgs = (imgs - tgt) % mods          # target subtracted once
        lists.append((imgs, [('leaf', i, masks)]))
    tree = []
    for lvl in range(h):
        new = []
        for j in range(0, len(lists), 2):
            (ia, ta), (ib, tb) = lists[j], lists[j + 1]
            imgs, (ai, bi) = merge((ia, None), (ib, None), lv[lvl], mods, cap, rng, stats)
            new.append((imgs, [('node', ta, tb, ai, bi)]))
        lists = new
    imgs, prov = lists[0]
    assert imgs.shape[0] == 0 or not np.any(imgs)       # root: all components killed
    def subset(p, k):
        if p[0] == 'leaf':
            _, i, masks = p; m = int(masks[k]); return [i * b + t for t in range(b) if m >> t & 1]
        _, ta, tb, ai, bi = p
        return subset(ta[0], int(ai[k])) + subset(tb[0], int(bi[k]))
    sols = []
    for k in range(imgs.shape[0]):
        s = subset(prov[0], k)
        if want_nonempty and not s: continue
        sols.append(s)
        if len(sols) >= 8: break
    return sols

def gf2_solve(cols, t):
    """cols: list of int bitmasks (E-images), t bitmask. Return a 0/1 combination (list of indices) or None."""
    basis = {}                                    # pivot bit -> (vector, combination bitmask)
    for i, v in enumerate(cols):
        c = 1 << i
        while v:
            hb = v.bit_length() - 1
            if hb in basis: bv, bc = basis[hb]; v ^= bv; c ^= bc
            else: basis[hb] = (v, c); break
    c = 0
    while t:
        hb = t.bit_length() - 1
        if hb not in basis: return None
        bv, bc = basis[hb]; t ^= bv; c ^= bc
    return [i for i in range(len(cols)) if c >> i & 1]

def main():
    fn, b, bH, r, trials, seed = sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
    rand = len(sys.argv) > 7 and sys.argv[7] == 'rand'
    meta, mat, tgt = load(fn)
    N = mat.shape[1]
    rng = np.random.default_rng(seed)
    cap = 4 * 2 ** b
    res = {'gen': [], 'split': []}
    for tr in range(trials):
        H, E, hb = pick_components(meta, bH, r, rng)
        comps = H + E
        mods = np.array([meta[i][0] ** meta[i][1] for i in comps], dtype=np.int64)
        cols = mat[comps].T.copy()                     # N x R
        if rand: cols = rng.integers(0, mods, size=cols.shape)
        t = tgt[comps].copy()
        Btot = hb + r
        # ---- GEN: smallest h with (h+1) b >= Btot + 1
        h = max(1, math.ceil((Btot + 2) / b) - 1)
        need = 2 ** h * b
        st = {'ops': 0}; t0 = time.time()
        if need <= N:
            perm = rng.permutation(N)[:need]
            sols = ktree(cols[perm], mods, t, b, h, cap, rng, st)
            ok = bool(sols)
            if ok:
                s = [int(perm[i]) for i in sols[0]]
                assert np.all((cols[s].sum(0) - t) % mods == 0)
        else:
            ok = None
        res['gen'].append((ok, need, st['ops'], time.time() - t0, h))
        # ---- SPLIT: blocks with k-tree on H only, then GF(2) on E
        nH = len(H); modsH = mods[:nH]; tH = t[:nH]
        assert not np.any(tH)
        hH = max(1, math.ceil((hb + 2) / b) - 1)
        per = 2 ** hH * b
        st2 = {'ops': 0}; t0 = time.time()
        perm = rng.permutation(N); used = 0; Ts = []
        kmax = min(r + 6, N // per)
        for blk in range(kmax):
            idx = perm[used:used + per]; used += per
            sols = ktree(cols[idx][:, :nH], modsH, tH, b, hH, cap, rng, st2, want_nonempty=True)
            if sols: Ts.append([int(idx[i]) for i in sols[0]])
            if len(Ts) >= r + 4: break
        Ev = []
        for T in Ts:
            v = 0
            for j in range(r):
                if int(cols[T, nH + j].sum()) & 1: v |= 1 << j
            Ev.append(v)
        tv = sum(1 << j for j in range(r) if t[nH + j] & 1)
        c = gf2_solve(Ev, tv)
        ok2 = c is not None and len(c) > 0
        if ok2:
            S = sorted(x for i in c for x in Ts[i]); assert len(S) == len(set(S))
            assert np.all((cols[S].sum(0) - t) % mods == 0)
        res['split'].append((ok2, used, st2['ops'], time.time() - t0, len(Ts)))
        print(f"trial {tr}: bits H {hb:.1f} + E {r} = {Btot:.1f}; GEN h={h} primes {need} ok={ok} ops={st['ops']} "
              f"| SPLIT blocks {len(Ts)} (hH={hH}, {per}/block) primes {used} ok={ok2} ops={st2['ops']}", flush=True)
    for k in ('gen', 'split'):
        rs = [x for x in res[k] if x[0] is not None]
        if not rs: print(f"{k}: not runnable (needs more primes than the pool)"); continue
        sr = np.mean([x[0] for x in rs]); pr = np.mean([x[1] for x in rs]); op = np.mean([x[2] for x in rs])
        tm = np.mean([x[3] for x in rs])
        print(f"SUMMARY {'rand' if rand else 'real'} b={b} bH={bH} r={r} {k}: success {sr:.2f} over {len(rs)}, "
              f"primes {pr:.0f}, list-ops {op:.3g} (log2 {math.log2(op):.1f}), time {tm:.2f}s")

if __name__ == "__main__":
    main()
