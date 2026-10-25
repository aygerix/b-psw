"""Own parser/builder for bbc instance files (adversarial FX check; independent of e1/lagr_bb.py and mkinst.py)."""
import math
def read_inst(fn):
    with open(fn) as f: L = f.read().split('\n')
    nP, nI, nR = map(int, L[0].split())
    items = []
    for l in L[1:1 + nI]:
        a = l.split(); items.append((int(a[0]), int(a[1]), int(a[2]), int(a[3]), float(a[4])))
    prim = []
    for l in L[1 + nI:1 + nI + nP]:
        a = l.split(); p, k = int(a[0]), int(a[1]); its = [int(x) for x in a[2:]]
        assert len(its) == k; prim.append((p, its))
    assert len(prim) == nP
    return nP, nI, nR, items, prim
def item_q(items):
    """recover q of every item from its cost: level 1 -> log2(q-1), level >= 2 -> log2(q); check primality/consistency"""
    qs = []
    for (r, s, e, par, c) in items:
        q = round(2 ** c) + 1 if e == 1 else round(2 ** c)
        assert abs(c - (math.log2(q - 1) if e == 1 else math.log2(q))) < 1e-9, (r, s, e, c)
        qs.append(q)
    return qs
def smallfac(x, T0):
    f = {}
    q = 3
    while x > 1 and q <= T0:
        if x % q == 0:
            e = 0
            while x % q == 0: x //= q; e += 1
            f[q] = e
        q += 2
    if x != 1: return None
    return f
def req_from_record(o1, o2, T0):
    f1 = smallfac(o1, T0); f2 = smallfac(o2, T0)
    if f1 is None or f2 is None: return None
    return sorted([(q, 0, a) for q, a in f1.items()] + [(q, 1, a) for q, a in f2.items()])
def inst_reqs(items, prim):
    qs = item_q(items)
    out = {}
    for p, its in prim:
        out[p] = sorted((qs[j], items[j][1], items[j][2]) for j in its)
    return out, qs
def check_items(items, qs):
    """item chain sanity: parent of (q,s,e) is (q,s,e-1); same row index <-> same q; cost formula"""
    key = {}
    for j, (r, s, e, par, c) in enumerate(items): key[(qs[j], s, e)] = j
    bad = 0; rowq = {}
    for j, (r, s, e, par, c) in enumerate(items):
        if rowq.setdefault(r, qs[j]) != qs[j]: bad += 1
        if e == 1 and par != -1: bad += 1
        if e > 1 and par != key.get((qs[j], s, e - 1), -5): bad += 1
    return bad
