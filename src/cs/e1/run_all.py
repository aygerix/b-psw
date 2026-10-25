"""E1 driver: Lagrangian B&B certificate for every (D, 2-adic class) of the pool, lambdas in decreasing order.
For each instance: M_ub(lambda) (certified), the hull point (bits, N) of the incumbent, nodes, time.
Envelope UB(b) = min over finished lambdas of lambda*b + M_ub(lambda); kill check UB(b) < (b+1)/3 for b in BGRID.
A lambda that hits NODECAP stops the descent for that instance (smaller lambdas are harder).
Usage: python3 run_all.py OUT.tsv lam1,lam2,... D1[,D2...] recfile1 [recfile2 ...]   [env NODECAP, MINPOOL, CLASSES]
"""
import sys, os, time
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lagr_bb as L

BGRID = [60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 275, 300]

def main():
    out = sys.argv[1]; lams = [float(x) for x in sys.argv[2].split(',')]
    Ds = [int(x) for x in sys.argv[3].split(',')]; files = sys.argv[4:]
    minpool = int(os.environ.get('MINPOOL', 1))
    only = os.environ.get('CLASSES')
    only = None if not only else {tuple(map(int, c.split(':'))) for c in only.split(',')}
    recs = {}
    for fn in files:
        for l in open(fn):
            a = l.split()
            if len(a) == 7: recs.setdefault(int(a[0]), {})[int(a[1])] = (int(a[3]), int(a[4]))
    fo = open(out, 'a')
    for D in Ds:
        cls = Counter(recs.get(D, {}).values())
        for c, cnt in sorted(cls.items(), key=lambda t: -t[1]):
            if cnt < minpool or (only and c not in only): continue
            pool, bad = L.load(D, c[0], c[1], files)
            inst = L.Inst(pool)
            res = []
            for lam in lams:
                r = L.bb(inst, lam)
                v, n, b, sel = r['best']
                res.append((lam, r['status'], r['M_ub'], n, b, r['nodes'], r['time']))
                fo.write(f"{D}\t{c[0]}\t{c[1]}\t{inst.n}\t{lam:.4f}\t{r['status']}\t"
                         f"{'' if r['M_ub'] is None else round(r['M_ub'], 6)}\t{n}\t{b:.3f}\t{r['nodes']}\t{r['time']:.1f}\t"
                         f"{' '.join(str(inst.P[i]) for i in sel)}\n")
                fo.flush()
                if r['status'] != 'DONE': break
            done = [(lam, M) for lam, st, M, *_ in res if st == 'DONE']
            env = {bb_: min(lam * bb_ + M for lam, M in done) if done else None for bb_ in BGRID}
            kill = all(env[bb_] is not None and env[bb_] < (bb_ + 1) / 3 for bb_ in BGRID if bb_ <= 275)
            print(f"D {D:4d} class {c} pool {inst.n:5d} bad {bad} | " +
                  " ".join(f"l={lam:.3f}:{st[0]}:{'-' if M is None else f'{M:.3f}'}@({b:.1f},{n}),{nd}n,{t:.0f}s"
                           for lam, st, M, n, b, nd, t in res) +
                  f" | UB(60,140,214~220,275) = {[None if env[x] is None else round(env[x], 1) for x in (60, 140, 220, 275)]}"
                  f" | kill1/3 {'YES' if kill else 'NO'}", flush=True)

if __name__ == "__main__":
    main()
