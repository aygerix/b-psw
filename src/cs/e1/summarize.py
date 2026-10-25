"""Summarise run_bbc.sh output: per instance the certified envelope UB(b) = min over finished lambda of
lambda*b + M_ub(lambda), the kill check UB(b) < (b+1)/3 for b in BGRID (b <= 275), and aggregates.
Usage: python3 summarize.py OUTFILE [more OUTFILES]   (later files override earlier ones per instance and lambda)"""
import sys
from collections import defaultdict

BGRID = [60, 80, 100, 120, 140, 160, 180, 200, 214, 220, 240, 260, 275, 300]

def parse(files):
    inst = {}
    for fn in files:
        cur = None
        for l in open(fn):
            a = l.split()
            if not a: continue
            if a[0] == '==':
                cur = (int(a[1]), int(a[2]), int(a[3])); inst.setdefault(cur, dict(n=int(a[4]), lam={}, rc=None))
            elif a[0] == 'lam' and cur:
                lam = float(a[1]); st = a[2]
                d = dict(status=st, M=None if st != 'DONE' else float(a[4]), best=float(a[6]), N=int(a[8]),
                         bits=float(a[10]), nodes=int(a[12]), time=float(a[18]))
                inst[cur]['lam'][lam] = d
            elif a[0] == 'rc' and cur:
                inst[cur]['rc'] = int(a[1])
    return inst

def main():
    inst = parse(sys.argv[1:])
    worst = {b: (0, None) for b in BGRID}
    notkilled, nodone = [], []
    lamstat = defaultdict(lambda: [0, 0])
    for key, v in sorted(inst.items()):
        done = [(lam, d['M']) for lam, d in v['lam'].items() if d['status'] == 'DONE']
        for lam, d in v['lam'].items(): lamstat[lam][d['status'] == 'DONE'] += 1
        if not done: nodone.append((key, v['n'])); continue
        env = {b: min(lam * b + M for lam, M in done) for b in BGRID}
        for b in BGRID:
            r = env[b] / (b + 1)
            if r > worst[b][0]: worst[b] = (r, key, env[b])
        if not all(env[b] < (b + 1) / 3 for b in BGRID if b <= 275): notkilled.append((key, v['n'], env))
    print(f"instances parsed: {len(inst)}; with >= 1 finished lambda: {len(inst) - len(nodone)}; none finished: {len(nodone)}")
    print("lambda: finished / not finished:", {k: tuple(reversed(v)) for k, v in sorted(lamstat.items(), reverse=True)})
    print("worst certified UB(b)/(b+1) over all instances with a finished lambda:")
    for b in BGRID:
        r, key, e = worst[b]
        print(f"   b = {b:3d}: {r:.4f}  (UB {e:.2f}, instance D,j1,j2 = {key})")
    print(f"instances NOT killed at 1/3 (b <= 275): {len(notkilled)}")
    for key, n, env in notkilled[:30]: print("   ", key, n, {b: round(env[b], 1) for b in (60, 140, 214, 275)})
    print(f"instances with no finished lambda: {nodone[:40]}")
    # best (smallest) lambda finished per D, for the big classes
    print("largest classes: finished lambdas and M_ub, hull points (bits, N):")
    for key, v in sorted(inst.items(), key=lambda t: -t[1]['n'])[:40]:
        print(f"   D {key[0]:4d} ({key[1]},{key[2]}) pool {v['n']:5d}: " + "; ".join(
            f"l={lam:.3f} {d['status'][0]} M={d['M'] if d['M'] is None else round(d['M'], 3)} ({d['bits']:.1f},{d['N']}) "
            f"{d['nodes']}n {d['time']:.0f}s" for lam, d in sorted(v['lam'].items(), reverse=True)))

if __name__ == "__main__":
    main()
