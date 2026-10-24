"""Driver: all k = 7 axis-class structures left by Lemma 23 (end classes singletons; middle as listed), real LM each.
B-type middle classes of size 5: every bipartite simple graph with 5 edges (members) between a "-" heights and b "+"
heights, no isolated heights, at least one shared "-" and one shared "+" height, up to relabelling.
Usage: python3 rq_k7.py starts
"""
import sys, itertools, subprocess, os

def canon(edges, a, b):
    best = None
    for pa in itertools.permutations(range(a)):
        for pb in itertools.permutations(range(b)):
            e = tuple(sorted((pa[x], pb[y]) for x, y in edges))
            if best is None or e < best: best = e
    return best

def patterns5():
    out = set()
    for a in range(1, 5):
        for b in range(1, 5):
            allE = [(x, y) for x in range(a) for y in range(b)]
            for E in itertools.combinations(allE, 5):
                if {x for x, _ in E} != set(range(a)) or {y for _, y in E} != set(range(b)): continue
                dl = [sum(1 for x, _ in E if x == i) for i in range(a)]
                dr = [sum(1 for _, y in E if y == j) for j in range(b)]
                if max(dl) < 2 or max(dr) < 2: continue
                out.add((a, b, canon(E, a, b)))
    return sorted(out)

def mstr(E): return ",".join(f"m{x + 1}p{y + 1}" for x, y in E)

S = "m1p1"
structs = [("A_a [1|2-|2+|1|1]", f"{S}|m1p1,m1p2|m1p1,m2p1|{S}|{S}"),
           ("A_b [1|2-|1|2+|1]", f"{S}|m1p1,m1p2|{S}|m1p1,m2p1|{S}"),
           ("A_c [1|2+|2-|1|1]", f"{S}|m1p1,m2p1|m1p1,m1p2|{S}|{S}"),
           ("C [1|1|3path|1|1]", f"{S}|{S}|m1p1,m2p1,m2p2|{S}|{S}")]
for a, b, E in patterns5():
    structs.append((f"B [1|5:{a}x{b}|1]", f"{S}|{mstr(E)}|{S}"))
starts = sys.argv[1] if len(sys.argv) > 1 else "100"
here = os.path.dirname(os.path.abspath(__file__))
print(f"{len(structs)} structures", flush=True)
for name, st in structs:
    print("==", name, flush=True)
    r = subprocess.run([sys.executable, os.path.join(here, "rq_struct.py"), st, starts, "7"], capture_output=True, text=True)
    print(r.stdout.strip(), flush=True)
