# HD at k = 3: for every shape (a; b1 <= b2 <= b3), a >= 1, not a pencil, D = 3a + sum b <= Dmax:
#   U_j = min(d_j, deg-bound(L_j), D - 2 d_j)  (A-reduction lemma + (f_j - 1) h_j | P_j - 1).
#   If sum U_j <= D/2: HD holds for the shape by counting (PROVED).  Otherwise list every target vector g (g_j <= U_j,
#   sum g = floor(D/2) + 1, up to the symmetry of equal b's) for the exact Singular test (gen_k3h.py).
import itertools, sys
Dmax = int(sys.argv[1])
def degL(b, j):
    i, l = [x for x in range(3) if x != j]
    return max(b[j] + max(b[i], b[l]), b[i] + b[l])
proved = []; todo = []
for a in range(1, Dmax):
    for b in itertools.combinations_with_replacement(range(0, Dmax), 3):
        D = 3*a + sum(b)
        if D > Dmax or sum(b) == 0: continue
        d = [a + x for x in b]
        U = [max(0, min(d[j], degL(b, j), D - 2*d[j])) for j in range(3)]
        if 2*sum(U) <= D: proved.append((a, b, U)); continue
        T = D//2 + 1; seen = set(); G = []
        for g in itertools.product(*[range(0, U[j]+1) for j in range(3)]):
            if sum(g) != T: continue
            key = tuple(sorted(zip(b, g)))
            if key in seen: continue
            seen.add(key); G.append(g)
        todo.append((a, b, U, G))
print(f"# Dmax={Dmax}: shapes proved by counting: {len(proved)}; shapes needing exact tests: {len(todo)}")
for a, b, U in proved: print(f"COUNT a={a} b={b} U={U} sumU={sum(U)} D={3*a+sum(b)}")
specs = []
for a, b, U, G in todo:
    print(f"TEST a={a} b={b} U={U} D={3*a+sum(b)} targets={G}")
    for g in G: specs.append(f"{a}/{''.join(map(str,b))}/{''.join(map(str,g))}")
open(f"hd3_specs_{Dmax}.txt", "w").write(";".join(specs))
print("# total exact jobs:", len(specs))
