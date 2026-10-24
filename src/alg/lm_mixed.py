"""Numerical search over C for Williams polynomial families with a given degree pattern.
Williams condition (Lemma 19 form): for each j, f_j^2 - 1 divides P_j - 1 with P_j = prod_{i != j} f_i.
Normalisation (removes the affine reparametrisation t -> a t + b): f_1 = t^d1 + 0*t^(d1-1) + lower terms.
Each found solution is classified: MONOMIAL if all members have the same root set (f_j = c_j u^e_j, Conjecture W'),
otherwise NON-MONOMIAL (a counterexample to W').
Usage: python3 lm_mixed.py d1,d2,...,dk seed trials"""
import numpy as np, sys
degs = [int(x) for x in sys.argv[1].split(',')]
rng = np.random.default_rng(int(sys.argv[2]))
trials = int(sys.argv[3])
k = len(degs)
nfree = [degs[0] - 1] + [d + 1 for d in degs[1:]]   # f_1: monic, no t^(d-1) term
nv = sum(nfree)

def polys(z):
    f = []; pos = 0
    for j, d in enumerate(degs):
        if j == 0:
            c = np.concatenate([[1.0 + 0j, 0.0], z[pos:pos + d - 1]]) if d >= 2 else np.array([1.0 + 0j, z[pos] if False else 0.0])
            if d == 1: c = np.array([1.0 + 0j, 0.0])
            pos += d - 1
        else:
            c = z[pos:pos + d + 1]; pos += d + 1
        f.append(np.array(c, dtype=complex))
    return f

def res(z):
    f = polys(z); out = []
    for j in range(k):
        P = np.array([1.0 + 0j])
        for i in range(k):
            if i != j: P = np.polymul(P, f[i])
        P = np.polysub(P, [1])
        M = np.polysub(np.polymul(f[j], f[j]), [1])
        if abs(M[0]) < 1e-12: return np.full(2 * sum(degs), 1e6, dtype=complex)
        M = M / M[0]
        q, r = np.polydiv(P, M)
        L = 2 * degs[j]
        r = np.concatenate([np.zeros(L - len(r)), r]) if len(r) < L else r[-L:]
        out.append(r)
    return np.concatenate(out)

def realres(x):
    z = x[:nv] + 1j * x[nv:]; r = res(z); return np.concatenate([r.real, r.imag])

def lm(x, it=400):
    lam = 1e-3; r = realres(x); c = r @ r
    for _ in range(it):
        J = np.empty((len(r), len(x))); h = 1e-7
        for m in range(len(x)):
            e = np.zeros(len(x)); e[m] = h; J[:, m] = (realres(x + e) - r) / h
        A = J.T @ J; g = J.T @ r
        while True:
            try: dx = np.linalg.solve(A + lam * np.diag(np.diag(A) + 1e-12), -g)
            except np.linalg.LinAlgError: return x, c
            xn = x + dx; rn = realres(xn); cn = rn @ rn
            if cn < c: x, r, c = xn, rn, cn; lam = max(lam / 3, 1e-12); break
            lam *= 4
            if lam > 1e12: return x, c
        if c < 1e-26: break
    return x, c

def classify(f):
    roots = [np.sort_complex(np.roots(fi)) for fi in f]
    base = set(np.round(r, 5) for fi in roots for r in fi)
    same = all(set(np.round(r, 5) for r in ri) == set(np.round(r, 5) for r in roots[0]) for ri in roots)
    return "MONOMIAL" if same else "NON-MONOMIAL"

found = {"MONOMIAL": 0, "NON-MONOMIAL": 0}; best = []
for trial in range(trials):
    x = rng.normal(size=2 * nv) * rng.choice([0.5, 1, 2])
    x, c = lm(x)
    z = x[:nv] + 1j * x[nv:]
    f = polys(z)
    lead = min(abs(fi[0]) for fi in f)
    dist = min((np.linalg.norm(f[i] - f[j]) for i in range(k) for j in range(i + 1, k) if degs[i] == degs[j]), default=1)
    best.append(c)
    if c < 1e-18 and lead > 1e-4:
        kind = classify(f); found[kind] += 1
        if kind == "NON-MONOMIAL" or found[kind] <= 2:
            print(kind, "c=%.1e" % c, "distinct" if dist > 1e-4 else "repeated", [np.round(fi, 5).tolist() for fi in f], flush=True)
print("pattern", degs, "trials", trials, "min residual %.2e" % min(best), "solutions", found, flush=True)
