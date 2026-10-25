"""FX seeded variant of ls300.py (see the comment above the seeded search).  Usage as ls300.py plus env SEEDCFG, MUS.
Adversarial local search (falsification attempt) for the E1 certificate, independent code.
Maximises N(C) - mu*bits(C) over configurations (each row: off, or side F/L with cap e) by iterated local search with
exact row-state moves; reports (a) for each lambda in LAMS the best N - lambda*bits seen (must be <= the certified M_ub),
(b) the primal frontier: best N seen with bits <= b (a lower bound on N_max(b)).
Also climbs from the minimal configuration of each of S randomly chosen single primes (S = -1: every prime).
Usage: python3 aud_ls.py D j1 j2 LAMS restarts seed S recfile1 [recfile2 ...]"""
import sys, math, random
import numpy as np
D, j1, j2 = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
LAMS = [float(x) for x in sys.argv[4].split(',')]; R = int(sys.argv[5]); seed = int(sys.argv[6]); S = int(sys.argv[7]); files = sys.argv[8:]
import os
T0 = int(os.environ.get("T0", 200))
PR = [q for q in range(3, T0 + 1) if all(q % d for d in range(2, int(q ** .5) + 1))]
def fac(x):
    f = {}
    for q in PR:
        while x % q == 0: x //= q; f[q] = f.get(q, 0) + 1
    assert x == 1
    return f
pool = {}
for fn in files:
    for l in open(fn):
        t = tuple(map(int, l.split()))
        if len(t) == 7 and t[0] == D and (t[3], t[4]) == (j1, j2): pool[t[1]] = (t[5], t[6])
P = sorted(pool); nP = len(P)
rows = sorted({q for p in P for q in list(fac(pool[p][0])) + list(fac(pool[p][1]))}); ri = {q: k for k, q in enumerate(rows)}
nR = len(rows)
need = np.zeros((2, nR, nP), dtype=np.int16)
for k, p in enumerate(P):
    for s in (0, 1):
        for q, a in fac(pool[p][s]).items(): need[s, ri[q], k] = a
E = need.max(axis=2)                                   # E[s, r]
def cst(q, e): return 0.0 if e == 0 else math.log2(q - 1) + (e - 1) * math.log2(q)
states = []                                            # per row: list of (s, e, bits, violvec)
for r, q in enumerate(rows):
    st = [(None, 0, 0.0, ((need[0, r] > 0).astype(np.int16) + (need[1, r] > 0)))]
    for s in (0, 1):
        for e in range(1, E[s, r] + 1):
            st.append((s, e, cst(q, e), ((need[s, r] > e).astype(np.int16) + (need[1 - s, r] > 0))))
    states.append(st)
SV = [np.stack([t[3] for t in st]) for st in states]  # (nstates, nP)
SB = [np.array([t[2] for t in st]) for st in states]
rng = random.Random(seed)
best = {lam: (-1e9, None) for lam in LAMS}
front = {}
def note(N, bits, cfg):
    for lam in LAMS:
        v = N - lam * bits
        if v > best[lam][0] + 1e-12: best[lam] = (v, (N, round(bits, 3), tuple(cfg)))
    b = int(math.ceil(bits - 1e-9))
    if N > front.get(b, (-1,))[0]: front[b] = (N, round(bits, 3), tuple(cfg))
def climb(cfg, mu):
    viol = sum(SV[r][cfg[r]] for r in range(nR)); bits = sum(SB[r][cfg[r]] for r in range(nR))
    improved = True
    while improved:
        improved = False
        for r in rng.sample(range(nR), nR):
            base = viol - SV[r][cfg[r]]
            Ns = ((base[None, :] + SV[r]) == 0).sum(axis=1)
            bs = bits - SB[r][cfg[r]] + SB[r]
            for k in range(len(Ns)): note(int(Ns[k]), float(bs[k]), cfg[:r] + [k] + cfg[r + 1:])
            obj = Ns - mu * bs
            k = int(np.argmax(obj))
            if obj[k] > obj[cfg[r]] + 1e-12:
                cfg[r] = k; viol = base + SV[r][k]; bits = float(bs[k]); improved = True
    return cfg
def single(k):
    cfg = []
    for r in range(nR):
        st = 0
        for idx, t in enumerate(states[r]):
            if t[0] is not None and need[t[0], r, k] == t[1] and need[1 - t[0], r, k] == 0: st = idx
        cfg.append(st)
    return cfg
# FX seeded search: start from a given configuration SEEDCFG = "q:s:e,q:s:e,..." (s = F or L), report the marginal gain of
# every single extra row item (q, s, e = 1) not in the seed, then climb from the seed (and R perturbed copies) at each mu in
# LAMS + [MUS...].  Everything else (states, note(), climb(), frontier, DUMP) is the ls300.py code.
seedspec = [t.split(':') for t in os.environ['SEEDCFG'].split(',')]
cfg0 = [0] * nR
for q, sd, e in seedspec:
    r = ri[int(q)]; sd = 0 if sd == 'F' else 1; e = int(e)
    k = next(idx for idx, t in enumerate(states[r]) if t[0] == sd and t[1] == e); cfg0[r] = k
viol0 = sum(SV[r][cfg0[r]] for r in range(nR)); N0 = int((viol0 == 0).sum()); b0 = float(sum(SB[r][cfg0[r]] for r in range(nR)))
note(N0, b0, cfg0)
print(f"seed: N {N0} bits {b0:.3f}")
gains = []
for r in range(nR):
    if cfg0[r]: continue
    for k, t in enumerate(states[r]):
        if t[0] is None or t[1] != 1: continue
        v = viol0 - SV[r][0] + SV[r][k]; N1 = int((v == 0).sum()); gains.append((N1 - N0, rows[r], 'FL'[t[0]], round(t[2], 3)))
gains.sort(reverse=True)
print("best single-row additions (gain, q, side, bits):", gains[:12])
MUS = [float(x) for x in os.environ.get('MUS', '0.2,0.22,0.25').split(',')]
for mu in MUS:
    for it in range(R + 1):
        c = list(cfg0)
        if it:
            for r in rng.sample(range(nR), rng.randint(1, 4)): c[r] = rng.randrange(len(states[r]))
        climb(c, mu)
print(f"D = {D}, class ({j1},{j2}): pool {nP}, rows {nR}, restarts {R}, seed {seed}")
for lam in LAMS:
    v, (N, b, cfg) = best[lam]
    desc = [(rows[r], 'FL'[states[r][k][0]], states[r][k][1]) for r, k in enumerate(cfg) if k]
    print(f"lam {lam}: best N - lam*bits = {v:.5f} at (N {N}, bits {b}); config {desc}")
import os
DUMP = [int(x) for x in os.environ.get('DUMP', '').split(',') if x]
run = (-1, 0); fr = {}
for b in range(0, 401):
    if b in front and front[b][0] > run[0]: run = front[b]
    fr[b] = run
print("primal frontier (b: N@bits):", " ".join(f"{b}:{fr[b][0]}@{fr[b][1]}" for b in (60, 100, 140, 180, 214, 227, 250, 275, 300, 350, 400)))
for b in DUMP:
    if len(fr[b]) != 3: continue          # no configuration with bits <= b was seen
    N, bits, cfg = fr[b]
    cap = np.zeros((2, nR), dtype=np.int16)
    for r, k in enumerate(cfg):
        s_, e_ = states[r][k][0], states[r][k][1]
        if s_ is not None: cap[s_, r] = e_
    ok = np.all(need <= cap[:, :, None], axis=(0, 1))
    used = [P[k] for k in range(nP) if ok[k]]
    E1 = [(rows[r], int(cap[0, r])) for r in range(nR) if cap[0, r]]; E2 = [(rows[r], int(cap[1, r])) for r in range(nR) if cap[1, r]]
    print(f"DUMP b<={b}: N {N} bits {bits} E1 {E1} E2 {E2} primes {used}")
