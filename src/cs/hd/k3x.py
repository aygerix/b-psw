# HD test 2: k = 3 real Korselt families f_j = 1 + c_j A B_j (A monic deg a = common part of all f_j - 1, B_j monic deg b_j),
# Korselt at roots of B_j (remainder of f_i f_l - 1 mod B_j), Lucas incidences: s_j singles (f_j(s) = -1, f_i f_l(s) = 1)
# and t triples (f_1 = f_2 = f_3 = -1).  Normalisation: A has no x^(a-1) term (translation; if a = 0 fix B_1's), c_1 = +-1.
# For each shape and pattern with G = sum s + 3 t targets, random-start LM; every converged candidate is re-checked from
# scratch (Korselt residual, actual gcd degrees by root matching) and reported if it violates HD (defect < D/2) or
# A2 (defect < 2 d_min).
import numpy as np, itertools, sys, time, json
from numpy.polynomial import polynomial as Pn
from scipy.optimize import least_squares
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
Dmax = int(sys.argv[1]) if len(sys.argv) > 1 else 7
NTR = int(sys.argv[3]) if len(sys.argv) > 3 else 60
def pmul(p, q): return np.convolve(p, q)        # coefficient arrays, highest degree first
def pval(p, x): return np.polyval(p, x)
def build(v, a, b, sgn):
    # unpack: A (a-1 free coeffs below x^(a-1)), B_j (b_j coeffs), c2, c3
    i = 0
    if a >= 1:
        A = np.concatenate(([1.0, 0.0], v[i:i + a - 1])) if a >= 2 else np.array([1.0, 0.0]); i += max(a - 1, 0)
        if a == 1: A = np.array([1.0, 0.0])
    else: A = np.array([1.0])
    Bs = []
    for j in range(3):
        if b[j] == 0: Bs.append(np.array([1.0])); continue
        if a == 0 and j == 0:
            Bj = np.concatenate(([1.0, 0.0], v[i:i + b[j] - 1])); i += b[j] - 1
        else:
            Bj = np.concatenate(([1.0], v[i:i + b[j]])); i += b[j]
        Bs.append(Bj)
    c = [sgn, v[i], v[i + 1]]; i += 2
    f = []
    for j in range(3):
        g = c[j] * pmul(A, Bs[j]); g = g.copy(); g[-1] += 1.0; f.append(g)
    return A, Bs, c, f, i
def nfree(a, b): return max(a - 1, 0) + sum(b) - (1 if (a == 0 and b[0] >= 1) else 0) + 2
def resid(w, a, b, sgn, pat):
    s, t = pat
    A, Bs, c, f, i = build(w, a, b, sgn)
    out = []
    for j in range(3):
        if b[j] == 0: continue
        o = [l for l in range(3) if l != j]
        P = pmul(f[o[0]], f[o[1]]).copy(); P[-1] -= 1.0
        r = np.polydiv(P, Bs[j])[1]
        r = np.concatenate((np.zeros(b[j] - len(r)), r)) if len(r) < b[j] else r[-b[j]:]
        out.extend(list(r))
    sig = w[i:]
    q = 0
    for j in range(3):
        o = [l for l in range(3) if l != j]
        for _ in range(s[j]):
            x = sig[q]; q += 1
            out.append(pval(f[j], x) + 1.0); out.append(pval(f[o[0]], x) * pval(f[o[1]], x) - 1.0)
    for _ in range(t):
        x = sig[q]; q += 1
        for j in range(3): out.append(pval(f[j], x) + 1.0)
    return np.array(out)
def roots_mult(p, tol=1e-5):
    r = np.roots(p); cl = []
    for z in r:
        for C in cl:
            if abs(C[0] - z) < tol * (1 + abs(z)): C[1] += 1; break
        else: cl.append([z, 1])
    return cl
def actual_defect(f):
    F = pmul(pmul(f[0], f[1]), f[2])
    Fm = F.copy(); Fm[-1] -= 1; Fp = F.copy(); Fp[-1] += 1
    # Korselt check: F - 1 mod (f_j - 1) relative residual
    kor = 0.0
    for j in range(3):
        g = f[j].copy(); g[-1] -= 1.0; g = g / g[0]
        r = np.polydiv(Fm / max(1.0, np.max(np.abs(Fm))), g)[1]; kor = max(kor, np.max(np.abs(r)))
    df = []
    for j in range(3):
        h = f[j].copy(); h[-1] += 1.0
        good = 0
        for z, m in roots_mult(h):
            # multiplicity of z in F + 1: count derivatives vanishing
            mm = 0; Q = Fp.copy(); sc = np.max(np.abs(Fp))
            while mm < m and abs(np.polyval(Q, z)) < 1e-7 * sc * (1 + abs(z)) ** len(Q):
                mm += 1; Q = np.polyder(Q)
            good += min(m, mm)
        df.append(len(f[j]) - 1 - good)
    return kor, df
shapes = []
for a in range(0, 4):
    for b in itertools.combinations_with_replacement(range(0, 5), 3):
        D = 3 * a + sum(b)
        if D > Dmax or D < 3: continue
        if a == 0 and b[0] == 0: continue      # constant member
        shapes.append((a, b))
t0 = time.time(); viol = []
for (a, b) in shapes:
    d = [a + bj for bj in b]; D = sum(d); dmin = min(d)
    targets = sorted(set([D // 2 + 1, D - 2 * dmin + 1]))
    nf = nfree(a, b)
    pats = []
    for t in range(0, dmin + 1):
        for s in itertools.product(*[range(0, d[j] - t + 1) for j in range(3)]):
            G = sum(s) + 3 * t
            if G in targets: pats.append((s, t))
    for (s, t) in pats:
        npar = nf + sum(s) + t
        for sgn in (1.0, -1.0):
            for tr in range(NTR):
                w0 = rng.normal(size=npar) * rng.choice([0.7, 1.5, 3.0])
                try: r = least_squares(resid, w0, args=(a, b, sgn, (s, t)), method='trf', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=400)
                except Exception: continue
                if not np.all(np.isfinite(r.fun)) or np.max(np.abs(r.fun)) > 1e-10: continue
                A, Bs, c, f, i = build(r.x, a, b, sgn)
                if min(abs(c[1]), abs(c[2])) < 1e-4 or max(abs(c[1]), abs(c[2])) > 1e4: continue
                L = max(len(x) for x in f)
                pad = [np.concatenate((np.zeros(L - len(x)), x)) for x in f]
                if min(np.max(np.abs(pad[p] - pad[q])) for p in range(3) for q in range(p)) < 1e-4: continue
                sig = r.x[i:]
                if len(sig) > 1 and min(abs(sig[p] - sig[q]) for p in range(len(sig)) for q in range(p)) < 1e-5: continue
                kor, df = actual_defect(f)
                if kor > 1e-7: continue
                tot = sum(df)
                tag = []
                if 2 * tot < D: tag.append("HD-VIOLATION")
                if tot < 2 * dmin: tag.append("A2-VIOLATION")
                if tag:
                    rec = dict(shape=[a, list(b)], deg=d, pat=[list(s), t], sgn=sgn, defects=df, total=tot, D=D, tags=tag,
                               f=[list(map(float, x)) for x in f], sig=list(map(float, sig)), lead=[float(x[0]) for x in f])
                    viol.append(rec); print(json.dumps(rec), flush=True)
                    break
    print(f"shape a={a} b={b} deg={d} D={D}: patterns {len(pats)} done, t={time.time()-t0:.0f}s, violations so far {len(viol)}", flush=True)
print("TOTAL violations", len(viol))
