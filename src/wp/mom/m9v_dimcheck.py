#!/usr/bin/env python3
"""m9v (verifier): INDEPENDENT re-implementation of the m9 Krawczyk certificates and of the local-dimension claim.

No code is shared with cubic_moments.gp / m9_kraw.gp.  The four boundary equations are evaluated directly from the
truncated Laurent series (no symbolic polynomial expansion):
   x_i = v s^3 / (1 + (1+u) s + (u+2v) s^2),  linear member x_0 = s/2,  G = prod (1 + x_i) mod s^k,
   V(t) = sum_{j<k} G_j t^(k-1-j),  E1 = V(-1) - 2^(1-k),  E2 = V'(-1) - (1-k) 2^(1-k),  E3, E4: same with u -> 1-u.
Jacobians by forward-mode automatic differentiation.  Rigorous part: exact rational interval arithmetic (Fractions),
so the Krawczyk verdict involves no floating point.  Floating point (mpmath) is used only to produce the centre and
the preconditioner Y, which are then rationalised.

Usage:
  m9v_dimcheck.py gen  <m> <lin> <solfile> <unknowns comma-separated, e.g. u9,v1,v5,v9>
  m9v_dimcheck.py sym  <solfile>          (reflection-symmetric slice, 11 cubics = 5 pairs + centre)
"""
import sys
import itertools
from fractions import Fraction as Fr
import mpmath

mpmath.mp.dps = 220


class Iv:
    __slots__ = ('lo', 'hi')

    def __init__(self, lo, hi=None):
        self.lo = lo
        self.hi = lo if hi is None else hi

    def __add__(self, o):
        o = toiv(o)
        return Iv(self.lo + o.lo, self.hi + o.hi)

    __radd__ = __add__

    def __neg__(self):
        return Iv(-self.hi, -self.lo)

    def __sub__(self, o):
        o = toiv(o)
        return Iv(self.lo - o.hi, self.hi - o.lo)

    def __rsub__(self, o):
        return toiv(o) - self

    def __mul__(self, o):
        o = toiv(o)
        p = (self.lo * o.lo, self.lo * o.hi, self.hi * o.lo, self.hi * o.hi)
        return Iv(min(p), max(p))

    __rmul__ = __mul__

    def mag(self):
        return max(abs(self.lo), abs(self.hi))


def toiv(x):
    return x if isinstance(x, Iv) else Iv(x)


class Du:
    """value + gradient (generic base type: Iv, Fraction or mpf)"""
    __slots__ = ('a', 'g')

    def __init__(self, a, g):
        self.a = a
        self.g = g

    def __add__(self, o):
        if isinstance(o, Du):
            return Du(self.a + o.a, [x + y for x, y in zip(self.g, o.g)])
        return Du(self.a + o, list(self.g))

    __radd__ = __add__

    def __neg__(self):
        return Du(-self.a, [-x for x in self.g])

    def __sub__(self, o):
        return self + (-o)

    def __rsub__(self, o):
        return (-self) + o

    def __mul__(self, o):
        if isinstance(o, Du):
            return Du(self.a * o.a, [self.a * y + o.a * x for x, y in zip(self.g, o.g)])
        return Du(self.a * o, [x * o for x in self.g])

    __rmul__ = __mul__


def boundary(members, lin, k, zero, one, cst):
    """members: list of (u, v) (Du objects); returns [E1, E2, E3, E4]."""
    out = []
    for refl in (0, 1):
        G = [one] + [zero] * (k - 1)
        for (u, v) in members:
            uu = (one - u) if refl else u
            p = one + uu
            q = uu + 2 * v
            c = [one, -p]
            while len(c) < k:
                c.append(-p * c[-1] - q * c[-2])
            x = [zero] * k
            for j in range(3, k):
                x[j] = v * c[j - 3]
            H = [zero] * k
            for a in range(k):
                for b in range(k - a):
                    if b == 0:
                        H[a] = H[a] + G[a]
                    else:
                        H[a + b] = H[a + b] + G[a] * x[b]
            G = H
        if lin:
            H = list(G)
            for a in range(k - 1):
                H[a + 1] = H[a + 1] + G[a] * cst(Fr(1, 2))
            G = H
        V0 = zero
        V1 = zero
        for j in range(k):
            e = k - 1 - j
            V0 = V0 + G[j] * ((-1) ** e)
            if e >= 1:
                V1 = V1 + G[j] * (e * (-1) ** (e - 1))
        out += [V0 - cst(Fr(1, 2 ** (k - 1))), V1 - cst(Fr(1 - k, 2 ** (k - 1)))]
    return out


def ratround(x, D):
    return Fr(int(mpmath.nint(mpmath.mpf(x) * D)), D)


def mkzero(n, base):
    return Du(base(0), [base(0)] * n)


def krawczyk(build, n, x0, label):
    """build(X, base) -> list of n equations (Du with n-gradient) for unknown vector X (list of base numbers)."""
    # Newton in mpmath
    xt = [mpmath.mpf(x) for x in x0]
    for _ in range(60):
        E = build([Du(xt[j], [mpmath.mpf(1) if i == j else mpmath.mpf(0) for i in range(n)]) for j in range(n)], 'mp')
        Fv = mpmath.matrix([e.a for e in E])
        Jm = mpmath.matrix([[e.g[c] for c in range(n)] for e in E])
        dx = mpmath.lu_solve(Jm, Fv)
        xt = [xt[j] - dx[j] for j in range(n)]
    print(f"  [{label}] Newton correction from the numerical start: {mpmath.nstr(max(abs(xt[j] - x0[j]) for j in range(n)), 5)}")
    xq = [ratround(x, 10 ** 110) for x in xt]
    E = build([Du(mpmath.mpf(xq[j].numerator) / xq[j].denominator,
                  [mpmath.mpf(1) if i == j else mpmath.mpf(0) for i in range(n)]) for j in range(n)], 'mp')
    Jm = mpmath.matrix([[e.g[c] for c in range(n)] for e in E])
    print(f"  [{label}] det of the {n}x{n} Jacobian minor at the centre (float): {mpmath.nstr(mpmath.det(Jm), 12)}")
    Yf = Jm ** -1
    Y = [[ratround(Yf[r, c], 10 ** 40) for c in range(n)] for r in range(n)]
    # exact F at the centre
    Eq = build([Du(xq[j], [Fr(0)] * n) for j in range(n)], 'q')
    Fq = [e.a for e in Eq]
    YF = [sum(Y[r][l] * Fq[l] for l in range(n)) for r in range(n)]
    rho = Fr(1, 10 ** 60)
    B = [Iv(xq[j] - rho, xq[j] + rho) for j in range(n)]
    EB = build([Du(B[j], [Iv(Fr(1 if i == j else 0)) for i in range(n)]) for j in range(n)], 'iv')
    ok = True
    rows = []
    for r in range(n):
        sb = Fr(0)
        for c in range(n):
            acc = Iv(Fr(1 if r == c else 0))
            for l in range(n):
                acc = acc + Iv(-Y[r][l]) * EB[l].g[c]
            sb += acc.mag()
        rows.append(sb)
        if not (abs(YF[r]) + rho * sb < rho):
            ok = False
    print(f"  [{label}] max |Y F(x~)| = {float(max(abs(z) for z in YF)):.3e}, rho = 1e-60, "
          f"row sums |I - Y J(B)|: {[float(z) for z in rows]}")
    print(f"  [{label}] KRAWCZYK {'PASSED (unique zero in B; ||I - Y J(B)||_inf < 1 so J nonsingular on B)' if ok else 'FAILED'}")
    return ok, xq, B, rho


def readsol(fn):
    env = {}
    for line in open(fn):
        line = line.strip()
        if line.startswith('\\\\') or not line:
            continue
        name, val = line.rstrip(';').split('=', 1)
        env[name.strip()] = [mpmath.mpf(z) for z in val.strip().strip('[]').split(',')]
    return env


def domain_and_distinct(mem_iv):
    """mem_iv: list of (Iv u, Iv v). exact checks u^2 < 8v, (1-u)^2 < 8v on the box, pairwise distinct."""
    dom = True
    minm = None
    for (u, v) in mem_iv:
        u2 = max(u.lo ** 2, u.hi ** 2)
        w2 = max((1 - u.lo) ** 2, (1 - u.hi) ** 2)
        v8 = 8 * v.lo
        if not (u2 < v8 and w2 < v8 and v.lo > 0):
            dom = False
        mg = v8 - max(u2, w2)
        minm = mg if minm is None else min(minm, mg)
    dist = True
    for i, j in itertools.combinations(range(len(mem_iv)), 2):
        (ui, vi), (uj, vj) = mem_iv[i], mem_iv[j]
        if not (ui.hi < uj.lo or uj.hi < ui.lo or vi.hi < vj.lo or vj.hi < vi.lo):
            dist = False
    return dom, dist, minm


def base_of(tag):
    return {'mp': mpmath.mpf, 'q': Fr, 'iv': (lambda z: Iv(Fr(z)))}[tag]


def cst_of(tag):
    return {'mp': (lambda q: mpmath.mpf(q.numerator) / q.denominator), 'q': (lambda q: q), 'iv': (lambda q: Iv(q))}[tag]


def run_gen(m, lin, solfile, unk):
    k = m + lin
    sol = readsol(solfile)
    U0, V0 = sol['U0'], sol['V0']
    names = [f"u{i + 1}" for i in range(m)] + [f"v{i + 1}" for i in range(m)]
    X0 = list(U0) + list(V0)
    sel = [names.index(z) for z in unk]
    fix = [i for i in range(2 * m) if i not in sel]
    fq = {i: ratround(X0[i], 10 ** 8) for i in fix}
    print(f"==== k = {k}: {m} cubics, linear member {lin}; unknowns {unk}")
    print("  fixed rationals:", [str(fq[i]) for i in fix])

    def build(X, tag):
        base = base_of(tag)
        n = 4
        coord = {}
        cst = cst_of(tag)
        for i in fix:
            coord[i] = Du(cst(fq[i]), [base(0)] * n)
        for j, i in enumerate(sel):
            coord[i] = X[j]
        mem = [(coord[i], coord[m + i]) for i in range(m)]
        zero = mkzero(n, base)
        one = Du(base(1), [base(0)] * n)
        return boundary(mem, lin, k, zero, one, cst)

    ok, xq, B, rho = krawczyk(build, 4, [X0[i] for i in sel], f"k={k},m={m}")
    coordiv = {}
    for i in fix:
        coordiv[i] = Iv(fq[i])
    for j, i in enumerate(sel):
        coordiv[i] = B[j]
    dom, dist, minm = domain_and_distinct([(coordiv[i], coordiv[m + i]) for i in range(m)])
    print(f"  domain on box: {'HOLD' if dom else 'FAIL'} (min margin {float(minm):.6g}); distinct: {'YES' if dist else 'NO'}")
    # symmetry diagnostic of the certified zero (float)
    xc = {i: float(fq[i]) for i in fix}
    for j, i in enumerate(sel):
        xc[i] = float(xq[j])
    print("  centre of unknowns:", [f"{float(z):.15f}" for z in xq])
    print(f"  VERDICT {'VALID' if ok and dom and dist else 'NOT VALID'}")
    return ok and dom and dist


def run_sym(solfile):
    sol = readsol(solfile)
    U0, V0 = sol['U0'], sol['V0']
    m = 11
    p = 5
    W0 = [U0[j] - mpmath.mpf(1) / 2 for j in range(p)]
    PV0 = [V0[j] for j in range(p)]
    vc0 = V0[10]
    # numeric symmetry of the seed
    asym = max([abs(U0[j] + U0[j + p] - 1) for j in range(p)] + [abs(V0[j] - V0[j + p]) for j in range(p)] + [abs(U0[10] - 0.5)])
    print(f"==== symmetric slice, 11 cubics (5 pairs + centre u = 1/2); seed asymmetry {mpmath.nstr(asym, 3)}")
    names = [f"w{j + 1}" for j in range(p)] + [f"v{j + 1}" for j in range(p)] + ["vc"]
    X0 = W0 + PV0 + [vc0]

    def mk(sel, fq):
        def build(X, tag):
            base = base_of(tag)
            n = len(sel)
            coord = {}
            cst = cst_of(tag)
            for i in range(len(names)):
                if i in sel:
                    coord[i] = X[sel.index(i)]
                else:
                    coord[i] = Du(cst(fq[i]), [base(0)] * n)
            half = cst(Fr(1, 2))
            mem = [(coord[j] + half, coord[p + j]) for j in range(p)] + \
                  [((-coord[j]) + half, coord[p + j]) for j in range(p)] + \
                  [(Du(half, [base(0)] * n), coord[2 * p])]
            zero = mkzero(n, base)
            one = Du(base(1), [base(0)] * n)
            E = boundary(mem, 0, 11, zero, one, cst)
            return E[:2]
        return build

    # choose 2 unknowns maximising |det| of the numerical 2x2 minor
    fq0 = {i: ratround(X0[i], 10 ** 8) for i in range(len(names))}
    allsel = list(range(len(names)))
    fullb = mk(allsel, fq0)
    nn = len(names)
    E = fullb([Du(X0[j], [mpmath.mpf(1) if i == j else mpmath.mpf(0) for i in range(nn)]) for j in range(nn)], 'mp')
    print("  residuals of E1, E2 at seed:", [mpmath.nstr(e.a, 3) for e in E])
    # also the reflected pair at the seed, full 4 equations (must coincide with E1, E2)
    best = None
    for a, b in itertools.combinations(range(nn), 2):
        d = abs(E[0].g[a] * E[1].g[b] - E[0].g[b] * E[1].g[a])
        if best is None or d > best[0]:
            best = (d, [a, b])
    sel = best[1]
    print(f"  unknowns {[names[i] for i in sel]}, |det| numerical minor {mpmath.nstr(best[0], 6)}")
    fq = {i: fq0[i] for i in range(nn) if i not in sel}
    print("  fixed rationals:", {names[i]: str(fq[i]) for i in fq})
    ok, xq, B, rho = krawczyk(mk(sel, fq), 2, [X0[i] for i in sel], "sym11")
    co = {}
    for i in range(nn):
        co[i] = B[sel.index(i)] if i in sel else Iv(fq[i])
    half = Fr(1, 2)
    mem = [(co[j] + half, co[p + j]) for j in range(p)] + [((-co[j]) + half, co[p + j]) for j in range(p)] + [(Iv(half), co[2 * p])]
    dom, dist, minm = domain_and_distinct(mem)
    print(f"  domain on box: {'HOLD' if dom else 'FAIL'} (min margin {float(minm):.6g}); distinct: {'YES' if dist else 'NO'}")
    # all four equations at the exact centre (E3 = E1, E4 = E2 must hold identically on the slice)
    cq = {i: (xq[sel.index(i)] if i in sel else fq[i]) for i in range(nn)}
    memq = [(Du(cq[j] + half, []), Du(cq[p + j], [])) for j in range(p)] + \
           [(Du(half - cq[j], []), Du(cq[p + j], [])) for j in range(p)] + [(Du(half, []), Du(cq[2 * p], []))]
    E4 = boundary(memq, 0, 11, Du(Fr(0), []), Du(Fr(1), []), cst_of('q'))
    print("  exact check at the centre: E3 - E1 == 0:", E4[2].a - E4[0].a == 0, ", E4 - E2 == 0:", E4[3].a - E4[1].a == 0,
          f", |E1|, |E2| at centre: {float(abs(E4[0].a)):.2e}, {float(abs(E4[1].a)):.2e}")
    # rank of the full 4 x 22 Jacobian at the (float) symmetric point
    uu = [mpmath.mpf(cq[j].numerator) / cq[j].denominator for j in range(nn)]
    full = [uu[j] + 0.5 for j in range(p)] + [0.5 - uu[j] for j in range(p)] + [mpmath.mpf(0.5)] + uu[p:2 * p] + uu[p:2 * p] + [uu[2 * p]]
    X = [Du(full[j], [mpmath.mpf(1) if i == j else mpmath.mpf(0) for i in range(22)]) for j in range(22)]
    Ef = boundary([(X[i], X[11 + i]) for i in range(11)], 0, 11, mkzero(22, mpmath.mpf), Du(mpmath.mpf(1), [mpmath.mpf(0)] * 22), cst_of('mp'))
    Jf = mpmath.matrix([[e.g[c] for c in range(22)] for e in Ef])
    sv = mpmath.svd_r(Jf, compute_uv=False)
    print("  singular values of the full 4 x 22 Jacobian at the symmetric point (float):", [mpmath.nstr(z, 6) for z in sv])
    print(f"  VERDICT {'VALID: symmetric slice meets the feasible set (11 cubics)' if ok and dom and dist else 'NOT VALID'}")


if __name__ == "__main__":
    if sys.argv[1] == 'gen':
        run_gen(int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5].split(','))
    else:
        run_sym(sys.argv[2])
