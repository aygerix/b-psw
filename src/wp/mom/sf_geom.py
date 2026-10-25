# sf_geom.py -- audit of the "sign-free" extension of Theorem 10 (geometric steps).
# T1: oval lemma for NEGATIVE quadratics f = -al((t-x)^2 + b^2), al b^2 > 1 (f < -1 on R):
#     {Im t >= 0, |f| <= 1} lies in Im t <= y_top = sqrt(b^2 + 1/al), only top point x + i y_top, and f(top) = +1.
# T2: Step 2 covering for MIXED-sign families without real special points: at the topmost oval top a*
#     (over all members), every member whose top differs from a* has |f_i(a*)| > 1.  Random + adversarial.
# T3: Step 3, all members sharing the top a*: positive members al g - 1, negative members 1 - al g.
#     All-negative: at the "+" point of the largest al, |W| < 1.  Random + adversarial sup search.
# Usage: python3 sf_geom.py SEED
import sys, numpy as np
from scipy.optimize import minimize
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 1)

def quad(s, al, x, b):  # s = leading sign
    return lambda t: s * al * ((t - x)**2 + b**2)
def top(s, al, x, b):
    return x + 1j*np.sqrt(b**2 + 1.0/al)

# ---------- T1 ----------
viol = 0; topval_err = 0.0; nin = 0; maxratio = 0.0
for trial in range(4000):
    al = np.exp(rng.uniform(-3, 3)); x = rng.uniform(-3, 3)
    b = np.sqrt((1.0 + np.exp(rng.uniform(-6, 3))) / al)       # al b^2 > 1
    f = quad(-1, al, x, b); a = top(-1, al, x, b); yt = a.imag
    topval_err = max(topval_err, abs(f(a) - 1.0))               # must be +1 (a "-" special point)
    # sample near the oval: box [x - 2R, x + 2R] x [0, 2 yt]
    R = np.sqrt(np.sqrt(b**4 + 1/al**2)) + yt
    u = rng.uniform(-2*R, 2*R, 200); v = rng.uniform(0, 1.5*yt, 200)
    t = x + u + 1j*v; ins = np.abs(f(t)) <= 1.0
    nin += ins.sum()
    if ins.any():
        maxratio = max(maxratio, (v[ins] / yt).max())
        viol += int(((v[ins] > yt*(1 + 1e-12))).sum())
    # exact check of the defining inequality in (U, v): U^2 + 2U(b^2+v^2) + (b^2-v^2)^2 - c^2 <= 0 at the samples
print(f"T1 negative-quadratic oval: samples inside = {nin}, violations (Im > y_top) = {viol}, "
      f"max Im/y_top inside = {maxratio:.6f}, max |f(top) - 1| = {topval_err:.2e}")

# ---------- T2 ----------
def rand_member(sgn):
    al = np.exp(rng.uniform(-3, 3)); x = rng.uniform(-2, 2)
    b = np.sqrt((1.0 + np.exp(rng.uniform(-6, 3))) / al)
    return (sgn, al, x, b)
bad = 0; mincov = np.inf; tested = 0
for trial in range(40000):
    k = rng.integers(2, 12); sg = rng.choice([-1, 1], size=k)
    if abs(sg.sum()) == k: sg[0] = -sg[0]                       # force mixed signs
    mem = [rand_member(s) for s in sg]
    tops = [top(*m) for m in mem]; js = int(np.argmax([z.imag for z in tops])); a = tops[js]
    fj = quad(*mem[js]); ej = fj(a)                             # +-1
    for i, m in enumerate(mem):
        if i == js: continue
        val = abs(quad(*m)(a)); tested += 1
        mincov = min(mincov, val)
        if val <= 1.0: bad += 1
print(f"T2 mixed-sign Step 2 (random): pairs tested = {tested}, |f_i(a*)| <= 1 cases = {bad}, min |f_i(a*)| = {mincov:.6f}")

# adversarial: two members, i's top strictly below or equal height but different point; minimise |f_i(a*)|
def adv_obj(p, s1, s2):
    # member j*: (s1, al1, 0, b1); member i: (s2, al2, x2, b2) with top height <= top of j*
    la1, lb1, la2, x2, lb2 = p
    al1, al2 = np.exp(la1), np.exp(la2)
    b1 = np.sqrt((1 + np.exp(lb1)) / al1); b2 = np.sqrt((1 + np.exp(lb2)) / al2)
    a = top(s1, al1, 0.0, b1); ti = top(s2, al2, x2, b2)
    if ti.imag > a.imag: return 10.0 + (ti.imag - a.imag)       # infeasible: i would be topmost
    dist = abs(ti - a)
    return abs(quad(s2, al2, x2, b2)(a)) - 1.0 + 0.0*dist
best = []
for (s1, s2) in [(1, -1), (-1, 1), (-1, -1), (1, 1)]:
    bv = np.inf; arg = None
    for rs in range(200):
        p0 = rng.normal(size=5)
        r = minimize(adv_obj, p0, args=(s1, s2), method="Nelder-Mead", options=dict(maxiter=4000, xatol=1e-12, fatol=1e-14))
        if r.fun < bv: bv, arg = r.fun, r.x
    la1, lb1, la2, x2, lb2 = arg; al1, al2 = np.exp(la1), np.exp(la2)
    b1 = np.sqrt((1 + np.exp(lb1)) / al1); b2 = np.sqrt((1 + np.exp(lb2)) / al2)
    gap = abs(top(s2, al2, x2, b2) - top(s1, al1, 0.0, b1))
    best.append((s1, s2, bv, gap))
for s1, s2, bv, gap in best:
    print(f"T2 adversarial signs (j*={s1:+d}, i={s2:+d}): inf(|f_i(a*)| - 1) found = {bv:.3e}, |top_i - a*| there = {gap:.3e}")

# ---------- T3 ----------
# all tops shared, all negative: f_l = 1 - al_l g; at the "+" point of the largest al_i (g = 2/al_i) W = prod (1 - 2 al_l/al_i)
mx = 0.0
for trial in range(100000):
    k = rng.integers(2, 12); al = np.exp(rng.uniform(-4, 4, size=k))
    i = int(np.argmax(al)); W = np.prod([1 - 2*al[l]/al[i] for l in range(k) if l != i])
    mx = max(mx, abs(W))
print(f"T3 all-negative shared top (random, 1e5 families, k <= 11): max |W| = {mx:.6f}")
# adversarial sup of |W| for fixed k (ratios r_l in (0,1)): sup is 1 only at the boundary r -> 0 or r -> 1
for k in (2, 3, 5, 9):
    def negW(z):
        r = 1/(1 + np.exp(-z)); return -abs(np.prod(1 - 2*r))
    bv = 0; bz = None
    for rs in range(100):
        r = minimize(negW, rng.normal(size=k-1)*3, method="Nelder-Mead", options=dict(maxiter=5000))
        if -r.fun > bv: bv, bz = -r.fun, r.x
    rr = 1/(1 + np.exp(-bz))
    print(f"T3 adversarial k={k}: sup |W| found = {bv:.9f}, ratios al_l/al_max = {np.round(rr, 6)}")
# mixed shared top: F(a*) = -1 from a positive member and +1 from a negative one -- check numerically
g = lambda t, x, H: (t - x)**2 + H**2
x, H = 0.3, 1.1; a = x + 1j*H
print(f"T3 mixed shared top: positive member value at a* = {1.7*g(a,x,H) - 1:.3g}, negative member value = {1 - 0.4*g(a,x,H):.3g}")
