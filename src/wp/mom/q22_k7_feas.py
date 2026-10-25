"""(2,2) quartic pencils, k = 7: feasibility of the 8 exact boundary equations (q22_k7_B.txt, from q22_moments.gp) in the
actual unknowns r1, r2 in (-1, 1) and w_j = 1/c_j > 0 (7 members), with the reality constraints (K-_j, K+_j have no real
roots). Maximises the minimum reality margin subject to the equations. margin < 0 at every solution => infeasible.
Usage: python3 q22_k7_feas.py starts seed"""
import sys, re, numpy as np
from scipy.optimize import minimize, least_squares
txt = open('q22_k7_B.txt').read().strip().strip('[]')
exprs = [e.strip() for e in re.split(r',(?![^()]*\))', txt)]
# split on top-level commas
parts, depth, cur = [], 0, ''
for ch in txt:
    if ch == '(': depth += 1
    if ch == ')': depth -= 1
    if ch == ',' and depth == 0: parts.append(cur); cur = ''
    else: cur += ch
parts.append(cur)
assert len(parts) == 8, len(parts)
code = [compile(p.strip().replace('^', '**'), 'B', 'eval') for p in parts]
k = 7
def unpack(z):
    r1, r2 = np.tanh(z[0]), np.tanh(z[1]); w = np.exp(z[2:2 + k]); return r1, r2, w
def env(r1, r2, w):
    e1, e2 = r1 + r2, r1 * r2; Dl = (1 - r1**2) * (1 - r2**2)
    d = dict(e1=e1, e2=e2, lam=2 * e1 / Dl, mu=2 * (1 - e1**2 + e2) / Dl)
    for a in range(1, 10): d['p%d' % a] = np.sum(w**a)
    return d
def eqs(z):
    r1, r2, w = unpack(z); d = env(r1, r2, w)
    return np.array([eval(c, {}, d) for c in code])
def margins(z):
    r1, r2, w = unpack(z); d = env(r1, r2, w); lam, mu, e1, e2 = d['lam'], d['mu'], d['e1'], d['e2']; c = 1 / w
    dm = (lam - c * e1)**2 - 4 * c * (c * e2 + mu)          # disc of K-_j = c H + L
    dp = lam**2 - 4 * c * (mu + lam * e1 - c)                 # disc of K+_j = c G + L'
    sep = np.diff(np.sort(z[2:]))                              # distinct members
    return np.concatenate([-dm / (1 + c**2), -dp / (1 + c**2), [np.min(sep) if len(sep) else 1.0]])
starts, seed = int(sys.argv[1]), int(sys.argv[2]); rng = np.random.default_rng(seed)
sols = []
for s in range(starts):
    z0 = np.concatenate([rng.uniform(-2.5, 2.5, 2), rng.uniform(-4, 4, k)])
    ls = least_squares(eqs, z0, method='trf', max_nfev=3000, xtol=1e-14, ftol=1e-14)
    if np.max(np.abs(eqs(ls.x))) > 1e-8: continue
    # maximise min margin along the solution set
    t0 = np.concatenate([ls.x, [np.min(margins(ls.x))]])
    r = minimize(lambda t: -t[-1], t0, method='SLSQP', options={'maxiter': 500, 'ftol': 1e-13},
                 constraints=[{'type': 'eq', 'fun': lambda t: eqs(t[:-1])}, {'type': 'ineq', 'fun': lambda t: margins(t[:-1]) - t[-1]}])
    z = r.x[:-1] if np.max(np.abs(eqs(r.x[:-1]))) < 1e-8 else ls.x
    sols.append((np.min(margins(z)), np.min(margins(ls.x)), z))
sols.sort(key=lambda x: -x[0])
print(f"starts {starts}: solutions of the 8 equations found {len(sols)}")
for m, m0, z in sols[:8]:
    r1, r2, w = unpack(z); mg = margins(z)
    print(f"  best margin {m:+.4e} (at first hit {m0:+.4e})  r = ({r1:+.5f}, {r2:+.5f})  c = {np.round(np.sort(1/w), 5)}"
          f"  K- ok {int(np.sum(mg[:k] > 0))}/7  K+ ok {int(np.sum(mg[k:2*k] > 0))}/7  sep {mg[-1]:.2e}")
