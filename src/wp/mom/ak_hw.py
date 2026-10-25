# ak_hw.py -- test the H_w extremal principle for real cubic members.
# Coordinates V = 1/a, ut = c/a - 1/2.  Field n^w = (-3/2 - ut, -rho), rho = 2V + 1/8 - ut^2/2 (exact: Re f_Q(w_P) - 1 = a_Q <n^w_P, Q - P>).
# H_w has level curves = integral curves of dut/dV = -(3/2 + ut)/rho; H_w increases toward -n^w.
# (1) local condition L_w(V, ut) >= 0 over the domain;  (2) GLOBAL test: for random admissible P, sample admissible Q in the
#     covering ellipse E_P^w = {Q : |f_Q(w_P)| < 1} and check Q lies strictly on the high side of the level curve through P.
import sys, numpy as np
from scipy.integrate import solve_ivp
from ak_core import special, fval
rng = np.random.default_rng(int(sys.argv[1])); NP = int(sys.argv[2]); NQ = int(sys.argv[3]); which = sys.argv[4]  # 'w' or 'v'
sgn = 1.0 if which == 'w' else -1.0     # v-case = reflection ut -> -ut of the w-case
def Lw(V, ut):
    return 4*(2*V - (0.25 - ut/2)**2)*(2*V + 0.5 + ut)**2 - V*(1.5 + ut)*(6*V + 3/8 + 1.5*ut - ut**2/2)
def admissible(V, ut):
    return (V > 1/32) & ((np.abs(ut) + 0.5)**2 < 8*V)
# (1) local condition map on a grid of the domain (in the w-frame; the v-frame is ut -> -ut)
Vs = np.exp(np.linspace(np.log(1/32), np.log(20), 400)); fails = []
for V in Vs:
    U = np.linspace(-(np.sqrt(8*V) - 0.5), np.sqrt(8*V) - 0.5, 401)[1:-1]
    L = Lw(V, sgn*U)
    if np.any(L < 0): fails.append((1/V, U[L < 0].min(), U[L < 0].max()))
print("local condition fails (a, ut-range):", [(round(f[0], 3), round(f[1], 3), round(f[2], 3)) for f in fails[::max(1, len(fails)//15)]], " total rows", len(fails))
amax_fail = max([f[0] for f in fails]) if fails else None; amin_fail = min([f[0] for f in fails]) if fails else None
print("local failure occurs for a in [%s, %s]" % (amin_fail, amax_fail))
# (2) global test
def rhs(V, y):
    ut = sgn*y[0]; rho = 2*V + 1/8 - ut**2/2
    return [sgn*(-(1.5 + ut)/rho)]
def phi(P, V1):
    if abs(V1 - P[0]) < 1e-14: return P[1]
    s = solve_ivp(rhs, (P[0], V1), [P[1]], rtol=1e-10, atol=1e-12)
    return s.y[0, -1] if s.success else np.nan
viol = 0; tested = 0; worst = []
for ip in range(NP):
    while True:
        V = np.exp(rng.uniform(np.log(1/32), np.log(0.5))); ut = rng.uniform(-1, 1)*(np.sqrt(8*V) - 0.5)
        if admissible(V, ut): break
    a = 1/V; c = a*(ut + 0.5)
    w, v = special(np.array([a]), np.array([c])); zeta = (w if which == 'w' else v)[0]
    # ellipse in (a', c'): f_Q(zeta) = l(zeta) + p(zeta)(a' zeta + c') = Z with |Z| < 1
    pz = zeta*(zeta + 1); lz = 2*zeta + 1
    for iq in range(NQ):
        r = np.sqrt(rng.uniform()); th = rng.uniform(0, 2*np.pi); Z = r*np.exp(1j*th)
        y = (Z - lz)/pz                     # = a' zeta + c'
        a2 = y.imag/zeta.imag; c2 = y.real - a2*zeta.real
        if not (0 < a2 < 32): continue
        V2 = 1/a2; ut2 = c2/a2 - 0.5
        if not admissible(V2, ut2) or V2 > 0.6: continue
        if abs(V2 - V) + abs(ut2 - ut) < 1e-9: continue
        tested += 1
        ph = phi((V, ut), V2)
        high = (ut2 > ph) if which == 'w' else (ut2 < ph)
        if not high:
            viol += 1
            if len(worst) < 8: worst.append((round(a, 4), round(ut, 4), round(a2, 4), round(ut2, 4), round(float(abs(Z)), 4)))
print("global test (%s): %d covering pairs tested, %d violations; examples (a_P, ut_P, a_Q, ut_Q, |f_Q(zeta_P)|):" % (which, tested, viol), worst)
