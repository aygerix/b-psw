# Half-plane criterion C4 for real cubic Williams families f = 2t+1 + t(t+1)(a t + c), a > 0.
# (1) check identity  Re f_l(w_i) - 1 = a_l <n^w_i, P_l - P_i>,  Re(-f_l(v_i)) - 1 = -a_l <n^v_i, P_l - P_i>
#     with P = (V, ut) = (1/a, c/a - 1/2), u = c/a, n^w = (-u-1, R), n^v = (2-u, R), R = -2V - u(1-u)/2;
#     w_i = upper root of g_i = 2 + (t+1)(a_i t + c_i), v_i = upper root of h_i = 2 + t(a_i t + c_i).
# (2) random configurations: fraction with NO member triggering C4.
import numpy as np, sys
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
def f(a, c, t): return 2*t + 1 + t*(t+1)*(a*t + c)
def uroot(coef):
    r = np.roots(coef); return r[np.argmax(r.imag)]
def data(a, c):
    V, u = 1/a, c/a
    w = uroot([a, a + c, c + 2]); v = uroot([a, c, 2])
    R = -2*V - u*(1-u)/2
    return np.array([V, u - 0.5]), np.array([-u-1, R]), np.array([2-u, R]), w, v
def sample():
    while True:
        V = rng.exponential(0.5) + 1e-3; u = rng.normal(0.5, 1.5)
        if u*u < 8*V and (1-u)**2 < 8*V: return 1/V, u/V
err = 0
for _ in range(2000):
    (ai, ci), (al, cl) = sample(), sample()
    Pi, nw, nv, w, v = data(ai, ci); Pl = data(al, cl)[0]
    err = max(err, abs(f(al, cl, w).real - 1 - al*nw @ (Pl - Pi)), abs((-f(al, cl, v)).real - 1 + al*nv @ (Pl - Pi)))
print("identity max err", err)
for k in (3, 5, 7, 9):
    bad = 0; T = 20000
    for _ in range(T):
        pts = [sample() for _ in range(k)]
        D = [data(a, c) for a, c in pts]
        trig = False
        for i in range(k):
            dw = [D[i][1] @ (D[l][0] - D[i][0]) for l in range(k) if l != i]
            dv = [D[i][2] @ (D[l][0] - D[i][0]) for l in range(k) if l != i]
            if min(dw) > 0 or max(dv) < 0: trig = True; break
        bad += not trig
    print("k", k, "configs with no C4 trigger:", bad, "/", T)
