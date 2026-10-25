# ak_modsys.py -- does the MODULUS-ONLY system have real non-diagonal solutions?
#   R^-_j = sum_{i != j} log|f_i(w_j)| = 0,  R^+_j = sum_{i != j} log|f_i(v_j)| = 0   (2k eqs, 2k unknowns)
# Also (mode 'arg') the ARGUMENT-only system  Im sum log f_i(zeta) = 0 mod 2pi via sin/cos of total phase.
# usage: python3 ak_modsys.py k mode nstarts seed [lin]
import sys, numpy as np
from scipy.optimize import least_squares
from ak_core import ac_from, special, fval, logmat
k = int(sys.argv[1]); mode = sys.argv[2]; N = int(sys.argv[3]); seed = int(sys.argv[4])
lin = int(sys.argv[5]) if len(sys.argv) > 5 else 0
rng = np.random.default_rng(seed)
m = k - lin
def unpack(x):
    return ac_from(x[:m], x[m:])
def resid(x):
    a, c = unpack(x)
    Fw, Fv = logmat(a, c)
    if lin:
        w, v = special(a, c)
        lw = np.log(np.abs(2*w + 1)); lv = np.log(np.abs(2*v + 1))
        aw = np.angle(2*w + 1); av = np.angle(2*v + 1)
    else:
        lw = lv = aw = av = 0.0
    if mode == 'mod':
        return np.concatenate([np.log(np.abs(Fw)).sum(0) + lw, np.log(np.abs(Fv)).sum(0) + lv])
    else:  # arg: total phase must be 0 mod 2pi; use sin(phase/2) (vanishes iff phase in 2pi Z)
        pw = np.angle(Fw).sum(0) + aw; pv = np.angle(-Fv).sum(0) + av   # -Fv: f_i(v_j)/(-1)... product = 1 means prod_{i!=j} f_i(v_j) = 1
        pv = np.angle(Fv).sum(0) - np.angle(Fv[np.arange(m), np.arange(m)]) + av
        return np.concatenate([np.sin(pw/2), np.sin(pv/2)])
defl = (len(sys.argv) > 6 and sys.argv[6] == 'defl')
def P_of(x):
    a, c = unpack(x)
    return np.stack([x[:m], x[m:]], 1)
def resid_d(x):
    r = resid(x)
    if not defl: return r
    P = P_of(x)
    D2 = ((P[:, None, :] - P[None, :, :])**2).sum(-1)
    iu = np.triu_indices(m, 1)
    fac = np.prod(1.0 + 1.0/D2[iu])
    return r*fac
found = []
for it in range(N):
    x0 = np.concatenate([rng.normal(0, 2.5, m), rng.normal(0, 1.2, m)])
    try:
        r = least_squares(resid_d, x0, method='lm', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=4000)
    except Exception as e:
        continue
    x = r.x; res = np.max(np.abs(resid(x)))
    a, c = unpack(x)
    P = np.stack([np.log(a), c/np.sqrt(a)], 1)
    D = np.sqrt(((P[:, None, :] - P[None, :, :])**2).sum(-1)) + np.eye(m)*1e9
    dmin = D.min()
    if res < 1e-11 and dmin > 1e-3:
        found.append((res, dmin, a.copy(), c.copy()))
        print("SOL it=%d res=%.2e dmin=%.3g" % (it, res, dmin), "a=", np.array2string(a, precision=6), "c=", np.array2string(c, precision=6), flush=True)
    elif it < 5 or it % 50 == 0:
        print("it=%d res=%.2e dmin=%.3g" % (it, res, dmin), flush=True)
print("k=%d lin=%d mode=%s: %d genuine solutions in %d starts" % (k, lin, mode, len(found), N))
