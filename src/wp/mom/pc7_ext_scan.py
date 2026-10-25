# extremal-member sign scan for (2,2) quartic pencils
import numpy as np, sys
def shape(r1, r2):
    e1, e2 = r1 + r2, r1 * r2
    Dl = (1 - r1**2) * (1 - r2**2)
    lam, mu = 2 * e1 / Dl, 2 * (1 - e1**2 + e2) / Dl
    return e1, e2, lam, mu
def info(r1, r2, c):
    e1, e2, lam, mu = shape(r1, r2)
    # K- = c H + L = c t^2 + (lam - c e1) t + (c e2 + mu)
    Km = [c, lam - c * e1, c * e2 + mu]
    Kp = [c, lam, mu + lam * e1 - c]
    dm = Km[1]**2 - 4 * Km[0] * Km[2]; dp = Kp[1]**2 - 4 * Kp[0] * Kp[2]
    th_m = np.roots(Km); th_m = th_m[np.argmax(th_m.imag)]
    th_p = np.roots(Kp); th_p = th_p[np.argmax(th_p.imag)]
    G = lambda t: t * t - 1; H = lambda t: t * t - e1 * t + e2
    L = lambda t: lam * t + mu; Lp = lambda t: lam * t + mu + lam * e1
    Pm = G(th_m) * H(th_m); Pp = G(th_p) * H(th_p)
    return dm, dp, Pm, Pp
rs = [(-0.9, 0.5), (-0.5, 0.3), (0.1, 0.8), (-0.95, -0.9), (0.2, 0.3), (-0.3, 0.9), (0.5, 0.95), (-0.8, 0.8)]
for r1, r2 in rs:
    cs = np.exp(np.linspace(-8, 8, 40001))
    rows = []
    for c in cs:
        dm, dp, Pm, Pp = info(r1, r2, c)
        if dm < 0 and dp < 0:
            rows.append((c, np.sign(Pm.real), np.sign(Pp.real)))
    if not rows:
        print(r1, r2, 'no admissible c'); continue
    rows = np.array(rows)
    # region A (bottom ok): Re P- < 0 and Re P+ > 0 ; region B (top ok): Re P- > 0 and Re P+ < 0
    A = rows[(rows[:, 1] < 0) & (rows[:, 2] > 0), 0]
    B = rows[(rows[:, 1] > 0) & (rows[:, 2] < 0), 0]
    # transitions
    segs = []
    prev = None
    for c, a, b in rows:
        key = (a, b)
        if key != prev: segs.append([c, c, key]); prev = key
        else: segs[-1][1] = c
    print(f"r=({r1},{r2}) adm c in [{rows[0,0]:.4g},{rows[-1,0]:.4g}]  A:[{A.min() if len(A) else None},{A.max() if len(A) else None}] B:[{B.min() if len(B) else None},{B.max() if len(B) else None}]")
    for s in segs: print("   ", f"{s[0]:.5g}..{s[1]:.5g}", s[2])
