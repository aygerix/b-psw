# audC_c23.py -- numerical sanity check of the AUDITOR's proofs of C2 and C3 (new code, mpmath 40 digits).
# C2 proof: P = (V, u) chart; line direction d with d_V >= 0; i = member maximising the coordinate along d;
#   A = <n^w_i, d>; if A <= 0 use w_i, else v_i.  Claim: all |f_l(z)| > 1 at the chosen point.
# C3 proof (ray c = u0 a, optional linear member 2t+1):
#   u0 <= -1: i = argmax a, point w_i;  u0 > -1, no linear: i = argmin a, point w_i;
#   u0 > -1 with linear (>= 2 cubics): alpha_i = 8V_i + 2u0 - 2 at i = argmin a: if alpha_i >= 0 use w_i (all >= 1, cubics > 1)
#   else use w_j, j = argmax a (all < 1).
# Also: the line ak's strict claim "at argmax V some point has Re(factor) > 1 for ALL others" on VERTICAL lines with R_i = 0.
import sys, mpmath as mp
mp.mp.dps = 40
import random
rnd = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 7)
def f(a, c, t): return 2*t + 1 + t*(t+1)*(a*t + c)
def uroot(A, B, C):
    d = mp.sqrt(B*B - 4*A*C); r1 = (-B + d)/(2*A); r2 = (-B - d)/(2*A)
    return r1 if mp.im(r1) > 0 else r2
def wv(a, c): return uroot(a, a + c, c + 2), uroot(a, c, mp.mpf(2))
def inD(a, c): return a > 0 and c*c < 8*a and (a - c)**2 < 8*a
def rand_D():
    while True:
        a = mp.mpf(rnd.uniform(0.001, 32)) if rnd.random() < 0.5 else mp.e**mp.mpf(rnd.uniform(-7, mp.log(32)))
        lo, hi = a - mp.sqrt(8*a), mp.sqrt(8*a)
        c = lo + (hi - lo)*mp.mpf(rnd.random())
        if inD(a, c): return a, c
fails2 = 0; n2 = 0; minmod2 = mp.inf
for trial in range(int(sys.argv[2]) if len(sys.argv) > 2 else 3000):
    a0, c0 = rand_D(); kind = rnd.random()
    m = rnd.randint(2, 9)
    if kind < 0.15:   # vertical line (all a equal) -- includes R = 0 endpoints sometimes
        if rnd.random() < 0.5 and a0 < 2:
            V = 1/a0; uR = (1 + mp.sqrt(1 + 16*V))/2   # R = 0 point at u = uR
            cs = [uR*a0] + [uR*a0 - mp.mpf(rnd.uniform(0, 3)) for _ in range(m-1)]
        else:
            cs = [c0 + mp.mpf(rnd.gauss(0, 1)) for _ in range(m)]
        pts = [(a0, c) for c in cs]
    else:
        th = rnd.uniform(0, mp.pi); sc = rnd.choice([0.01, 0.3, 3.0])
        pts = []
        for _ in range(m):
            s = mp.mpf(rnd.gauss(0, 1))*sc; pts.append((a0 + s*mp.cos(th), c0 + s*mp.sin(th)))
    pts = [q for q in pts if inD(*q)]
    pts = list({(mp.nstr(a, 30), mp.nstr(c, 30)): (a, c) for a, c in pts}.values())
    if len(pts) < 2: continue
    n2 += 1
    P = [(1/a, c/a) for a, c in pts]
    # direction d of the line in the P chart, with d_V >= 0
    dV, du = P[1][0] - P[0][0], P[1][1] - P[0][1]
    if dV < 0 or (dV == 0 and du < 0): dV, du = -dV, -du
    nrm = mp.sqrt(dV*dV + du*du); dV, du = dV/nrm, du/nrm
    i = max(range(len(P)), key=lambda l: P[l][0]*dV + P[l][1]*du)
    V, u = P[i]; R = -2*V - u*(1 - u)/2
    A = -(1 + u)*dV + R*du
    w, v = wv(*pts[i]); z, sgn = (w, 1) if A <= 0 else (v, -1)
    mods = [abs(f(a, c, z)) for l, (a, c) in enumerate(pts) if l != i]
    mm = min(mods); minmod2 = min(minmod2, mm)
    if not mm > 1: fails2 += 1; print("C2 FAIL", pts, i, A, mods)
print("C2 proof check: %d collinear configs (15%% vertical), failures %d, min chosen |factor| = %s" % (n2, fails2, mp.nstr(minmod2, 12)))
# C3
fails3 = 0; n3 = 0
for trial in range(int(sys.argv[2]) if len(sys.argv) > 2 else 3000):
    u0 = mp.mpf(rnd.uniform(-3, 4)); m = rnd.randint(1, 8); lin = rnd.random() < 0.7
    pts = []
    for _ in range(40):
        a = mp.e**mp.mpf(rnd.uniform(-7, mp.log(32)))
        if inD(a, u0*a): pts.append((a, u0*a))
        if len(pts) >= m: break
    if len(pts) + lin < 2 or len(pts) < 1: continue
    n3 += 1
    allm = pts + ([(mp.mpf(0), mp.mpf(0))] if lin else [])
    ia = min(range(len(pts)), key=lambda l: pts[l][0]); ja = max(range(len(pts)), key=lambda l: pts[l][0])
    if len(pts) == 1:   # {cubic, 2t+1}: product = 2w+1 != 1 since w non-real
        w, _ = wv(*pts[0]); ok = abs(mp.im(2*w + 1)) > 0
    elif u0 <= -1:
        w, _ = wv(*pts[ja]); ok = all(abs(f(a, c, w)) > 1 for l, (a, c) in enumerate(allm) if l != ja)
    elif not lin:
        w, _ = wv(*pts[ia]); ok = all(abs(f(a, c, w)) > 1 for l, (a, c) in enumerate(allm) if l != ia)
    else:
        a, c = pts[ia]; alpha = 8/a + 2*u0 - 2
        if alpha >= 0:
            w, _ = wv(*pts[ia]); vals = [abs(f(aa, cc, w)) for l, (aa, cc) in enumerate(allm) if l != ia]
            ok = all(x >= 1 for x in vals) and max(vals) > 1
        else:
            w, _ = wv(*pts[ja]); vals = [abs(f(aa, cc, w)) for l, (aa, cc) in enumerate(allm) if l != ja]
            ok = all(x < 1 for x in vals)
    if not ok: fails3 += 1; print("C3 FAIL", u0, pts, lin)
print("C3 proof check: %d ray configs, failures %d" % (n3, fails3))
# the line's strict claim on a vertical line with R = 0 at both endpoints: a = 1, u = (1 +- sqrt 17)/2, plus a middle member
a = mp.mpf(1); cs = [(1 + mp.sqrt(17))/2, (1 - mp.sqrt(17))/2, mp.mpf(1)/2]
pts = [(a, c) for c in cs]; print("vertical example in D:", [inD(*q) for q in pts])
for i in range(3):
    w, v = wv(*pts[i]); V, u = 1/a, cs[i]/a; R = -2*V - u*(1-u)/2
    rw = [mp.nstr(mp.re(f(a, c, w)), 15) for l, (aa, c) in enumerate(pts) if l != i]
    rv = [mp.nstr(mp.re(-f(a, c, v)), 15) for l, (aa, c) in enumerate(pts) if l != i]
    mw = [mp.nstr(abs(f(a, c, w)), 15) for l, (aa, c) in enumerate(pts) if l != i]
    print("member %d: R_i = %s  Re f_l(w_i) = %s  Re(-f_l(v_i)) = %s  |f_l(w_i)| = %s" % (i, mp.nstr(R, 5), rw, rv, mw))
