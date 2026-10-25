# audC_red.py -- reduced form of the collinear + 2t+1 problem.
# Line L through cubic points meets the c-axis at Chat = (0, c0) (L not vertical, c0 != 0 since O not on L).
# For X in D and a member Y on L with a_Y < a_X: f_Y(z) - 1 = delta (f_Chat(z) - 1), delta = 1 - a_Y/a_X in (0,1);
# for a_Y > a_X, delta < 0.  f_Chat = 2t+1 + c0 t(t+1) (quadratic), f_O = 2t+1.
# Single-point data at X (closed forms, checked here against direct evaluation):
#   Rw = Re f_Chat(w_X) - 1 = -(1+u) + c0 R,  Rv = Re(-f_Chat(v_X)) - 1 = (u-2) - c0 R
#   Mw = |f_Chat(w_X)|^2 - 1 = (8V+2u-2) + 2 c0 u((1-u)/2 - 2V) + c0^2 2V(u+2V)
#   Mv = |f_Chat(v_X)|^2 - 1 = (8V-2u) + 2 c0 (1-u)(2V - u/2) + c0^2 2V(1-u+2V)
#   alpha = |2w+1|^2 - 1 = 8V+2u-2, beta = |2v+1|^2 - 1 = 8V-2u.
# T(X) "good as top (min-a) endpoint" := (Rw <= 0 & alpha >= 0) | (Rv <= 0 & beta >= 0)
# B(X) "good as bottom (max-a) endpoint" := (Mw <= 0 & alpha <= 0) | (Mv <= 0 & beta <= 0) | (Rw >= 0 & alpha >= 0) | (Rv >= 0 & beta >= 0)
# A line through Chat meets D in a segment; along it a increases monotonically. Robust coverage FAILS on that line iff
# some X1 with ~T precedes (smaller a) some X2 with ~B.  We scan lines through Chat (slopes mu) for each c0.
import sys, numpy as np
def fv(a, c, t): return 2*t + 1 + t*(t+1)*(a*t + c)
def up(A, B, C):
    d = np.sqrt(B*B - 4*A*C + 0j); r = (-B + d)/(2*A); return np.where(r.imag >= 0, r, np.conj(r))
def data(a, c, c0):
    V, u = 1/a, c/a; R = -2*V - u*(1-u)/2
    Rw = -(1+u) + c0*R; Rv = (u-2) - c0*R
    Mw = (8*V+2*u-2) + 2*c0*u*((1-u)/2 - 2*V) + c0**2*2*V*(u+2*V)
    Mv = (8*V-2*u) + 2*c0*(1-u)*(2*V - u/2) + c0**2*2*V*(1-u+2*V)
    return Rw, Rv, Mw, Mv, 8*V+2*u-2, 8*V-2*u
# closed-form check
rng = np.random.default_rng(0); err = 0
for _ in range(2000):
    a = 32*rng.random(); c = a - np.sqrt(8*a) + 2*np.sqrt(8*a)*rng.random() - a*rng.random()*0
    if not (c*c < 8*a and (a-c)**2 < 8*a): continue
    c0 = rng.normal(0, 5); w = up(a, a+c, c+2); v = up(a, c, 2.0)
    Rw, Rv, Mw, Mv, al, be = data(a, c, c0)
    err = max(err, abs(fv(0, c0, w).real - 1 - Rw), abs((-fv(0, c0, v)).real - 1 - Rv), abs(abs(fv(0, c0, w))**2 - 1 - Mw), abs(abs(fv(0, c0, v))**2 - 1 - Mv), abs(abs(2*w+1)**2 - 1 - al), abs(abs(2*v+1)**2 - 1 - be))
print("closed-form max error:", err)
def T(d): Rw, Rv, Mw, Mv, al, be = d; return (Rw <= 0) & (al >= 0) | (Rv <= 0) & (be >= 0)
def B(d): Rw, Rv, Mw, Mv, al, be = d; return (Mw <= 0) & (al <= 0) | (Mv <= 0) & (be <= 0) | (Rw >= 0) & (al >= 0) | (Rv >= 0) & (be >= 0)
c0s = [float(x) for x in sys.argv[1].split(",")] if len(sys.argv) > 1 else [0.5, 1, 2, 3.9, 4, 4.1, 6, 10, 30, -0.5, -2, -5]
NS = int(sys.argv[2]) if len(sys.argv) > 2 else 4000; NA = 4000
for c0 in c0s:
    bad = 0; nlines = 0; nT = 0; nB = 0; npts = 0; ex = None
    # lines through (0, c0): c = c0 + mu a. D lies in 0 < a < 32; c in (a - sqrt(8a), sqrt(8a)) => mu range bounded per a.
    for mu in np.linspace(-6, 6, NS):
        a = np.linspace(1e-6, 32, NA)[1:-1]; c = c0 + mu*a
        m = (c*c < 8*a) & ((a-c)**2 < 8*a)
        if m.sum() < 2: continue
        a, c = a[m], c[m]; nlines += 1
        d = data(a, c, c0); t = T(d); b = B(d); npts += len(a); nT += (~t).sum(); nB += (~b).sum()
        # failure: exists idx1 < idx2 with ~t[idx1] and ~b[idx2]
        first_notT = np.argmax(~t) if (~t).any() else None
        last_notB = len(b) - 1 - np.argmax((~b)[::-1]) if (~b).any() else None
        if first_notT is not None and last_notB is not None and first_notT < last_notB:
            bad += 1
            if ex is None: ex = (mu, a[first_notT], c[first_notT], a[last_notB], c[last_notB])
    print("c0 = %6.3f: lines %d, points %d, notT %d, notB %d, lines with robust failure %d, example %s" % (c0, nlines, npts, nT, nB, bad, ex))
