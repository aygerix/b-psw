# audC_eab.py -- numerical sanity check of the auditor's lemmas E-A, E-B and the corollary (collinear cubics + 2t+1).
# E-A: i (smaller a) in cl N  => EXT(i,w) or EXT(i,v).
# E-B: i in W, a_j >= a_i, E^v_i(Q_j) < 0  => j in W and c_j < c_i.   (mirror: i in V ...)
# Lemma R: a > 2 in D => R < 0.   H > 0, K > 0 on W.
import sys, numpy as np
rng = np.random.default_rng(int(sys.argv[1]))
def fv(a, c, t): return 2*t + 1 + t*(t+1)*(a*t + c)
def up(A, B, C):
    r = np.roots([A, B, C]); return r[np.argmax(r.imag)]
def inD(a, c): return 0 < a < 32 and c*c < 8*a and (a-c)**2 < 8*a
def region(a, c): return "NWVB"[(a - c > 4) + 2*(c > 4)]
def samp():
    while True:
        a = 32*rng.random()**rng.choice([1, 0.5, 2]); lo, hi = a - np.sqrt(8*a), np.sqrt(8*a); c = lo + (hi - lo)*rng.random()
        if inD(a, c): return (a, c)
fa = fb = fr = fh = 0; na = nb = nr = 0
for _ in range(int(sys.argv[2])):
    Q1, Q2 = samp(), samp()
    if rng.random() < 0.3: s = rng.choice([1e-3, 1e-2, 1e-1]); Q2 = (Q1[0] + s*(Q2[0]-Q1[0]), Q1[1] + s*(Q2[1]-Q1[1]))
    Qi, Qj = (Q1, Q2) if Q1[0] <= Q2[0] else (Q2, Q1)
    a, c = Qi; w = up(a, a + c, c + 2); v = up(a, c, 2)
    Ew = fv(Qj[0], Qj[1], w).real - 1; Ev = (-fv(Qj[0], Qj[1], v)).real - 1
    al = abs(2*w+1) - 1; be = abs(2*v+1) - 1
    if c <= 4 and a - c <= 4:
        na += 1
        if not ((Ew >= 0 and al >= -1e-15) or (Ev >= 0 and be >= -1e-15)): fa += 1
    if region(*Qi) == "W" and Ev < 0:
        nb += 1
        if not (region(*Qj) == "W" and Qj[1] < Qi[1]): fb += 1
    if a > 2:
        nr += 1; R = -2/a - (c/a)*(1 - c/a)/2
        if R >= 0: fr += 1
    if region(a, c) == "W":
        H = 4*a*a - 6*a*c - c*c*(a - c); K = 4*a*a + a*c*(a - c)
        if H <= 0 or K <= 0: fh += 1
print("E-A: %d pairs, failures %d;  E-B: %d pairs, failures %d;  Lemma R: %d points, failures %d;  H,K>0 on W failures %d" % (na, fa, nb, fb, nr, fr, fh))
