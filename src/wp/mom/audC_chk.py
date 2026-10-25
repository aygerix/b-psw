# audC_chk.py -- evaluate the 8 robust endpoint tests (as in audC_rob.py) at given endpoint pairs, 40 digits.
import sys, mpmath as mp
mp.mp.dps = 40
def f(a, c, t): return 2*t + 1 + t*(t+1)*(a*t + c)
def ur(A, B, C):
    d = mp.sqrt(B*B - 4*A*C); r = (-B + d)/(2*A); return r if mp.im(r) > 0 else (-B - d)/(2*A)
def tests(Qi, Qj):
    out = {}
    for tag, (A, B) in (("i", (Qi, Qj)), ("j", (Qj, Qi))):
        a, c = A; w = ur(a, a + c, c + 2); v = ur(a, c, mp.mpf(2))
        for zn, F in (("w", lambda Q: f(Q[0], Q[1], w)), ("v", lambda Q: -f(Q[0], Q[1], v))):
            Z = F(B) - 1; FO = F((0, 0))
            out["EXT"+tag+zn] = (mp.re(Z) >= 0 and abs(FO) >= 1)
            out["INT"+tag+zn] = (abs(1 + Z) <= 1 and abs(FO) <= 1)
            out["mod"+tag+zn] = (mp.nstr(abs(1+Z), 8), mp.nstr(abs(FO), 8), mp.nstr(mp.re(Z), 6))
    return out
pairs = eval(sys.argv[1])
for Qi, Qj in pairs:
    Qi = tuple(mp.mpf(x) for x in Qi); Qj = tuple(mp.mpf(x) for x in Qj)
    t = tests(Qi, Qj)
    print("Qi", [mp.nstr(x, 10) for x in Qi], "Qj", [mp.nstr(x, 10) for x in Qj], "passing:", [k for k, v in t.items() if v is True])
    print("   (|cubic factor|, |O factor|, Re Z):", {k: v for k, v in t.items() if k.startswith("mod")})
