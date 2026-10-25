#!/usr/bin/env python3
"""xn line: adaptive cover of an S-range by RIGOROUS Farkas certificates (xn_farkas9.certify), k = 9.

Works in R = sqrt(8 (S - (m-1)/32)) - 1/2 (the r-coordinate of the largest admissible member, rational endpoints).
For every interval tries the flag sets in order (simplest first), bisects on failure.
Usage: xn_cover9.py lin R0 R1 step [minstep [certfile]]      (R0 = 0 is S = m/32)
certfile: one gp vector per certified interval  [lin, Ra, Rb, y, [[kind, param, z], ...]]  (exact rationals)
Output lines 'OK lin R-interval S-interval flags message' or 'FAIL ...'.
"""
import sys
from fractions import Fraction as Fr
import xn_farkas9 as X

FLAGSETS = [["CS"], ["CS", "PHI"], ["CS", "PHI", "POS"], ["CS", "PHI", "POS", "N3L"], ["CS", "PHI", "POS", "N3L", "N3U"]]


def try_interval(lin, Ra, Rb, m):
    Sa, Sb = X.S_of_R(Ra, m), X.S_of_R(Rb, m)
    for fl in FLAGSETS:
        try:
            ok, msg = X.certify_R(lin, Ra, Rb, fl)
        except Exception as e:  # numerical trouble -> treat as failure
            ok, msg = False, f"exception {e}"
        if ok:
            return True, fl, msg
    return False, None, msg


CERTF = None


def main():
    global CERTF
    lin = int(sys.argv[1])
    if len(sys.argv) > 6:
        CERTF = open(sys.argv[6], "w")
    R0, R1, step = Fr(sys.argv[2]), Fr(sys.argv[3]), Fr(sys.argv[4])
    minstep = Fr(sys.argv[5]) if len(sys.argv) > 5 else Fr(1, 256)
    m = 9 - lin
    stack = []
    R = R0
    while R < R1:
        stack.append((R, min(R + step, R1)))
        R += step
    stack.reverse()
    while stack:
        Ra, Rb = stack.pop()
        ok, fl, msg = try_interval(lin, Ra, Rb, m)
        Sa, Sb = X.S_of_R(Ra, m), X.S_of_R(Rb, m)
        if ok and CERTF is not None:
            L = X.LAST_CERT
            zs = ", ".join(f'["{nm}", {pr}, {zz}]' for nm, pr, zz in L["z"])
            CERTF.write(f'[{L["lin"]}, {L["Ra"]}, {L["Rb"]}, [{", ".join(str(q) for q in L["y"])}], [{zs}]]\n')
            CERTF.flush()
        if ok:
            print(f"OK   lin={lin} R=[{Ra},{Rb}] S=[{float(Sa):.6f},{float(Sb):.6f}] flags={'+'.join(fl)} | {msg.split(': eps')[1] if ': eps' in msg else msg}", flush=True)
        elif Rb - Ra > minstep:
            Rm = (Ra + Rb) / 2
            stack.append((Rm, Rb))
            stack.append((Ra, Rm))
        else:
            print(f"FAIL lin={lin} R=[{Ra},{Rb}] S=[{float(Sa):.6f},{float(Sb):.6f}] | {msg}", flush=True)


if __name__ == "__main__":
    main()
