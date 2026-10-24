"""W5c: single-row character bias. For each row q^e and each prime l | phi(q^e) (l <= 50), compare the number of
off-support pool primes whose discrete log is == 0 mod l (p is an l-th power residue mod q^e) with the matched
model's expectation. Off support means q does not divide p -+ 1 on the row's own side; such p are uniform among
units x != +-1 (mod q) under the model.
Usage: python3 rowbias.py system.txt
"""
import sys, math
from structtest import read_system, vq

def main():
    ps, coords = read_system(sys.argv[1])
    out, chi2, dof = [], 0.0, 0
    for q, e, ph, g, side, t, row in coords:
        m = q ** e; s = 1 if side == 'F' else -1
        off = [i for i, p in enumerate(ps) if vq(p - s, q) == 0]
        if len(off) < 200: continue
        # model: x uniform among units mod q^e with x mod q not in {1, -1}; count l-th powers exactly
        units = [x for x in range(1, m) if x % q and (x - 1) % q and (x + 1) % q]
        lg = {}; x = 1
        for k in range(ph): lg[x] = k; x = x * g % m
        for l in [l for l in range(2, 51) if ph % l == 0 and all(l % d for d in range(2, int(l ** .5) + 1))]:
            pexp = sum(1 for x in units if lg[x] % l == 0) / len(units)
            obs = sum(1 for i in off if row[i] % l == 0)
            exp = pexp * len(off); sd = math.sqrt(len(off) * pexp * (1 - pexp))
            if sd == 0: continue
            z = (obs - exp) / sd
            chi2 += z * z; dof += 1
            out.append((abs(z), z, q, e, side, l, obs, exp, len(off)))
    out.sort(reverse=True)
    print(f"rows x primes tested: {dof}; sum z^2 = {chi2:.1f} (expect {dof} +- {math.sqrt(2 * dof):.1f})")
    for a, z, q, e, side, l, obs, exp, n in out[:15]:
        print(f"  q^e = {q}^{e} ({side}) l = {l}: {obs} of {n} off-support primes are l-th powers, expected {exp:.1f} (z {z:+.2f})")

if __name__ == "__main__":
    main()
