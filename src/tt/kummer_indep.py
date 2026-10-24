"""Two-torus test (T3): are the index events of 2 (in F_p^*) and of gamma = phi/phibar (in the norm-1 torus T_p)
independent for primes p inert in Q(sqrt 5)?

For each inert p <= X and each odd prime q <= QMAX:
  F-event q:  q | p - 1 and q | (p - 1)/ord_p(2)   (2 is a q-th power residue mod p)
  L-event q:  q | p + 1 and q | (p + 1)/z(p)       (U_{(p+1)/q} == 0 mod p, i.e. gamma is a q-th power in T_p)
Kummer theory predicts P(F-event q | q | p-1) = 1/q and P(L-event q | q | p+1) = 1/q, independent across all
(q1, q2) pairs. A positive correlation would be a 'two-torus' lead: index events of 2 and gamma occurring together.
Also tabulates the 2-adic data (p mod 8, v2 z(p)) against the odd F-events.
Usage: python3 kummer_indep.py X [QMAX]
"""
import sys, math
from collections import defaultdict

def primes_upto(n):
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = bytearray(len(range(i * i, n + 1, i)))
    return [i for i in range(n + 1) if s[i]]

def fib(k, p):
    """(F_k, F_{k+1}) mod p by fast doubling."""
    a, b = 0, 1
    for bit in bin(k)[2:]:
        c = a * ((2 * b - a) % p) % p
        d = (a * a + b * b) % p
        if bit == '1': a, b = d, (c + d) % p
        else: a, b = c, d
    return a, b

def main():
    X = int(sys.argv[1]); QMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 23
    qs = [q for q in primes_upto(QMAX) if q > 2 and q != 5]
    nF = defaultdict(int); hF = defaultdict(int); nL = defaultdict(int); hL = defaultdict(int)
    nFL = defaultdict(int); hFL = defaultdict(lambda: [0, 0, 0, 0])
    v2tab = defaultdict(int)
    ninert = 0
    for p in primes_upto(X):
        if p < 7 or p % 5 not in (2, 3): continue
        ninert += 1
        Fe, Le = {}, {}
        for q in qs:
            if (p - 1) % q == 0:
                Fe[q] = pow(2, (p - 1) // q, p) == 1
                nF[q] += 1; hF[q] += Fe[q]
            elif (p + 1) % q == 0:
                Le[q] = fib((p + 1) // q, p)[0] == 0
                nL[q] += 1; hL[q] += Le[q]
        for q1, e1 in Fe.items():
            for q2, e2 in Le.items():
                nFL[(q1, q2)] += 1; hFL[(q1, q2)][2 * e1 + e2] += 1
        # 2-adic: v2 of the gamma-index restricted to the 2-part of p+1
        m = p + 1; t = (m & -m).bit_length() - 1
        j = 0
        while j < t and fib(m >> (j + 1), p)[0] == 0: j += 1   # j = v2 of the index of gamma
        key = (p % 8, j, tuple(sorted(q for q, e in Fe.items() if e))[:1])
        v2tab[(p % 8, j, bool(any(Fe.values())))] += 1
    print(f"inert primes 7 <= p <= {X}: {ninert}")
    print("single events: q, P(F-event | q|p-1) * q, P(L-event | q|p+1) * q   (Kummer: both 1)")
    for q in qs:
        print(f"  q={q:3d}: F {hF[q] / nF[q] * q:.4f} (n={nF[q]}), L {hL[q] / nL[q] * q:.4f} (n={nL[q]})")
    print("pairs (q1 F-side, q2 L-side): observed joint events vs independence; chi2 of the 2x2 table")
    tot_chi = 0.0; dof = 0; worst = []
    for (q1, q2), n in sorted(nFL.items()):
        c = hFL[(q1, q2)]
        pf = sum(c[2:]) / n; pl = (c[1] + c[3]) / n
        exp = [n * (1 - pf) * (1 - pl), n * (1 - pf) * pl, n * pf * (1 - pl), n * pf * pl]
        if min(exp) < 5: continue
        chi = sum((o - e) ** 2 / e for o, e in zip(c, exp)); tot_chi += chi; dof += 1
        worst.append((chi, q1, q2, c[3], exp[3], n))
    worst.sort(reverse=True)
    print(f"  {dof} pairs with expected cells >= 5: sum chi2 = {tot_chi:.1f} (expect {dof} +- {math.sqrt(2 * dof):.1f})")
    for chi, q1, q2, o, e, n in worst[:6]:
        print(f"    ({q1},{q2}): joint {o} vs {e:.1f} expected (n={n}), chi2 {chi:.2f}")
    print("2-adic: (p mod 8, v2 of gamma-index, some odd F-event) -> count")
    for k in sorted(v2tab): print("  ", k, v2tab[k])

if __name__ == "__main__":
    main()
