"""Selector line, check of Lemma S1 and the exact Chebotarev densities of Kummer selectors (inert primes, D = 5).

Kummer selector K(theta, m), theta = 2^a gamma^b (gamma = -phi^2, the Lucas element alpha/beta for P = 1, Q = -1):
  pO_K splits completely in K(zeta_m, theta^(1/m))  <=>  p^2 == 1 (mod m) and theta^((p^2-1)/m) == 1 in F_{p^2}.
Lemma S1 (PROVED) predicts, for inert p: K(theta, m) holds iff for every q^e || m:
  q^e | p - 1 and v_q(i_2) >= e - v_q(a),   or   q^e | p + 1 and v_q(i_g) >= e - v_q(b),
with i_2 = (p-1)/ord_p 2, i_g = (p+1)/ord(gamma).  (Mixed generators never couple the two sides at one q.)
Exact density among inert primes: prod over q^e || m of (1/phi(q^e)) (q^-(e - v_q a)^+ + q^-(e - v_q b)^+).
The script tests S1 deterministically (every p) and compares the densities; it also evaluates the advisor's
correlated selector S = (A and B) or (not A and not B), A = {3 | i_2}, B = {7 | i_g}.
Usage: python3 sel_check.py X
"""
import sys
from sympy import primerange, factorint, n_order

def mul(x, y, p):          # elements a + b s of F_p[s]/(s^2 - 5)
    return ((x[0] * y[0] + 5 * x[1] * y[1]) % p, (x[0] * y[1] + x[1] * y[0]) % p)

def pw(x, e, p):
    r = (1, 0)
    while e:
        if e & 1: r = mul(r, x, p)
        x = mul(x, x, p); e >>= 1
    return r

def ordT(g, p, fac):       # order of g (norm 1) dividing p + 1, fac = factorint(p + 1)
    o = p + 1
    for q, e in fac.items():
        for _ in range(e):
            if pw(g, o // q, p) == (1, 0): o //= q
            else: break
    return o

def vq(x, q):
    v = 0
    while x % q == 0: x //= q; v += 1
    return v

def phi_pp(q, e): return (q - 1) * q ** (e - 1)

def main():
    X = int(sys.argv[1])
    sels = [(1, 0, 3), (0, 1, 3), (1, 1, 3), (1, 1, 21), (3, 1, 9), (1, 3, 9), (1, 2, 63), (2, 1, 7 * 11), (1, 1, 9 * 7),
            (3, 7, 21)]
    hit = {s: 0 for s in sels}; mism = {s: 0 for s in sels}
    nA = nB = nAB = nS = 0; ninert = 0
    for p in primerange(7, X):
        if p % 5 not in (2, 3): continue
        ninert += 1
        inv2 = (p + 1) // 2
        gam = ((-3 * inv2) % p, (-inv2) % p)        # -(3 + s)/2
        i2 = (p - 1) // n_order(2, p)
        ig = (p + 1) // ordT(gam, p, factorint(p + 1))
        A = (p - 1) % 3 == 0 and i2 % 3 == 0
        B = (p + 1) % 7 == 0 and ig % 7 == 0
        nA += A; nB += B; nAB += A and B
        if A == B: nS += 1
        for (a, b, m) in sels:
            if (p * p - 1) % m: direct = False
            else:
                th = mul(pw((2, 0), a, p), pw(gam, b, p), p)
                direct = pw(th, (p * p - 1) // m, p) == (1, 0)
            pred = True
            for q, e in factorint(m).items():
                if (p - 1) % q ** e == 0: ok = vq(i2, q) >= e - vq(a, q) if a else True
                elif (p + 1) % q ** e == 0: ok = vq(ig, q) >= e - vq(b, q) if b else True
                else: ok = False
                pred = pred and ok
            hit[(a, b, m)] += direct; mism[(a, b, m)] += direct != pred
    print(f"inert primes 7 <= p < {X}: {ninert}")
    print("selector (a, b, m): observed density * exact prediction^-1, mismatches with Lemma S1 (must be 0)")
    for (a, b, m) in sels:
        d = 1.0
        for q, e in factorint(m).items():
            ta = q ** -max(0, e - (vq(a, q) if a else 99)); tb = q ** -max(0, e - (vq(b, q) if b else 99))
            d *= (ta + tb) / phi_pp(q, e)
        print(f"  theta = 2^{a} gamma^{b}, m = {m:3d}: obs {hit[(a, b, m)] / ninert:.6f}  exact {d:.6f}  ratio "
              f"{hit[(a, b, m)] / ninert / d:.4f}  (count {hit[(a, b, m)]}), mismatches {mism[(a, b, m)]}")
    print(f"advisor selector: P(A) = {nA / ninert:.5f} (1/6 = {1 / 6:.5f}), P(B) = {nB / ninert:.5f} (1/42 = {1 / 42:.5f}), "
          f"P(A and B) = {nAB / ninert:.6f} (1/252 = {1 / 252:.6f}), P(S) = {nS / ninert:.4f} (206/252 = {206 / 252:.4f}), "
          f"P(A and B | S) = {nAB / nS:.5f} (1/206 = {1 / 206:.5f}); absolute density of A and B is unchanged by S")

if __name__ == "__main__":
    main()
