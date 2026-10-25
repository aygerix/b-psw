"""Validate the E2 exact criterion against src/bpsw.py on every odd squarefree composite n <= N0.
Fermat half:  strong_fermat2(n)  <=>  all v2(ord_p 2) equal (= a)  and  odd(ord_p) | n-1 for all p | n.
Lucas half:   strong_lucas(n)    <=>  [D = Selfridge D of n, gcd(n, Q) = 1]  all v2(z_D(p)) equal (= b)  and
              odd(z_D(p)) | n+1 for all p | n.   (The 2-adic conditions a <= v2(n-1), b <= v2(n+1) are NOT imposed:
              the claim says they are automatic; any n where they would matter shows up as a mismatch.)
usage: python3 aud_criterion.py N0 PATH_TO_SRC
"""
import sys, math
N0 = int(sys.argv[1]); sys.path.insert(0, sys.argv[2])
from bpsw import strong_fermat2, strong_lucas, selfridge_D
spf = list(range(N0 + 1))
for i in range(2, math.isqrt(N0) + 1):
    if spf[i] == i:
        for j in range(i * i, N0 + 1, i):
            if spf[j] == j: spf[j] = i
def fac(m):
    r = []
    while m > 1:
        q = spf[m]; r.append(q)
        while m % q == 0: m //= q
    return r
def order(g, p, N):          # order of g mod p dividing N
    o = N
    for q in fac(N):
        while o % q == 0 and pow(g, o // q, p) == 1: o //= q
    return o
def lucU(k, Q, p):
    U, V, Qk = 0, 2, 1; D = 1 - 4 * Q; inv2 = (p + 1) // 2
    for bit in bin(k)[2:]:
        U, V = U * V % p, (V * V - 2 * Qk) % p; Qk = Qk * Qk % p
        if bit == '1':
            U, V = (U + V) * inv2 % p, (D * U + V) * inv2 % p; Qk = Qk * Q % p
    return U
zc = {}
def zD(D, p):
    key = (D, p)
    if key in zc: return zc[key]
    Q = (1 - D) // 4
    if Q % p == 0: zc[key] = None; return None
    kr = pow(D % p, (p - 1) // 2, p); kr = -1 if kr == p - 1 else kr
    N = p - kr
    z = N
    for q in fac(N):
        while z % q == 0 and lucU(z // q, Q, p) == 0: z //= q
    assert lucU(z, Q, p) == 0
    zc[key] = z; return z
oc = {}
v2 = lambda x: (x & -x).bit_length() - 1
odd = lambda x: x >> v2(x)
mis = 0; nF = nL = nB = 0; tested = 0; need2a = need2b = 0
for n in range(9, N0 + 1, 2):
    ps = fac(n)
    if len(ps) < 2 or math.prod(ps) != n: continue      # squarefree composites only
    tested += 1
    os_ = []
    for p in ps:
        if p not in oc: oc[p] = order(2, p, p - 1)
        os_.append(oc[p])
    cF = len({v2(o) for o in os_}) == 1 and all((n - 1) % odd(o) == 0 for o in os_)
    if cF and not (n - 1) % (1 << v2(os_[0])) == 0: need2a += 1
    tF = strong_fermat2(n)
    D = selfridge_D(n); Q = (1 - D) // 4
    zs = [zD(D, p) for p in ps]
    if any(z is None for z in zs): cL = False
    else:
        cL = len({v2(z) for z in zs}) == 1 and all((n + 1) % odd(z) == 0 for z in zs)
        if cL and not (n + 1) % (1 << v2(zs[0])) == 0: need2b += 1
    tL = strong_lucas(n)
    nF += tF; nL += tL; nB += tF and tL
    if cF != tF or cL != tL:
        mis += 1
        if mis <= 20: print("MISMATCH", n, ps, "F", cF, tF, "L", cL, tL, "D", D)
print(f"N0 = {N0}: squarefree odd composites tested {tested}; strong psp(2) {nF}, strong Lucas psp {nL}, both {nB}; "
      f"criterion mismatches {mis}; cases where the 2-adic side condition would fail: Fermat {need2a}, Lucas {need2b}")
