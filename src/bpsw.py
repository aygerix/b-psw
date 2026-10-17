"""Exact Baillie-PSW check following the problem statement verbatim.
Strong Fermat base 2 + strong Lucas (Selfridge Method A, P=1)."""
import sys, json
from math import isqrt

def jacobi(a, n):
    assert n > 0 and n % 2 == 1
    a %= n
    r = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                r = -r
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            r = -r
        a %= n
    return r if n == 1 else 0

def strong_fermat2(n):
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2; s += 1
    x = pow(2, d, n)
    if x == 1 or x == n - 1:
        return True
    for _ in range(s - 1):
        x = x * x % n
        if x == n - 1:
            return True
    return False

def selfridge_D(n):
    if isqrt(n) ** 2 == n:
        return None          # no D exists for perfect squares
    D = 5
    while True:
        if jacobi(D, n) == -1:
            return D
        D = -D - 2 if D > 0 else -D + 2

def lucas_UV(k, P, Q, n):
    """U_k, V_k, Q^k mod n via binary method on (U,V) with U_{2m}=U_m V_m, V_{2m}=V_m^2-2Q^m."""
    D = P * P - 4 * Q
    U, V, Qk = 0, 2, 1          # index 0
    inv2 = (n + 1) // 2
    for bit in bin(k)[2:]:
        U, V = U * V % n, (V * V - 2 * Qk) % n
        Qk = Qk * Qk % n
        if bit == '1':
            U, V = (P * U + V) * inv2 % n, (D * U + P * V) * inv2 % n
            Qk = Qk * Q % n
    return U, V, Qk

def strong_lucas(n):
    D = selfridge_D(n)
    if D is None:
        return False
    P, Q = 1, (1 - D) // 4
    e, t = n + 1, 0
    while e % 2 == 0:
        e //= 2; t += 1
    U, V, Qk = lucas_UV(e, P, Q, n)
    if U == 0 or V == 0:
        return True
    for _ in range(t - 1):
        V = (V * V - 2 * Qk) % n
        Qk = Qk * Qk % n
        if V == 0:
            return True
    return False

def bpsw_pass(n):
    return n > 1 and n % 2 == 1 and strong_fermat2(n) and strong_lucas(n)

def check_submission(path):
    o = json.load(open(path))
    f1, f2 = int(o["factor1"]), int(o["factor2"])
    assert f1 > 1 and f2 > 1
    n = f1 * f2
    return n % 2 == 1 and strong_fermat2(n) and strong_lucas(n)

if __name__ == "__main__":
    print(check_submission(sys.argv[1]))
