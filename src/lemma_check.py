"""Check Lemma 1 (exact per-prime characterization, squarefree n, D=5) against direct tests,
and the 2-adic facts for inert primes."""
import random
from sympy import isprime, factorint, primerange, n_order, divisors
from bpsw import jacobi, strong_fermat2, strong_lucas, selfridge_D

def v2(x):
    return (x & -x).bit_length() - 1

def odd(x):
    return x >> v2(x)

def fib_entry(p):
    """z(p) = order of gamma = alpha/beta: least k with F_k = 0 mod p (p != 2,5)."""
    m = p - jacobi(5, p)
    def F(k):  # fast doubling mod p
        a, b = 0, 1
        for bit in bin(k)[2:]:
            c = a * ((2 * b - a) % p) % p
            d = (a * a + b * b) % p
            a, b = (d, (c + d) % p) if bit == '1' else (c, d)
        return a
    for d in sorted(divisors(m)):
        if F(d) == 0:
            return d

def predicted(fs):
    n = 1
    for p in fs: n *= p
    if jacobi(5, n) != -1:
        return None  # D != 5: lemma not applicable
    s, t = v2(n - 1), v2(n + 1)
    j1 = {v2(n_order(2, p)) for p in fs}
    j2 = {v2(fib_entry(p)) for p in fs}
    if len(j1) != 1 or len(j2) != 1: return False
    if j1.pop() > s or j2.pop() > t: return False
    return all((n - 1) % odd(n_order(2, p)) == 0 and (n + 1) % odd(fib_entry(p)) == 0 for p in fs)

# (a) 2-adic facts for inert primes
bad = 0
for p in primerange(7, 200000):
    if jacobi(5, p) != -1: continue
    z = fib_entry(p)
    if p % 4 == 3 and v2(z) != v2(p + 1): bad += 1
    if p % 4 == 1 and v2(z) != 0: bad += 1
    if p % 8 in (3, 5) and v2(n_order(2, p)) != v2(p - 1): bad += 1
print("2-adic fact violations (inert p<2e5):", bad)

# (b) lemma vs direct test on random squarefree n with D=5, enriched with structured cases
random.seed(1)
P = list(primerange(7, 3000))
agree = disagree = passes = 0
for trial in range(40000):
    k = random.choice([2, 3, 3, 4, 5])
    fs = sorted(random.sample(P, k))
    pr = predicted(fs)
    if pr is None: continue
    n = 1
    for p in fs: n *= p
    direct = strong_fermat2(n) and strong_lucas(n)
    agree += (pr == direct); disagree += (pr != direct); passes += direct
print("random: agree", agree, "disagree", disagree, "passes", passes)

# (c) the half-tests separately on all strong-Lucas psp and spsp2 below 1e6 with D=5 (squarefree)
def half_check(n, fs, which):
    s, t = v2(n - 1), v2(n + 1)
    if which == 'F':
        j = {v2(n_order(2, p)) for p in fs}
        return len(j) == 1 and max(j) <= s and all((n - 1) % odd(n_order(2, p)) == 0 for p in fs)
    j = {v2(fib_entry(p)) for p in fs}
    return len(j) == 1 and max(j) <= t and all((n + 1) % odd(fib_entry(p)) == 0 for p in fs)
cnt = mism = 0
for n in range(9, 1000000, 2):
    if isprime(n) or n % 5 == 0: continue
    f = factorint(n)
    if max(f.values()) > 1: continue
    fs = list(f)
    if jacobi(5, n) == -1:
        a = strong_lucas(n); b = half_check(n, fs, 'L')
        cnt += 1; mism += (a != b)
    a = strong_fermat2(n); b = half_check(n, fs, 'F')
    mism += (a != b)
print("half-test checks on squarefree n<1e6: lucas cases", cnt, "mismatches", mism)
