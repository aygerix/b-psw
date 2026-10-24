"""Independent check of a BPSW subset-product instance written by instance.py (stdlib only, shares no code with it).

For every listed prime p, with claimed o = ord_p 2 and z = z(p) (Fibonacci rank), it checks:
  - p is prime (deterministic Miller-Rabin, valid below 3.3e24);
  - p == 3 (mod 8) and (5/p) = -1;
  - o and z are exact (2^o == 1 and 2^(o/q) != 1 for each prime q | o; the same for F_z mod p);
  - v2(o) = 1, v2(z) = 2, odd(o) | L1, odd(z) | L2, gcd(L1, L2) = 1, and p does not divide L1 L2;
  - L1 and L2 equal the products of their listed prime powers.
With a subset (comma-separated indices) it also checks the congruences and runs the exact test in src/bpsw.py.
Usage: python3 verify_instance.py instance.json [i1,i2,...]
"""
import sys, json, math, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

def is_prime(n):
    if n < 2: return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41):
        if n % q == 0: return n == q
    d, s = n - 1, 0
    while d % 2 == 0: d //= 2; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41):
        x = pow(a, d, n)
        if x in (1, n - 1): continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1: break
        else: return False
    return True

def jacobi(a, n):
    a %= n; r = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5): r = -r
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3: r = -r
        a %= n
    return r if n == 1 else 0

def prime_factors(m):
    fs, q = [], 2
    while q * q <= m:
        if m % q == 0:
            fs.append(q)
            while m % q == 0: m //= q
        q += 1
    if m > 1: fs.append(m)
    return fs

def fib(n, p):
    # F_n mod p via the matrix [[1,1],[1,0]]^n
    def mul(A, B):
        return [[(A[0][0]*B[0][0] + A[0][1]*B[1][0]) % p, (A[0][0]*B[0][1] + A[0][1]*B[1][1]) % p],
                [(A[1][0]*B[0][0] + A[1][1]*B[1][0]) % p, (A[1][0]*B[0][1] + A[1][1]*B[1][1]) % p]]
    R, M = [[1, 0], [0, 1]], [[1, 1], [1, 0]]
    while n:
        if n & 1: R = mul(R, M)
        M = mul(M, M); n >>= 1
    return R[0][1]

def v2(x): return (x & -x).bit_length() - 1

def main():
    inst = json.load(open(sys.argv[1]))
    L1, L2 = int(inst["L1"]), int(inst["L2"])
    assert L1 == math.prod(int(q) ** e for q, e in inst["L1_factors"].items())
    assert L2 == math.prod(int(q) ** e for q, e in inst["L2_factors"].items())
    assert math.gcd(L1, L2) == 1 and L1 % 2 and L2 % 2
    for k, (ps, os_, zs) in enumerate(inst["primes"]):
        p, o, z = int(ps), int(os_), int(zs)
        assert is_prime(p), p
        assert p % 8 == 3 and jacobi(5, p) == -1, p
        assert (p - 1) % o == 0 and pow(2, o, p) == 1, p
        assert all(pow(2, o // q, p) != 1 for q in prime_factors(o)), p
        assert (p + 1) % z == 0 and fib(z, p) == 0, p
        assert all(fib(z // q, p) != 0 for q in prime_factors(z)), p
        assert v2(o) == 1 and v2(z) == 2, p
        assert L1 % (o >> 1) == 0 and L2 % (z >> 2) == 0, p
        assert math.gcd(p, L1 * L2) == 1, p                      # unit condition (p must not divide the moduli)
    print(f"verified {len(inst['primes'])} primes; L1 {L1.bit_length()} bits, L2 {L2.bit_length()} bits, "
          f"gcd(L1, L2) = 1; claimed log2|G| = {inst['log2_G']}")
    if len(sys.argv) > 2:
        idx = [int(i) for i in sys.argv[2].split(",")]
        assert len(set(idx)) == len(idx) and len(idx) % 2 == 1 and len(idx) >= 3
        n = math.prod(int(inst["primes"][i][0]) for i in idx)
        c1, c2 = n % L1 == 1, n % L2 == L2 - 1
        from bpsw import bpsw_pass
        print(f"subset of {len(idx)} primes: n has {n.bit_length()} bits; n == 1 mod L1: {c1}; n == -1 mod L2: {c2}; "
              f"BPSW test passes: {bpsw_pass(n)}")

if __name__ == "__main__":
    main()
