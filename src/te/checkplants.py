"""checkplants.py -- independent (pure Python) re-check of a plant file written by mkplants.gp (format: SPEC.md).

Usage: python3 checkplants.py D X FILE [FILE ...]

For every "T p1 p2 p3 S A B" line it re-derives everything from the text alone:
  - p1 < p2 < p3 <= X, each a prime with Jacobi(D, p) = -1 and p not dividing 2QD (Q = (1-D)/4);
  - no prime is used by two plants; every planted prime has exactly one "O p M" line and vice versa;
  - M | p^2 - 1 and (W1) M | n/p - 1 at all three primes, n = p1 p2 p3 S;
  - M = gcd(p^2-1, n/p-1) (header mrule=gcd) or its lcm-form part (header mrule=lcmform);
  - S = 1 with A = B = 1, or S = product of the 1-2 distinct split primes listed in the comment (Sfac=...),
    A, B nontrivial, (W2) A | n-1, B | n+1, gcd(A, B) | 2;
  - the class thresholds named in the preceding comment line ("class=...").
Exit status 0 iff no errors.
"""
import sys
from math import gcd

MR_BASES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41)   # deterministic for n < 3.3e24


def is_prime(n):
    if n < 2:
        return False
    for b in MR_BASES:
        if n % b == 0:
            return n == b
    if n >= 3_317_044_064_679_887_385_961_981:
        raise ValueError("is_prime: n too large for the deterministic base set")
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for b in MR_BASES:
        x = pow(b, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


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


def v2(x):
    return (x & -x).bit_length() - 1


def lcmform(p, M):
    cut = max(0, v2(M) - max(v2(p - 1), v2(p + 1)))
    return M >> cut


def lcm(*xs):
    r = 1
    for x in xs:
        r = r * x // gcd(r, x)
    return r


def check(D, X, path):
    Q = (1 - D) // 4
    bad = []
    O, T, mrule = {}, [], "gcd"
    cls, sfac = None, None
    for ln, line in enumerate(open(path), 1):
        w = line.split()
        if not w:
            continue
        if w[0] == "#":
            for t in w[1:]:
                if t == "mrule=lcmform":
                    mrule = "lcmform"
                if t.startswith("class="):
                    cls = t[6:]
                if t.startswith("Sfac="):
                    sfac = t[5:]
            continue
        if w[0] == "O":
            p, M = int(w[1]), int(w[2])
            if p in O:
                bad.append(f"line {ln}: duplicate override {p}")
            O[p] = M
        elif w[0] == "T":
            T.append((ln, tuple(int(x) for x in w[1:7]), cls, sfac))
            cls, sfac = None, None
        else:
            bad.append(f"line {ln}: unparsed")

    def inert(p):
        return p > 2 and (2 * Q * D) % p != 0 and is_prime(p) and jacobi(D, p) == -1

    def split(q):
        return q > 2 and (2 * Q * D) % q != 0 and is_prime(q) and jacobi(D, q) == 1

    used, counts = set(), {}
    for ln, (p1, p2, p3, S, A, B), c, sf in T:
        err = lambda m: bad.append(f"line {ln} T {p1} {p2} {p3} {S} {A} {B} [{c}]: {m}")
        n = p1 * p2 * p3 * S
        if not (3 <= p1 < p2 < p3 <= X):
            err("order/range")
        for p in (p1, p2, p3):
            if not inert(p):
                err(f"{p} not an inert prime")
            if p in used:
                err(f"{p} reused")
            used.add(p)
            if p not in O:
                err(f"no override for {p}")
                continue
            M, g = O[p], gcd(p * p - 1, n // p - 1)
            if (p * p - 1) % M:
                err(f"M({p}) does not divide p^2-1")
            if (n // p - 1) % M:
                err(f"W1 fails at {p}")
            if M != (g if mrule == "gcd" else lcmform(p, g)):
                err(f"M({p}) does not follow mrule={mrule}")
        if S == 1:
            if (A, B) != (1, 1):
                err("S = 1 needs A = B = 1")
        else:
            qs = [int(x) for x in (sf or "").split("*")] if sf and sf != "1" else []
            prod = 1
            for q in qs:
                prod *= q
            if not (1 <= len(qs) <= 2 and prod == S and len(set(qs)) == len(qs) and all(split(q) for q in qs)):
                err(f"S is not the product of the listed distinct split primes ({sf})")
            if not (3 <= A < n - 1 and 3 <= B < n + 1):
                err("A or B trivial")
            if (n - 1) % A or (n + 1) % B:
                err("W2 fails")
            if 2 % gcd(A, B):
                err("gcd(A,B) does not divide 2")
        if c is None:
            err("no class comment")
            continue
        base = c.split("/")[-1]
        if c.startswith("d/") != (S != 1):
            err("class d iff S > 1 violated")
        if any(p not in O for p in (p1, p2, p3)):
            continue
        M2, M3 = O[p2], O[p3]
        idx3 = (p3 * p3 - 1) // M3
        L23 = lcm(M2, M3)
        idxmin = min(10**8, (p3 * p3 - 1) // 48)
        if base == "a":
            ok = idx3 <= 1000
        elif base == "b":
            ok = idx3 >= idxmin and p2 <= 100 * L23
        elif base in ("c", "c!"):
            ok = idx3 >= idxmin and p2 >= min(1000, p2 // 48) * L23 and (base == "c" or L23 <= 12)
        else:
            ok = False
        if not ok:
            err(f"class {c} thresholds fail (idx3={idx3}, p2/lcm(M2,M3)={p2 // L23})")
        counts[c] = counts.get(c, 0) + 1
    if set(O) != used:
        bad.append(f"overrides for unplanted primes: {sorted(set(O) - used)[:10]}")
    big = sum(1 for _, t, _, _ in T if t[0] * t[1] * t[2] * t[3] >= 2**128)
    nonlcm = sum(1 for p, M in O.items() if lcmform(p, M) != M)
    passes = len({t[3:] for _, t, _, _ in T})
    print(f"{path}: D={D} X={X} mrule={mrule} plants={len(T)} errors={len(bad)} classes={dict(sorted(counts.items()))} "
          f"passes(S,A,B)={passes} n>=2^128={big} non-lcm-form overrides={nonlcm}")
    for b in bad[:20]:
        print("  ", b)
    return len(bad)


if __name__ == "__main__":
    D, X = int(sys.argv[1]), int(float(sys.argv[2]))
    nbad = sum(check(D, X, f) for f in sys.argv[3:])
    sys.exit(1 if nbad else 0)
