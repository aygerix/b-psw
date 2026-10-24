"""Build the explicit BPSW subset-product instance from a Pomerance pool (D = 5, class (1,2)).

Takes a subpool S' (e.g. the lambda = 1 min-cut optimum from mincut.py) and writes JSON with the primes, their odd orders, L1, L2 and the claim:

  For any S in S' with |S| odd and >= 3, if prod(S) == 1 (mod L1) and prod(S) == -1 (mod L2),
  then n = prod(S) is composite and passes the strong Baillie-PSW test (Selfridge D = 5, P = 1, Q = -1).

Proof (Lemma 1 specialised). Every p in S' is a prime == 3 (mod 8) with (5/p) = -1, v2(ord_p 2) = 1 and
v2 z(p) = 2, odd(ord_p 2) | L1 and odd(z(p)) | L2.
 - (5/n) = (-1)^|S| = -1, so D = 5 is selected. n is squarefree, hence not a square.
 - Fermat: n == 3 (mod 8), so n - 1 = 2d with d odd and odd(ord_p 2) | d. Then 2^d has order exactly 2 mod p,
   so 2^d == -1 mod every p, hence mod n.
 - Lucas: n + 1 = 4d' with d' odd and odd(z(p)) | d'. The rank z(p) = 4 z' is the order of gamma = alpha/beta,
   so gamma^(2d') = -1 mod p, i.e. V_(2d') == 0 mod every p, hence mod n (strong test with r = 1 < s = 2).
Usage: python3 instance.py subpool.txt out.json [Tmax]
"""
import sys, json, math
from profile2 import primes_upto, factor_smooth, order_mod, zrank, v2, bits

def main():
    """Write the instance for a given subpool (one prime per line), e.g. the lambda = 1 output of mincut.py."""
    P = primes_upto(int(sys.argv[3]) if len(sys.argv) > 3 else 5000)
    sel = []
    for p in sorted(int(l) for l in open(sys.argv[1]) if l.strip()):
        if p < 5: continue
        fm, fp = factor_smooth(p - 1, P), factor_smooth(p + 1, P)
        o, z = order_mod(2, p, fm), zrank(p, fp)
        assert (v2(o), v2(z)) == (1, 2), p
        F = factor_smooth(o >> 1, P); L = factor_smooth(z >> 2, P)
        assert all(q % 4 == 1 for q in F) and all(q % 4 == 3 for q in L), p
        sel.append((p, o, z))
    L1f, L2f = {}, {}
    for p, o, z in sel:
        for q, e in factor_smooth(o >> 1, P).items(): L1f[q] = max(L1f.get(q, 0), e)
        for q, e in factor_smooth(z >> 2, P).items(): L2f[q] = max(L2f.get(q, 0), e)
    L1 = math.prod(q ** e for q, e in L1f.items()); L2 = math.prod(q ** e for q, e in L2f.items())
    assert math.gcd(L1, L2) == 1
    bad = [p for p, _, _ in sel if (L1 * L2) % p == 0]
    assert not bad, f"unit condition violated by {bad} (run mincut.py, which removes them)"
    log2G = sum(bits(q, e) for q, e in list(L1f.items()) + list(L2f.items()))
    out = {
        "description": "Explicit subset-product instance for a Baillie-PSW counterexample (D = 5, class (1,2)). "
                       "Any S subset of 'primes' with |S| >= 3, prod(S) == 1 mod L1 and prod(S) == -1 mod L2 "
                       "gives a composite n = prod(S) passing strong BPSW. |S| is then odd because 3 | L2 and every "
                       "listed prime is == -1 mod 3. Proof in src/pom/instance.py. No listed prime divides L1*L2.",
        "pool_file": sys.argv[1], "N": len(sel), "log2_G": round(log2G, 4),
        "expected_solutions_log2": round(len(sel) - 1 - log2G, 2),
        "L1": str(L1), "L2": str(L2),
        "L1_factors": {str(q): e for q, e in sorted(L1f.items())},
        "L2_factors": {str(q): e for q, e in sorted(L2f.items())},
        "primes": [[str(p), str(o), str(z)] for p, o, z in sel],
    }
    json.dump(out, open(sys.argv[2], "w"), indent=0)
    print(f"instance: N = {len(sel)}, log2|G| = {log2G:.4f}, log2 E(solutions) = {len(sel) - 1 - log2G:.2f}, "
          f"L1 {L1.bit_length()} bits ({len(L1f)} prime powers), L2 {L2.bit_length()} bits ({len(L2f)} prime powers)")

if __name__ == "__main__":
    main()
