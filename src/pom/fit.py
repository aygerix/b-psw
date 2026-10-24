"""Free additions to an instance: primes whose odd orders fit its existing moduli (no change to log2|G|).

A candidate p fits instance (L1, L2) if it is prime, p == 3 (mod 8), (5/p) = -1, v2(ord_p 2) = 1, v2 z(p) = 2,
odd(ord_p 2) | L1, odd(z(p)) | L2 and gcd(p, L1 L2) = 1. Then p can join the pool and the instance claim still
holds verbatim.
Usage: python3 fit.py instance.json candidates.txt out_augmented.json [Tmax]
"""
import sys, json, math
from profile2 import primes_upto, factor_smooth, order_mod, zrank, v2

def main():
    inst = json.load(open(sys.argv[1]))
    L1, L2 = int(inst["L1"]), int(inst["L2"])
    P = primes_upto(int(sys.argv[4]) if len(sys.argv) > 4 else 5000)
    have = {int(r[0]) for r in inst["primes"]}
    add, seen, tried = [], set(), 0
    for line in open(sys.argv[2]):
        if not line.strip(): continue
        p = int(line)
        if p in have or p in seen or p < 5: continue
        seen.add(p); tried += 1
        if p % 8 != 3 or (L1 * L2) % p == 0: continue
        fm, fp = factor_smooth(p - 1, P), factor_smooth(p + 1, P)
        if max(fm) > P[-1] or max(fp) > P[-1]: continue          # p -+ 1 not smooth enough to factor here
        o, z = order_mod(2, p, fm), zrank(p, fp)
        if (v2(o), v2(z)) != (1, 2): continue
        if L1 % (o >> 1) or L2 % (z >> 2): continue
        add.append([str(p), str(o), str(z)])
    inst["primes"] = inst["primes"] + add
    inst["N"] = len(inst["primes"])
    inst["expected_solutions_log2"] = round(inst["N"] - 1 - inst["log2_G"], 2)
    inst["description"] += f" Augmented with {len(add)} further primes fitting the same moduli (src/pom/fit.py)."
    json.dump(inst, open(sys.argv[3], "w"), indent=0)
    print(f"candidates tried {tried}, fitting {len(add)}; augmented N = {inst['N']}, log2|G| = {inst['log2_G']}, "
          f"log2 E = {inst['expected_solutions_log2']}, ratio {inst['N'] / inst['log2_G']:.4f}")

if __name__ == "__main__":
    main()
