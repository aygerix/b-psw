"""Export class pools for the C enumerator (c3/enum.c).

For each Selfridge D and each 2-adic class (j1, j2) of the pool built by exact.build (inert primes and
admissible split primes with p^2 - 1 y-smooth, Lemma 12(b) filter), write

  CLASS D j1 j2 k nq
  Q  m_1 ... m_nq                (moduli: q^emax for every odd q in the class supports, and 2^(max(j1,j2)+1))
  then k lines, one per prime:
  P  idx eps  r_1 .. r_nq  nF (qi e)*  nL (qi e)*  conflict-list
where r_i = p mod m_i, (qi e) are requirements q_qi^e | odd(ord) (F side) or odd(z) (L side), and the
conflict list holds the indices j of primes that cannot coexist with p (opposite-side support overlap).
Usage: python3 export.py y pairs_file out_file D1 [D2 ...]
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from collections import defaultdict
from sympy import factorint
import exact

def main():
    y, f, out = int(sys.argv[1]), sys.argv[2], sys.argv[3]
    args = sys.argv[4:]
    if '--mode' in args:
        i = args.index('--mode'); exact.MODE = args[i + 1]; del args[i:i + 2]
    Ds = [int(a) for a in args]
    xs = [int(l) for l in open(f) if l.strip()]
    fo = open(out, 'w')
    for D in Ds:
        rows = exact.build(xs, D)
        by = defaultdict(list)
        for r in rows:
            by[r['cls']].append(r)
        for c, R in sorted(by.items()):
            j1, j2 = c
            emax = defaultdict(int)
            for r in R:
                r['F'] = factorint(r['O1']); r['L'] = factorint(r['O2'])
                for q, e in list(r['F'].items()) + list(r['L'].items()):
                    emax[q] = max(emax[q], e)
            qs = sorted(emax)
            mods = [q ** emax[q] for q in qs] + [1 << (max(j1, j2) + 1)]
            qi = {q: i for i, q in enumerate(qs)}
            fo.write(f"CLASS {D} {j1} {j2} {len(R)} {len(mods)}\n")
            fo.write("Q " + " ".join(map(str, mods)) + "\n")
            for i, r in enumerate(R):
                conf = [j for j, s in enumerate(R) if j != i and ((r['S1'] & s['S2']) or (r['S2'] & s['S1']))]
                F = " ".join(f"{qi[q]} {e}" for q, e in sorted(r['F'].items()))
                L = " ".join(f"{qi[q]} {e}" for q, e in sorted(r['L'].items()))
                fo.write(f"P {r['p']} {r['eps']} " + " ".join(str(r['p'] % m) for m in mods) +
                         f" {len(r['F'])} {F} {len(r['L'])} {L} {len(conf)} " + " ".join(map(str, conf)) + "\n")
        print(f"D={D}: {len(rows)} primes, {len(by)} classes, max class {max(map(len, by.values())) if by else 0}",
              flush=True)
    fo.close()

if __name__ == "__main__":
    main()
