#!/usr/bin/env python3
"""Compare the adversarial per-prime forced records (adv_forced_records.txt) with the original fic.gp records.

Original record (23 fields): tag ef bits(p) bits(No) logidx Elog Vlog hq ehq o3 o10 o100 o1000 p3 p10 p100 p1000
                             U200 P200 U1500 P1500 U1e4 P1e4   (p itself is not recorded)
Matching key: (tag, ef, bits(p)); within a key the records are compared as sorted multisets.
Usage: cmp_orig.py adv_forced_records.txt orig1.out [orig2.out ...]
"""
import sys
from collections import defaultdict


def close(a, b, tol=2e-4):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def main():
    mine = defaultdict(list)
    for line in open(sys.argv[1]):
        w = line.split()
        tag, x, p, bp = w[0], int(w[1]), int(w[2]), int(w[3])
        L, EL, VL, hq, ehq = float(w[7]), float(w[8]), float(w[9]), int(w[10]), float(w[11])
        ot = [int(v) for v in w[12].split(',')]
        pt = [float(v) for v in w[13].split(',')]
        U = [int(v) for v in w[14].split(',')]
        P = [float(v) for v in w[15].split(',')]
        mine[(tag, x, bp)].append((L, EL, VL, hq, ehq, ot, pt, U, P, p))
    orig = defaultdict(list)
    for fn in sys.argv[2:]:
        for line in open(fn):
            if line.startswith('#'):
                continue
            f = line.split()
            if len(f) != 23 or f[0] not in ('FF', 'LF'):
                continue
            orig[(f[0], int(f[1]), int(f[2]))].append(
                (float(f[4]), float(f[5]), float(f[6]), int(f[7]), float(f[8]), [int(v) for v in f[9:13]],
                 [float(v) for v in f[13:17]], [int(f[17]), int(f[19]), int(f[21])],
                 [float(f[18]), float(f[20]), float(f[22])]))
    keys_o, keys_m = set(orig), set(mine)
    n_o = sum(len(v) for v in orig.values())
    n_m = sum(len(v) for v in mine.values())
    print('records: original %d, adversarial %d' % (n_o, n_m))
    only_o = sorted(keys_o - keys_m)
    only_m = sorted(keys_m - keys_o)
    print('keys only in original: %d %s' % (len(only_o), only_o[:20]))
    print('keys only in adversarial: %d %s' % (len(only_m), only_m[:20]))
    nbad = 0
    nsame = 0
    for k in sorted(keys_o & keys_m):
        a = sorted(orig[k], key=lambda r: (r[0], r[1]))
        b = sorted([r[:9] for r in mine[k]], key=lambda r: (r[0], r[1]))
        if len(a) != len(b):
            print('count mismatch', k, len(a), len(b))
            nbad += 1
            continue
        for ra, rb in zip(a, b):
            ok = (close(ra[0], rb[0]) and close(ra[1], rb[1]) and close(ra[2], rb[2]) and ra[3] == rb[3]
                  and close(ra[4], rb[4]) and ra[5] == rb[5] and all(close(x, y) for x, y in zip(ra[6], rb[6]))
                  and ra[7] == rb[7] and all(close(x, y) for x, y in zip(ra[8], rb[8])))
            if ok:
                nsame += 1
            else:
                nbad += 1
                if nbad <= 10:
                    print('MISMATCH', k, '\n  orig', ra, '\n  adv ', rb)
    print('matched records identical (all 18 statistic fields): %d; mismatched: %d' % (nsame, nbad))
    for k in only_m:
        for r in mine[k]:
            print('  adversarial-only record', k, 'p =', r[9], 'logidx %.4g' % r[0])


if __name__ == '__main__':
    main()
