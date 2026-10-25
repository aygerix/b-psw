#!/usr/bin/env python3
"""restrict.py MODE IN.out OUT.out ORIG1.out [ORIG2.out ...]
mode 'orig':  mark as FAIL every exponent that the original run failed to factor (#FAILPHI lines), i.e. reproduce the
              original exponent set exactly.
mode 'extra': keep ONLY those exponents (the ones the original lost), mark every other exponent FAIL.
REC lines are kept unchanged (the checker ignores RECs whose exponent is FAIL)."""
import sys

mode, src, dst = sys.argv[1], sys.argv[2], sys.argv[3]
failed = set()
for fn in sys.argv[4:]:
    for line in open(fn):
        if line.startswith('#FAILPHI'):
            w = line.split()
            failed.add(({'e': 'F', 'f': 'L'}[w[1]], int(w[2])))
n_keep = n_drop = 0
with open(dst, 'w') as out:
    for line in open(src):
        w = line.split()
        if w and w[0] == 'PHI':
            key = (w[1], int(w[2]))
            drop = (key in failed) if mode == 'orig' else (key not in failed)
            if drop and w[3] != 'FAIL':
                line = 'PHI %s %s FAIL\n' % (w[1], w[2])
                n_drop += 1
            elif w[3] != 'FAIL':
                n_keep += 1
        out.write(line)
print('%s: kept %d factored exponents, masked %d (original failures: %d)' % (mode, n_keep, n_drop, len(failed)))
