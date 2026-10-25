"""Compare census3 records with p in [a, b) against an audit record file (exact 7-tuples).
Usage: python3 cmp_range.py c3file auditfile a b"""
import sys
c3f, af, a, b = sys.argv[1], sys.argv[2], int(float(sys.argv[3])), int(float(sys.argv[4]))
def rd(fn, filt):
    s = set()
    for l in open(fn):
        t = tuple(map(int, l.split()))
        if len(t) == 7 and (not filt or a <= t[1] < b): s.add(t)
    return s
A = rd(c3f, True); B = rd(af, False)
print(f"[{a}, {b}): census3 {len(A)} records, audit {len(B)} records; census3-only {len(A - B)}, audit-only {len(B - A)}; "
      f"split {sum(1 for t in B if t[2])}, inert {sum(1 for t in B if not t[2])}")
for t in sorted(A - B)[:10]: print("  census3-only", t)
for t in sorted(B - A)[:10]: print("  audit-only", t)
