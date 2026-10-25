# Build a validation-only variant of anchor3.c: in fermat/lucas modes the final predicate becomes
# "(W1) with the one-torus M holds at all three primes" and branch A drops the Fermat prefilter/sieve/lanes.
# Its hit set must equal the brute-force weakord (fermat) / weakz (lucas) lists.
import sys
src = open(sys.argv[1]).read()
reps = [
 ('    case MODE_FERMAT: make_n(ts->n, p1, p2, p3); if(sprp2_mpz(&ts->b, ts->n)) record_hit("FERMAT", ts, p1, p2, p3, 0, branch); break;\n',
  '    case MODE_FERMAT: case MODE_LUCAS: if(w1_holds(get_M(i1), p2, p3) && w1_holds(get_M(i2), p1, p3) && w1_holds(get_M(i3), p1, p2)){ make_n(ts->n, p1, p2, p3); record_hit("W1ONE", ts, p1, p2, p3, 0, branch); } break;\n'),
 ('    case MODE_LUCAS: make_n(ts->n, p1, p2, p3); if(slprp_fixed_mpz(&ts->b, ts->n, (long)gD)) record_hit("LUCAS", ts, p1, p2, p3, 0, branch); break;\n', ''),
 ('  int usesieve = gSieve && (gMode == MODE_BPSW || gMode == MODE_FERMAT) && len >= 4 && step;\n', '  int usesieve = gSieve && (gMode == MODE_BPSW) && len >= 4 && step;\n'),
 ('  int lanes = PC.S1 && (gMode == MODE_BPSW || gMode == MODE_FERMAT);\n', '  int lanes = PC.S1 && (gMode == MODE_BPSW);\n'),
 ('    if(gMode == MODE_BPSW || gMode == MODE_FERMAT){\n      int nl = build_nm1(p3, m, ts->e);\n', '    if(gMode == MODE_BPSW){\n      int nl = build_nm1(p3, m, ts->e);\n'),
]
for a, b in reps:
    c = src.count(a)
    if c != 1: sys.exit("pattern count %d: %r" % (c, a[:70]))
    src = src.replace(a, b)
open(sys.argv[2], 'w').write(src)
print("ok")
