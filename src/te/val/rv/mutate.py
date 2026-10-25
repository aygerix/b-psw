# review-only: inject a known bug into anchor3_hook.c to check that the harness detects it (power check)
import sys
src = open(sys.argv[1]).read()
muts = {
 'sieve': ('for(u64 d = d0; d < len; d += l){ if(x != t) dead[d] = 1;', 'for(u64 d = d0; d < len; d += l){ if(x == t) dead[d] = 1;'),
 'm128': ('p = (u128)a1 * b0 + t0; t0 = (u64)p; q = (p >> 64) + (u128)a1 * b1 + t1; t1 = (u64)q; t2 += (u64)(q >> 64);', 'p = (u128)a1 * b0 + t0; t0 = (u64)p; q = (p >> 64) + (u128)a1 * b1 + t1; t1 = (u64)q;'),
 'runcS': ('    r = PC.S1 ? inv : mulmod128(inv, s3, M3);', '    r = inv;'),
 'lastp1': ('    for(; p1 < p2; p1 += L){', '    for(; p1 + L < p2 || (p1 < p2 && L > p2); p1 += L){'),
}
a, b = muts[sys.argv[3]]
if src.count(a) != 1: sys.exit('pattern count %d' % src.count(a))
open(sys.argv[2], 'w').write(src.replace(a, b)); print('ok', sys.argv[3])
