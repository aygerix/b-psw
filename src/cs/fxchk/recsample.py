import os
import random, sys
B = os.path.expanduser('~/b-psw/src/cs/')
random.seed(20261025)
out = open('recsample.gp', 'w'); out.write('{R = ['); first = True
for fn, k in ((B + 'fx/chm200_pool_xD.txt', 500), (B + 'fx/c3_xD_1e10.txt', 400), (B + 'fx/c3_T300_1e10.txt', 300), (B + 'e1/chm200_pool.txt', 200)):
    L = [l.split() for l in open(fn) if len(l.split()) == 7]
    # bias: include the largest primes and the class (2,0) of D = 37 / -15
    S = random.sample(L, k) + sorted(L, key=lambda a: -int(a[1]))[:50]
    for a in S:
        out.write(('' if first else ',') + '[' + ','.join(a) + ']'); first = False
out.write(']};\n'); out.close()
