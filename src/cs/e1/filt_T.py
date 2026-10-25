import sys
T = int(sys.argv[1]); PR = [q for q in range(3, 201) if all(q % d for d in range(2, int(q ** .5) + 1))]
def sm(x):
    for q in PR:
        if q > T: break
        while x % q == 0: x //= q
    return x == 1
for fn in sys.argv[2:]:
    for l in open(fn):
        a = l.split()
        if sm(int(a[5])) and sm(int(a[6])): sys.stdout.write(l)
