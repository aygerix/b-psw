"""Per-row exchange rate of the real instance. For row q^e: c_q = #{p : q | odd(ord_p 2) or q | odd(z(p))} (the primes
lost if the row is dropped), bits_q = log2 phi(q^e), rho_q = bits_q / log2(N / (N - c_q)) (group bits paid per doubling
of N). A row lowers the k-tree deficit (list size 2^b) iff rho_q < b. First-order (single-row) view of the min cut.
Usage: python3 rowrate.py instance.json"""
import sys, json, math
d = json.load(open(sys.argv[1])); N = d['N']
rows = [(int(q), e) for q, e in list(d['L1_factors'].items()) + list(d['L2_factors'].items())]
cnt = {q: 0 for q, _ in rows}
for ps, os_, zs in d['primes']:
    o, z = int(os_) >> 1, int(zs) >> 2
    for q in cnt:
        if o % q == 0 or z % q == 0: cnt[q] += 1
out = []
for q, e in rows:
    bits = math.log2((q - 1) * q ** (e - 1)); c = cnt[q]
    rho = bits / math.log2(N / (N - c)) if c < N else 0.0
    out.append((rho, q, e, c, bits))
out.sort()
for thr in (30, 50, 60, 100, 200, 400):
    sel = [x for x in out if x[0] < thr]
    print(f"rows with rho < {thr:3d}: {len(sel):3d} of {len(out)}, carrying {sum(x[4] for x in sel):7.1f} of {sum(x[4] for x in out):.1f} bits")
print("smallest rho:", [(round(r), q, c) for r, q, e, c, b in out[:8]])
print("median rho:", round(out[len(out)//2][0]), " largest rows (q > 1000): median rho",
      round(sorted(x[0] for x in out if x[1] > 1000)[len([x for x in out if x[1] > 1000])//2]),
      " mean support fraction", round(sum(x[3] for x in out if x[1] > 1000) / N / len([x for x in out if x[1] > 1000]), 4))
