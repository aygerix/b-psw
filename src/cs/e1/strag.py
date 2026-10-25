"""Completeness audit of the E1 pool (HEURISTIC extrapolation; the counts themselves are COMPUTED).
Per D and decade: CHM records (p - 1, p + 1 both 200-smooth), census3 records (p < 1e10, true orders 200-smooth),
index stragglers = census3 records not in CHM.  Extrapolated stragglers above 1e10 = sum over decades >= 1e10 of
CHM(decade) * s, with s = stragglers/CHM over [1e8, 1e10).
Usage: python3 strag.py census3_file chm_file"""
import sys, math
from collections import defaultdict

def main():
    c3, chm = defaultdict(set), defaultdict(set)
    for fn, dst in ((sys.argv[1], c3), (sys.argv[2], chm)):
        for l in open(fn):
            a = l.split(); dst[int(a[0])].add(int(a[1]))
    for D in sorted(c3, key=lambda d: (abs(d), d)):
        dec_c3, dec_chm, dec_st = defaultdict(int), defaultdict(int), defaultdict(int)
        for p in c3[D]:
            k = int(math.log10(p)); dec_c3[k] += 1
            if p not in chm[D]: dec_st[k] += 1
        for p in chm[D]: dec_chm[int(math.log10(p))] += 1
        lo = sum(dec_chm[k] for k in (8, 9)); st = sum(dec_st[k] for k in (8, 9))
        s = st / lo if lo else 0
        above = sum(v for k, v in dec_chm.items() if k >= 10)
        missing_chm = [p for p in chm[D] if p < 10 ** 10 and p not in c3[D]]
        print(f"D = {D}: census3 {len(c3[D])}, CHM {len(chm[D])} (CHM above 1e10: {above}); CHM primes < 1e10 absent "
              f"from census3: {len(missing_chm)}; stragglers < 1e10: {sum(dec_st.values())}; s(1e8..1e10) = {s:.3f}; "
              f"extrapolated stragglers above 1e10 ~ {s * above:.0f}")
        print("   decade: CHM / census3 / stragglers  " + "  ".join(
            f"1e{k}:{dec_chm[k]}/{dec_c3[k]}/{dec_st[k]}" for k in range(0, 20) if dec_chm[k] or dec_c3[k]))

if __name__ == "__main__":
    main()
