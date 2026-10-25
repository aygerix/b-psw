# audC_census.py -- two-sided single-point test at the ns census of real '-' half solutions (k = 7, 9, 11).
# For each stored solution (members (a_i, c_i), optional 2t+1): at each member i and z in {w_i, v_i}, moduli of the
# other factors; status 'ge' (all >= 1), 'le' (all <= 1) or 'straddle'.  A '-' half solution has prod = 1 at every w_i.
# If every (i, z) straddles, the configuration is UNCOVERED by two-sided single-point tests (positive control).
import sys, numpy as np, glob
def fv(a, c, t): return 2*t + 1 + t*(t+1)*(a*t + c)
def up(A, B, C):
    r = np.roots([A, B, C]); return r[np.argmax(r.imag)]
for fn in sorted(glob.glob("ns_out/census_k*_l*.npz")):
    d = np.load(fn); lin = int(fn.split("_l")[1][0])
    A = d['a']; C = d['c']; rem = d['remp'] if 'remp' in d.files else np.zeros(len(A))
    # distinct solutions (up to permutation): round sorted member lists
    seen = {}; 
    for j in np.argsort(rem):
        a = np.real(A[j]); c = np.real(C[j])
        if np.max(np.abs(np.imag(A[j]))) > 1e-9 or np.max(np.abs(np.imag(C[j]))) > 1e-9: continue
        key = tuple(sorted(zip(np.round(a, 5), np.round(c, 5))))
        if key not in seen: seen[key] = (a, c, rem[j])
    print("== %s: %d stored, %d distinct real; keys %s" % (fn, len(A), len(seen), d.files))
    for key, (a, c, r) in list(seen.items())[:12]:
        k = len(a); stat = []; prodw = []
        for i in range(k):
            w = up(a[i], a[i] + c[i], c[i] + 2); v = up(a[i], c[i], 2.0)
            for nm, z in (("w", w), ("v", v)):
                L = [np.log(abs(fv(a[l], c[l], z))) for l in range(k) if l != i] + ([np.log(abs(2*z+1))] if lin else [])
                s = "ge" if min(L) >= 0 else ("le" if max(L) <= 0 else "S")
                stat.append(nm + s)
                if nm == "w": prodw.append(abs(sum(L)))
        npass = sum(1 for s in stat if not s.endswith("S"))
        print("  rem %.2e  max|sum log| at w = %.1e  passing (i,z) tests: %d  v-statuses: %s" % (r, max(prodw), npass, [s for s in stat if s[0] == "v"]))
