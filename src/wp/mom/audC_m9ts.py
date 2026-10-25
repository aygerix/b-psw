# audC_m9ts.py -- two-sided single-point tests at stored m9 boundary / level-2 points (P-chart data U = u_i, V = v_i).
# For each member i and z in {w_i, v_i}: status of the other factors (all >= 1 / all <= 1 / straddle).
import re, numpy as np
def fv(a, c, t): return 2*t + 1 + t*(t+1)*(a*t + c)
def up(A, B, C):
    r = np.roots([A, B, C]); return r[np.argmax(r.imag)]
def lists(txt): return [[float(x) for x in l.split(",")] for l in re.findall(r"\[([-0-9eE., ]+)\]", txt)]
confs = []
for fn in ("m9_sol_9_0.gp", "m9_sol_9_1.gp", "m9_sol_11_0.gp", "m9_sol_11_1.gp", "m9_l2cert_11_0.gp"):
    try: L = lists(open(fn).read())
    except Exception as e: print(fn, e); continue
    L = [l for l in L if len(l) >= 8]
    for k in range(0, len(L) - 1, 2): confs.append((fn, L[k], L[k+1]))
txt = open("m9_l2cert_9.gp").read(); L = [l for l in lists(txt) if len(l) >= 8]
for k in range(0, len(L) - 1, 2): confs.append(("m9_l2cert_9.gp", L[k], L[k+1]))
for fn, U, V in confs:
    m = len(U)
    lin = (1 if fn.endswith("_1.gp") else 0) if fn.startswith("m9_sol") else (0 if fn == "m9_l2cert_11_0.gp" else (1 if m == 8 else 0))
    a = [1/v for v in V]; c = [u/v for u, v in zip(U, V)]
    inD = all(cc*cc < 8*aa and (aa - cc)**2 < 8*aa for aa, cc in zip(a, c))
    st = []
    for i in range(m):
        w = up(a[i], a[i] + c[i], c[i] + 2); v = up(a[i], c[i], 2.0)
        for nm, z in (("w", w), ("v", v)):
            Lg = [np.log(abs(fv(a[l], c[l], z))) for l in range(m) if l != i] + ([np.log(abs(2*z+1))] if lin else [])
            st.append(nm + ("ge" if min(Lg) >= 0 else ("le" if max(Lg) <= 0 else "S")))
    ok = [s for s in st if not s.endswith("S")]
    print("%-18s m=%d lin=%d inD=%s  passing two-sided tests: %d %s" % (fn, m, lin, inD, len(ok), ok))
