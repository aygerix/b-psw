import sys, numpy as np
d = np.load(sys.argv[1])
rr, sep, gap, mg, mina = d['rr'], d['sep'], d['gap'], d['mg'], d['mina']
print("n", len(rr))
print("sep quantiles", np.quantile(sep, [0, .1, .5, .9]))
print("gap quantiles", np.quantile(gap, [0, .1, .5, .9]))
print("margin quantiles", np.nanquantile(mg, [0, .1, .5, .9]) if np.isfinite(mg).any() else None)
print("min|a| quantiles", np.quantile(mina, [0, .1, .5, .9]))
print("max|a| quantiles", np.quantile(np.abs(d['a']).max(axis=1), [0, .1, .5, .9]))
o = np.argsort(rr)[:int(sys.argv[2]) if len(sys.argv) > 2 else 5]
for i in o:
    print(f"rr {rr[i]:.3e} rem {d['rem'][i]:.3e} sep {sep[i]:.2e} gap {gap[i]:.2e} mg {mg[i]:.2e}")
    print("  a", np.round(d['a'][i].real, 5)); print("  c", np.round(d['c'][i].real, 5))
