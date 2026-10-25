"""Merge bbc_fx runs of the sub-instances of one split instance into one run_bbc.sh-format record (for g80.py).
M_ub(lambda) = max over sub-instances (each configuration is feasible, with the same value, in the sub-instance matching
its side pattern on the split rows); printed only if every sub-instance finished (DONE and rc 0), else 'rc 1'.
Usage: python3 merge_split.py PREFIX "D j1 j2 pool" sub1,sub2,...   (files PREFIX_<sub>.out)"""
import sys
pre, head, subs = sys.argv[1], sys.argv[2].split(), sys.argv[3].split(',')
print("== " + " ".join(head) + " 0 0 0")
Ms, nodes, t, lam = [], 0, 0.0, None
for s in subs:
    try: L = open(f"{pre}_{s}.out").read().split('\n')
    except OSError: print("rc 1"); sys.exit()
    la = [l for l in L if l.startswith('lam')]; rc = [l for l in L if l.startswith('rc')]
    if not la or not rc or rc[-1].split()[1] != '0' or la[-1].split()[2] != 'DONE': print("rc 1"); sys.exit()
    a = la[-1].split(); lam = a[1]; Ms.append(float(a[4])); nodes += int(a[12]); t += float(a[18])
print(f"lam {lam} DONE M_ub {max(Ms):.6f} best {max(Ms):.6f} N 0 bits 0.000 nodes {nodes} exact 0 pruned 0 time {t:.2f} "
      f"sb 0 sbclosures 0 split {len(subs)} subM {','.join(f'{m:.6f}' for m in Ms)}")
print("rc 0")
