"""Turn DUMP lines of aud_ls.py into PARI calls of verify_cfg.gp.  Usage: python3 dump2gp.py D j1 j2 dumpfile > x.gp"""
import sys, re, ast
D, j1, j2, fn = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
print('read("verify_cfg.gp");')
for l in open(fn):
    if not l.startswith('DUMP'): continue
    m = re.match(r"DUMP b<=(\d+): N (\d+) bits ([\d.]+) E1 (\[.*?\]) E2 (\[.*?\]) primes (\[.*\])", l.strip())
    b, N, bits, E1, E2, L = m.groups()
    E1, E2, L = ast.literal_eval(E1), ast.literal_eval(E2), ast.literal_eval(L)
    gv = lambda E: '[' + ', '.join(f'[{q}, {e}]' for q, e in E) + ']'
    print(f'print("claimed b <= {b}: N {N}, bits {bits}"); vcfg({D}, {j1}, {j2}, {gv(E1)}, {gv(E2)}, {L});')
