"""FX (1)/(4): which D can the Selfridge rule of src/bpsw.py select?  (Corrected after the first version of this check
REFUTED the naive claim 'only prime |D| and -15': a skipped earlier D' may have Jacobi symbol 0, not +1.)
THEOREM (proof in the notes).  Selfridge returns the first odd k >= 5 with (n/k) = -1 (D = k* = (-1)^((k-1)/2) k, and
(k*/n) = (n/k) by reciprocity).  For composite k every prime factor q >= 5 of k was tested earlier, so (n/q) in {0, +1};
hence (n/k) = -1 forces gcd(n, k) = 1, (n/3) = -1 and v3(k) odd.  If (n/3) = -1 then k = 27 has (n/27) = -1, so the
selected k is <= 27.  Therefore the selected D is k* for a prime k >= 5, or -15 [(n/3) = -1, (n/5) = +1, (n/q) != -1 for
q = 7, 11, 13], or 21 [(n/3) = -1, 5 | n, 7 not | n], or -27 [(n/3) = -1, 35 | n].  For D = 21 (Q = -5) and D = -27 (Q = 7)
some prime of n divides Q, so U_k == V_k == 1 mod that prime for k >= 1 and the strong Lucas test FAILS.
So the D of any n passing the strong Lucas test of src/bpsw.py is k* (k prime >= 5) or -15.
Check: every odd non-square n in [3, X] with bpsw.selfridge_D verbatim; also runs bpsw.strong_lucas on every n whose D is
21 or -27 (must be False)."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from bpsw import selfridge_D, jacobi, strong_lucas
from collections import Counter
X = int(float(sys.argv[1]))
def isprime(k): return k > 1 and all(k % d for d in range(2, int(k ** .5) + 1))
def J(n, q): return jacobi(n, q)
cnt = Counter(); viol = Counter(); lucas_pass_bad = 0
for n in range(3, X + 1, 2):
    D = selfridge_D(n)
    if D is None: continue
    cnt[D] += 1
    k = abs(D)
    if isprime(k) and k >= 5: pass
    elif D == -15:
        if not (J(n, 3) == -1 and J(n, 5) == 1 and all(J(n, q) != -1 for q in (7, 11, 13))): viol['-15 cond'] += 1
    elif D == 21:
        if not (J(n, 3) == -1 and n % 5 == 0 and n % 7 != 0): viol['21 cond'] += 1
        if strong_lucas(n): lucas_pass_bad += 1
    elif D == -27:
        if not (J(n, 3) == -1 and n % 35 == 0): viol['-27 cond'] += 1
        if strong_lucas(n): lucas_pass_bad += 1
    else: viol['unexpected D'] += 1
    if n % 3 and J(n, 3) == -1 and k > 27: viol['(n/3)=-1 but |D|>27'] += 1
print(f"odd non-square n <= {X}: violations {dict(viol)} ; strong Lucas passes with D in {{21, -27}}: {lucas_pass_bad}")
print("selected D (count):", sorted(cnt.items(), key=lambda t: abs(t[0])))
