#!/usr/bin/env python3
"""Williams / Korselt families whose special points are all rational (theta in {-2,-1,0,1,2}).
Members f in Q[Y]: f - 1 = c*A, f + 1 = c*B with A, B monic, roots in {-2..2}, disjoint root sets, B - A = 2/c
constant (so c = 2/(B-A)). Williams family: lcm(A_j) | n - 1 and lcm(B_j) | n + 1 (n = prod f_j), with
root sets of all A's disjoint from all B's. Korselt-only mode drops the B condition.
This is the 'torsion-only' corner of the zero-sum world (cyclotomic special points of order 1,2,3,4,6).
Usage: python3 ratw.py DEGMAX KMAX [korselt]"""
import sys, itertools
from fractions import Fraction
import sympy as sp
Y = sp.symbols('Y')
PTS = [-2, -1, 0, 1, 2]
DEG, KMAX = int(sys.argv[1]), int(sys.argv[2]); KOR = len(sys.argv) > 3

def polys(deg, pts):
    for ms in itertools.combinations_with_replacement(pts, deg):
        yield ms, sp.Poly(sp.prod([Y - t for t in ms]), Y)

members = []
for d in range(1, DEG + 1):
    for sa in range(1, len(PTS)):
        for SA in itertools.combinations(PTS, sa):
            SB = [t for t in PTS if t not in SA]
            for ma, A in polys(d, list(SA)):
                if set(ma) != set(SA): continue
                for sb in range(1, len(SB) + 1):
                    for SBs in itertools.combinations(SB, sb):
                        for mb, B in polys(d, list(SBs)):
                            if set(mb) != set(SBs): continue
                            diff = B - A
                            if diff.degree() > 0 or diff.is_zero: continue
                            k = Fraction(int(sp.numer(diff.as_expr())), int(sp.denom(diff.as_expr())))
                            c = Fraction(2) / k
                            if c <= 0: continue
                            f = sp.Poly(1 + sp.Rational(c.numerator, c.denominator) * A.as_expr(), Y)
                            members.append(dict(f=f, A=A, B=B, sa=set(ma), sb=set(mb), deg=d,
                                                irr=sp.Poly(f.as_expr(), Y).is_irreducible))
print(f"members (deg <= {DEG}): {len(members)}, irreducible: {sum(m['irr'] for m in members)}")
for m in members[:12]:
    print("  ", sp.factor(m['f'].as_expr()), "| f-1 roots", sorted(m['sa']), "| f+1 roots", sorted(m['sb']))
M = [m for m in members if m['irr']]

def compatible(a, b):
    return not (a['sa'] & b['sb'] or a['sb'] & b['sa']) if not KOR else True

found = []
def check(S):
    n = sp.prod([M[i]['f'].as_expr() for i in S])
    LA = sp.lcm_list([M[i]['A'].as_expr() for i in S])
    if sp.rem(sp.expand(n - 1), LA, Y) != 0: return False
    if KOR: return True
    LB = sp.lcm_list([M[i]['B'].as_expr() for i in S])
    return sp.rem(sp.expand(n + 1), LB, Y) == 0

def rec(start, S):
    if len(S) >= 2 and check(S):
        found.append(list(S))
    if len(S) == KMAX: return
    for j in range(start + 1, len(M)):
        if all(compatible(M[i], M[j]) for i in S):
            rec(j, S + [j])
rec(-1, [])
print(f"{'Korselt' if KOR else 'Williams'} families with k <= {KMAX}: {len(found)}")
for S in found[:30]:
    print("  ", [str(M[i]['f'].as_expr()) for i in S])
