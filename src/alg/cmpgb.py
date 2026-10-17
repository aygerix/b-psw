# Exact comparison of reduced Groebner bases printed by Singular and by Calyx.
# Each polynomial is made monic w.r.t. grevlex in the ring's generator order, then the
# two sets are compared over Q. Usage: cmpgb.py SINGULAR_OUT CALYX_OUT
import sys
from sympy import symbols, Poly, QQ
from sympy.parsing.sympy_parser import parse_expr

def block(lines, header, size=None):
    i = next(k for k, l in enumerate(lines) if l.startswith(header))
    if size is None:
        size = int(lines[i + 1]); i += 1
    return lines[i + 1:i + 1 + size]

def canon(polys, gens):
    loc = {str(g): g for g in gens}
    out = []
    for s in polys:
        p = Poly(parse_expr(s.replace("^", "**"), local_dict=loc), *gens, domain=QQ)
        lc = p.coeffs(order="grevlex")[0]
        p = p.quo_ground(lc)
        out.append(tuple(p.terms(order="grevlex")))
    return out

S_ = [l.strip() for l in open(sys.argv[1])]
C_ = [l.strip() for l in open(sys.argv[2])]
R = symbols("y z c1 c2 c3 x"); S = symbols("c1 c2 c3")
ok_all = True
for name, sh, ch, gens, ssize in [("GB(I) in R", "== GB(I) dp", "== GB(I) grevlex", R, None),
                                 ("GB(eliminate(I,y*z)) in R", "== GB(K) dp", "== GB(K) grevlex", R, None),
                                 ("GB(C) in S", "== GB(C) in S", "== GB(C) in S", S, None),
                                 ("minimal prime (Singular minAssGTZ) vs GB(Q)", "-- prime 1", "== GB(Q) grevlex", S, 2)]:
    a = canon(block(S_, sh, ssize), gens); b = canon(block(C_, ch), gens)
    same = set(a) == set(b) and len(a) == len(b)
    ok_all &= same
    print(f"{name}: singular {len(a)} elements, calyx {len(b)} elements, identical as monic sets: {same}")
    if not same:
        for t in sorted(set(a) - set(b)): print("  only singular:", Poly.from_dict(dict(t), *gens).as_expr())
        for t in sorted(set(b) - set(a)): print("  only calyx:   ", Poly.from_dict(dict(t), *gens).as_expr())
print("ALL IDENTICAL" if ok_all else "DIFFERENCES FOUND")
