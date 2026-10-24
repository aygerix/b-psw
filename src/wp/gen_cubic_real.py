"""Generate Singular input for REAL CUBIC Williams families (positive leading coefficients), any odd k.

Structure (Lemma 24 in the notebook): by Lemma 6 every member shares the real special points. After an affine change of t:
    f_j = a_j t^3 + b_j t^2 + (b_j - a_j + 2) t + 1,   so f_j(0) = 1, f_j(-1) = -1,
    f_j - 1 = t g_j,  g_j = a_j t^2 + b_j t + (b_j - a_j + 2),
    f_j + 1 = (t + 1) h_j,  h_j = a_j t^2 + (b_j - a_j) t + 2.
Complex special points cannot be shared by distinct members (g_i ~ g_j forces f_i = f_j via g(-1) = 2; h_i ~ h_j via h(0) = 2).
Williams conditions: g_j | (F - 1)/t and h_j | (F + 1)/(t + 1) in C[t] (the real points 0, -1 are simple and common).
Remainders are computed over the parameter field, numerators are mapped to a polynomial ring, then std over GF(p).
Distinct members: 1 - y prod_{i<j} ((a_i - a_j) + c (b_i - b_j)) with a random constant c; a_j != 0 via the same y.
Usage: python3 gen_cubic_real.py k prime c > file.sing
"""
import sys
k = int(sys.argv[1]); P = int(sys.argv[2]); c = int(sys.argv[3])
A = [f"a{j}" for j in range(1, k + 1)]; B = [f"b{j}" for j in range(1, k + 1)]
params = ",".join(A + B)
L = []
L.append(f'ring R0 = ({P},{params}),(t),dp;')
L.append('poly F = 1;')
for j in range(k):
    L.append(f'poly f{j+1} = {A[j]}*t^3 + {B[j]}*t^2 + ({B[j]} - {A[j]} + 2)*t + 1; F = F*f{j+1};')
L.append('poly Em = (F - 1)/t; poly Ep = (F + 1)/(t + 1);')
L.append('ideal eqs;')
for j in range(k):
    a, b = A[j], B[j]
    L.append(f'poly g = t^2 + ({b}/{a})*t + ({b} - {a} + 2)/{a}; poly h = t^2 + (({b} - {a})/{a})*t + 2/{a};')
    L.append('poly rg = reduce(Em, std(g)); poly rh = reduce(Ep, std(h));')
    L.append('eqs = eqs, numerator(leadcoef(jet(rg,1) - jet(rg,0))), numerator(leadcoef(jet(rg,0) + 0*t)),'
             ' numerator(leadcoef(jet(rh,1) - jet(rh,0))), numerator(leadcoef(jet(rh,0) + 0*t));')
L.append(f'ring R1 = {P},({params},y),dp;')
L.append('ideal eqs = imap(R0, eqs);')
dist = "*".join(f"(({A[i]} - {A[j]}) + {c}*({B[i]} - {B[j]}))" for i in range(k) for j in range(i + 1, k))
prodA = "*".join(A)
L.append(f'ideal I = eqs, 1 - y*{prodA}*{dist};')
L.append('option(redSB); int t0 = timer; ideal G = std(I);')
L.append(f'print("k = {k}, real cubic ansatz, distinct members, GF({P}), c = {c}: dim = " + string(dim(G)) + ", size = " + string(size(G)));')
L.append('if (size(G) == 1 && G[1] == 1) { print("UNIT IDEAL: no family with distinct members (over the algebraic closure)"); }')
L.append('else { print("vdim = " + string(vdim(G))); print(G); }')
L.append('quit;')
print("\n".join(L))
