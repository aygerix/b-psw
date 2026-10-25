"""Singular input for REAL CUBIC Williams families, Laurent formulation (session 8; see NOTEBOOK "Laurent form").
f_j = 2t + 1 + t(t+1)(a_j t + c_j); g_j = (f_j - 1)/t = a_j t^2 + (a_j + c_j) t + (c_j + 2); h_j = (f_j + 1)/(t+1) = a_j t^2 + c_j t + 2.
Since non-real special points are never shared (Lemma 24c), Williams <=> F - 1 = t prod g_j W_- and F + 1 = (t+1) prod h_j W_+,
with W_- monic of degree k-1 in t whose t^(k-2), t^(k-3) coefficients vanish, and W_+ the same in s = t + 1.
Unknowns: a_j, c_j, w_0..w_(k-4), v_0..v_(k-4); equations: all t-coefficients of the two identities.
Exclusions: prod a_j != 0 and pairwise distinctness via ONE linear form (a_i - a_j) + L (c_i - c_j); run with two
values of L to cover every distinct pair (a pair killed by both forms has a_i = a_j and c_i = c_j).
Usage: python3 gen_cubic_laurent.py k prime L > file.sing
"""
import sys
k, P, L = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
A = [f"a{j}" for j in range(1, k + 1)]; C = [f"c{j}" for j in range(1, k + 1)]
W = [f"w{i}" for i in range(k - 3)]; V = [f"v{i}" for i in range(k - 3)]
vars_ = A + C + W + V + ["y"]
print(f"ring R = {P}, (t, {', '.join(vars_)}), (dp(1), dp({len(vars_)}));")
print("poly F = 1; poly G = 1; poly H = 1;")
for j in range(1, k + 1):
    print(f"poly f{j} = 2*t + 1 + t*(t+1)*(a{j}*t + c{j});")
    print(f"F = F*f{j}; G = G*(a{j}*t^2 + (a{j} + c{j})*t + (c{j} + 2)); H = H*(a{j}*t^2 + c{j}*t + 2);")
wm = f"t^{k-1}" + "".join(f" + w{i}*t^{i}" for i in range(k - 3))
wp = f"(t+1)^{k-1}" + "".join(f" + v{i}*(t+1)^{i}" for i in range(k - 3))
print(f"poly Wm = {wm}; poly Wp = {wp};")
print("poly E1 = F - 1 - t*G*Wm; poly E2 = F + 1 - (t+1)*H*Wp;")
print("ideal I = coeffs(E1, t), coeffs(E2, t);")
dist = "*".join(f"((a{i} - a{j}) + {L}*(c{i} - c{j}))" for i in range(1, k + 1) for j in range(1, i))
print(f"I = I, 1 - y*{'*'.join(A)}*{dist};")
print(f"ring S = {P}, ({', '.join(vars_)}), dp;")
print("ideal J = imap(R, I);")
print("option(redSB); ideal B = std(J);")
print('if (B[1] == 1) { "UNIT IDEAL: no (complex, hence no real) cubic Williams family in normal form, k = ' + str(k) + ', L = ' + str(L) + '"; }')
print('else { "dim = ", dim(B), " size = ", size(B); B; }')
print("quit;")
