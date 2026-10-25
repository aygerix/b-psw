"""Singular input: the "-" HALF of the real cubic Williams system in Laurent form (session 8, conjecture H-).
f_j = 2t + 1 + t(t+1)(a_j t + c_j), g_j = a_j t^2 + (a_j + c_j) t + (c_j + 2).  Equations: all t-coefficients of
F - 1 - t prod g_j W_-, W_- = t^(k-1) + sum_{i <= k-4} w_i t^i (shape forced, see NOTEBOOK).  Saturation:
prod a_j != 0, prod (c_j + 2) != 0 (real members have g_j(0) = f_j'(0) > 0), pairwise distinct via (a_i - a_j) + L (c_i - c_j).
A unit ideal (for two values of L) proves: no real cubic family in normal form satisfies even the "-" half, for this k.
Usage: python3 gen_cubic_half.py k prime L > file.sing
"""
import sys
k, P, L = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
A = [f"a{j}" for j in range(1, k + 1)]; C = [f"c{j}" for j in range(1, k + 1)]
W = [f"w{i}" for i in range(k - 3)]
vars_ = A + C + W + ["y"]
print(f"ring R = {P}, (t, {', '.join(vars_)}), (dp(1), dp({len(vars_)}));")
print("poly F = 1; poly G = 1;")
for j in range(1, k + 1):
    print(f"F = F*(2*t + 1 + t*(t+1)*(a{j}*t + c{j})); G = G*(a{j}*t^2 + (a{j} + c{j})*t + (c{j} + 2));")
wm = f"t^{k-1}" + "".join(f" + w{i}*t^{i}" for i in range(k - 3))
print(f"poly Wm = {wm};")
print("ideal I = coeffs(F - 1 - t*G*Wm, t);")
dist = "*".join(f"((a{i} - a{j}) + {L}*(c{i} - c{j}))" for i in range(1, k + 1) for j in range(1, i))
cp2 = "*".join(f"(c{j} + 2)" for j in range(1, k + 1))
print(f"I = I, 1 - y*{'*'.join(A)}*{cp2}*{dist};")
print(f"ring S = {P}, ({', '.join(vars_)}), dp;")
print("ideal J = imap(R, I);")
print("option(redSB); ideal B = std(J);")
print('if (B[1] == 1) { "UNIT IDEAL: the minus half alone has no solution, k = ' + str(k) + ', L = ' + str(L) + '"; }')
print('else { "dim = ", dim(B), " vdim = ", vdim(B), " size = ", size(B); }')
print("quit;")
