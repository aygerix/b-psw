\\ audC_vert.gp -- exact steps of the vertical all-B subcase (INT(t,v)).
f(A,C) = 2*x + 1 + x*(x+1)*(A*x + C);
nrm(phi, q) = polresultant(q, phi)/polcoeff(q, 2)^poldegree(phi);
\\ (1) Mv_t at a member with the same a: N_h(f_(a,cb)) - 1 = -2 k R + 2 k^2 V N',  k = cb - ct
ht = (f(a,ct) + 1)/(x+1); V = 1/a; u = ct/a; R = -2*V - u*(1-u)/2; Np = 1 - u + 2*V; k = cb - ct;
print("(1) ", simplify(nrm(f(a,cb), ht) - 1 - (-2*k*R + 2*k^2*V*Np)));
print("(2) -R/(V N') - (4a + ct(a-ct))/(2(a-ct+2)) = ", simplify(-R/(V*Np) - (4*a + ct*(a-ct))/(2*(a - ct + 2))));
Ft = (4*a + c*(a-c))/(2*(a - c + 2));
print("(3) numerator of Ft - (c - 4): ", simplify(numerator(Ft - (c - 4))*1), "  [expect -(c^2 - (a+12)c + 12a + 16) up to sign/const]");
print("    discriminant in c: ", factor(poldisc(c^2 - (a+12)*c + 12*a + 16, c)));
\\ (4) case c_lo = a - S, S = sqrt(8a), a = S^2/8, d = a - c
e4 = subst(subst(numerator(Ft - (c - (a - S))), c, a - d), a, S^2/8);
print("(4) numerator (c_lo = a - S) in d: ", simplify(e4), "  vs d^2 + (S^2/8 - 2S + 4) d + S^2/2 - 4S times const: ratio ", simplify(e4/(d^2 + (S^2/8 - 2*S + 4)*d + S^2/2 - 4*S)));
Dl = (S^2/8 - 2*S + 4)^2 - 4*(S^2/2 - 4*S);
print("    Delta(S) = ", Dl, " ; real roots: ", polroots(Dl*1.0));
\\ sign of Delta on [4(1+sqrt3), 16]: Sturm count of roots in the interval (use rational enclosure 10.92 < 4(1+sqrt3) < 10.93)
print("    roots of Delta in [10.92, 16]: ", polsturm(Dl, [1092/100, 16]), " ; Delta(12) = ", subst(Dl, S, 12));
\\ (5) breakpoints: sqrt(8a) = a - 4 at a = 8 + 4 sqrt3 ; check (a - 4)^2 - 8a = a^2 - 16 a + 16
print("(5) a^2 - 16a + 16 roots: ", polroots(a^2 - 16*a + 16), "; S at a = 8+4sqrt3: ", sqrt(8*(8 + 4*sqrt(3))), " = 4(1+sqrt3) = ", 4*(1+sqrt(3)));
