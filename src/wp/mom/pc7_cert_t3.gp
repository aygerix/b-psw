\\ D-basis vs B-basis: rank of the functionals at a double root (random w's), and data dump in D-basis
read("pc7_cert_shape.gp");
RR = 3/10;
subR(P) = substvec(P, [e1, e2, lam, mu], [2*RR, RR^2, 4*RR/(1 - RR^2)^2, 2*(1 - 3*RR^2)/(1 - RR^2)^2]);
XX = vector(7, j, eval(Str("z", j)));
WS = vector(7, j, (random(1000) + 1)/337);
BS = apply(subR, Bof(XX)); DS = apply(subR, Dof(XX));
JB = matrix(8, 7, i, j, substvec(deriv(BS[i], XX[j]), XX, WS));
JD = matrix(8, 7, i, j, substvec(deriv(DS[i], XX[j]), XX, WS));
print("double root r = 3/10: rank dB/dw = ", matrank(JB), ", rank dD/dw = ", matrank(JD));
BV = vector(8, i, substvec(BS[i], XX, WS));
print("B3 - (B2 - r B4) at random w (0 means B-encoding has only 3 independent V-functionals): ", BV[3] - (BV[2] - RR*BV[4]));
L = [[3/10, 3/10], [1/10, 1/10], [0, 0], [-1/2, -1/2], [29/100, 31/100], [1/10, 3/10], [-3/10, 3/10], [-7/10, 97/100]];
dump("pc7_cert_t3.dat", L);
