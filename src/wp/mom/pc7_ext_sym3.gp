Dl = 1 - e1^2 + 2*e2 + e2^2;
lam = 2*e1/Dl; mu = 2*(1 - e1^2 + e2)/Dl;
twoRe(s, q) = lam*(s^3 - 3*q*s) + mu*(s^2 - 2*q) - lam*s - 2*mu;
sm = (c*e1 - lam)/c; qm = (c*e2 + mu)/c;
sp = -lam/c; qp = (mu + lam*e1 - c)/c;
RX = twoRe(sm, qm); RY = twoRe(sp, qp) + 4;
dm = sm^2 - 4*qm; dp = sp^2 - 4*qp;
\\ c = z/Dl
zz(P) = subst(P, c, z/Dl);
X1 = zz(RX) * z^3 / Dl; Y1 = zz(RY) * z^3 / Dl;
print("z^3 RX/Dl = ", simplify(X1));
print("z^3 RY/Dl = ", simplify(Y1));
print("z^2 dm = ", simplify(zz(dm)*z^2));
print("z^2 dp = ", simplify(zz(dp)*z^2));
print("(X1+Y1)/2 = ", simplify((X1+Y1)/2), "   (X1-Y1)/2 = ", simplify((X1-Y1)/2));
