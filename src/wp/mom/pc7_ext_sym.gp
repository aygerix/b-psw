\\ extremal (bottom) member conditions for (2,2) pencils, symbolic in c, e1, e2
Dl = 1 - e1^2 + 2*e2 + e2^2;
lam = 2*e1/Dl; mu = 2*(1 - e1^2 + e2)/Dl;
twoRe(s, q) = lam*(s^3 - 3*q*s) + mu*(s^2 - 2*q) - lam*s - 2*mu;
sm = (c*e1 - lam)/c; qm = (c*e2 + mu)/c;
sp = -lam/c; qp = (mu + lam*e1 - c)/c;
RX = twoRe(sm, qm); RY = twoRe(sp, qp) + 4;
dm = sm^2 - 4*qm; dp = sp^2 - 4*qp;
fa(P) = factor(numerator(P));
print("RX num: ", fa(RX)); print("RX den: ", factor(denominator(RX)));
print("RY num: ", fa(RY)); print("RY den: ", factor(denominator(RY)));
print("dm num: ", fa(dm)); print("dm den: ", factor(denominator(dm)));
print("dp num: ", fa(dp)); print("dp den: ", factor(denominator(dp)));
