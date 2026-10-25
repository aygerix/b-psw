read("q22_moments.gp");
{ my(B = bnd(7, e1, e2, lam, mu)); for(i = 1, #B, print("B", i, " vars ", variables(B[i]), " deg p: ", vector(7, a, poldegree(B[i], eval(Str("p", a)))))); write("q22_k7_B.txt", B); }
