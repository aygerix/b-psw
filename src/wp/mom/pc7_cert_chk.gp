\\ check Bof (3-body form) against bnd() of q22_moments.gp at random rational data
read("q22_moments.gp");
read("pc7_cert_dec.gp");
{
  my(k = 7, B = bnd(k, e1, e2, lam, mu), ok = 1);
  for(trial = 1, 3,
    my(E1 = random(200)/97 - 1, E2 = random(200)/89 - 1, LA = random(200)/31 - 3, MU = random(200)/37 - 2);
    my(ws = vector(k, j, (random(1000) + 1)/113));
    my(Bn = vector(8, i, my(Q = substvec(B[i], [e1, e2, lam, mu], [E1, E2, LA, MU]));
       for(a = 1, 2*k, Q = subst(Q, eval(Str("p", a)), sum(j = 1, k, ws[j]^a))); Q));
    my(Bd = substvec(Bof(ws), [e1, e2, lam, mu], [E1, E2, LA, MU]));
    print("trial ", trial, ": max |diff| = ", vecmax(abs(Bn - Bd)), "  (|B| ~ ", vecmax(abs(Bn)) * 1., ")");
    if(Bn != Bd, ok = 0));
  print(if(ok, "AGREE EXACTLY", "MISMATCH"));
}
