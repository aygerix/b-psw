read("pc7_cert_bb.gp");
read("pc7_cert_alpha1.gp");
SHAPE = [-3/10, 4/5]; LO = 6/25; HI = 8/5;
foreach([4, 8, 12], n, N = n; print("== N = ", N); bb());
