\\ AUDIT: doubling chain. Solve R_j(a) = prod_{i != j}(1 - 2a_i/a_j) = 1, j = 2..k-1, a_1 = 1, a_k = M.
\\ Method (own code): fixed-point iteration of T_j(a) = 2a_{j-1}/(1 - 1/R'_j) at huge M, then continuation in M down to the
\\ target with analytic-Jacobian Newton. Reports #good, R_1, R_k, max |a_j/2^(j-1) - 1| over j = 2..k-1.
default(realprecision, 160);
Rj(s, j) = prod(i=1, #s, if(i==j, 1, 1 - 2*s[i]/s[j]));
Rp(s, j) = prod(i=1, #s, if(i==j || i==j-1, 1, 1 - 2*s[i]/s[j]));
full(x, M) = concat([1.], concat(x, [M]));
newton(x, M) = {
  my(k = #x + 2, s, F, J);
  for(t = 1, 100,
    s = full(x, M);
    F = vector(k-2, r, Rj(s, r+1) - 1)~;
    if (normlp(F) < 1e-70, return([x, 1]));
    J = matrix(k-2, k-2, r, c, my(j = r+1, cc = c+1, R = Rj(s, j));
          if (cc != j, R*(-2/s[j])/(1 - 2*s[cc]/s[j]),
                       R*sum(i=1, k, if(i==j, 0, (2*s[i]/s[j]^2)/(1 - 2*s[i]/s[j]))))); 
    x = x - (matsolve(J, F))~);
  [x, 0];
}
chain(k, Mt) = {
  my(M = max(Mt, 1e30), x = vector(k-2, j, 2.^j), s, ok);
  for(t = 1, 200, s = full(x, M); x = vector(k-2, r, my(j = r+1); 2*s[j-1]/(1 - 1/Rp(s, j))));
  [x, ok] = newton(x, M);
  while (M > Mt, M = max(Mt, M/1.5); [x, ok] = newton(x, M); if (!ok, return([0, M])));
  [x, ok];
}
{
  foreach([3,5,7,9,11,13,15,21], k,
    foreach(concat([2^(k+3), 10^5], if(k >= 15, [10^9], [])), M,
      my(r = chain(k, M));
      if (r[1] == 0, printf("k=%2d M=%-11s  continuation FAILED at M=%.4g\n", k, M, r[2]); next);
      my(s = full(r[1], M), R = vector(k, j, Rj(s, j)), good = sum(j=1, k, abs(R[j]-1) < 1e-60));
      my(dev = vecmax(vector(k-2, j, abs(s[j+1]/2^j - 1))), incr = prod(j=1, k-1, s[j] < s[j+1]));
      printf("k=%2d M=%-11s #good=%2d (middle ones: %d)  R1-1=%.3g  Rk-1=%.3g  max|a_j/2^(j-1)-1| (j=2..k-1)=%.3g  incr=%d\n",
             k, M, good, prod(j=2,k-1, abs(R[j]-1) < 1e-60), R[1]-1, R[k]-1, dev, incr)));
}
