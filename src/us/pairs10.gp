read("rows_1e10.gp");
\\ list every class-consistent, cross-coprime subset (size >= 2) of admissible split primes < 1e10, with its product
{foreach([5,-7,-11,13,-15,17,-19,-23,29,-31], D, my(R=select(r->r[1]==D && adm(r), ROWS), out=List());
  my(cl=Set(apply(r->[valuation(r[3],2),valuation(r[4],2)],R)));
  foreach(cl, c, my(L=select(r->[valuation(r[3],2),valuation(r[4],2)]==c, R), k=#L);
    forsubset(k, S, if(#S<2, next); my(A=1,B=1,m=1); for(i=1,#S, A=lcm(A,oddp(L[S[i]][3])); B=lcm(B,oddp(L[S[i]][4])); m*=L[S[i]][2]);
      if(gcd(A,B)==1, listput(out,[m, #S, c]))));
  printf("D=%d: multi-prime split parts %d; > 1e16: %d %s\n", D, #out, #select(t->t[1]>10^16, out), Vec(select(t->t[1]>10^16, out))));}
