\\ Shallue-Webster's two remaining admissible split primes for D = 5 below 2^40 (beyond our 1e10 scan):
\\ verify orders, then run the 2-factor last-prime search.
zrank(p,D,Q)={my(mm=p-kronecker(D,p), f=factor(mm)); for(i=1,#f~, my(r=f[i,1]);
  for(e=1,f[i,2], if(lucUV(mm/r,Q,D,p)[1]==0, mm/=r, break))); mm};
{foreach([26509131221, 422013019339], p, my(o=znorder(Mod(2,p)), z=zrank(p,5,-1), j1=valuation(o,2), j2=valuation(z,2));
  printf("p=%d prime=%d (5/p)=%d ord=%d z=%d class=(%d,%d) gcd(odd)=%d\n", p, isprime(p), kronecker(5,p), o, z, j1, j2, gcd(oddp(o),oddp(z)));
  print("  testm: ", testm(5, p, j1, j2)));}
