\\ Post-process "UNFACTORED part" lines of lastprime.gp runs: recompute G(m) for the reported (D, m),
\\ factor it COMPLETELY, and run bpsw(m q) on every prime q | G(m) not dividing m.
\\ Input: global FIX = [[D, m], ...]. Same conventions as lastprime.gp (class from m's primes, q split iff (D/m) = -1).
default(parisizemax, 4*10^9);
read("../fam/bpsw.gp");
zrank(p,D,Q)={my(m=p-kronecker(D,p), f=factor(m)); for(i=1,#f~, my(r=f[i,1]);
  for(e=1,f[i,2], if(lucUV(m/r,Q,D,p)[1]==0, m/=r, break))); m};
oddp(x)=x>>valuation(x,2);
Uexact(E,Q)=([1,-Q;1,0]^E)[2,1];
fix1(D,m)={my(Q=(1-D)/4, P=factor(m)[,1], j1=-1, j2=-1, kr=kronecker(D,m), E1, E2, G, F, nc=0);
  foreach(P, p, my(o=znorder(Mod(2,p)), z=zrank(p,D,Q)); if(j1<0, j1=valuation(o,2); j2=valuation(z,2),
    if(valuation(o,2)!=j1 || valuation(z,2)!=j2, error("mixed class at m=",m))));
  E1=(1<<j1)*oddp(m-1); E2=(1<<j2)*oddp(if(kr==-1, m+1, m-1));
  G=gcd((1<<E1)-1, Uexact(E2,Q)); F=factor(G);
  for(i=1,#F~, my(q=F[i,1]); if(m%q==0, next); if(!isprime(q), error("factor not prime")); nc++;
    if(bpsw(m*q), print("BPSW HIT D=",D," m=",m," q=",q)));
  printf("D=%d m=%d: log10 G %.1f, %d prime factors, %d candidates tested, largest %d digits\n",
    D, m, log(G)/log(10), #F~, nc, #digits(vecmax(F[,1])))};
foreach(FIX, t, fix1(t[1], t[2]));
