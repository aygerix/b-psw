\\ Last-prime search for BPSW pseudoprimes (Lemma 14 consequence).
\\ For n = m*q (q prime, q not | m) with Selfridge D and 2-adic class (j1, j2) of the primes of m:
\\   ord_q 2 | 2^j1 * odd(m-1)  and  z(q) | 2^j2 * odd(m+1) (q split) or 2^j2 * odd(m-1) (q inert),
\\ where q is split iff (D/m) = -1 (since (D/n) = -1). Hence q | G = gcd(2^E1 - 1, U_E2(1,Q)).
\\ Every BPSW psp with Selfridge D in DLIST and a prime factor q with n/q <= M0 is found (take m = n/q),
\\ unless a gcd is skipped (E > EMAX) or left unfactored; both are counted and printed.
\\ m prime: class of m itself. Usage: gp -q lastprime.gp
default(parisizemax, 4*10^9);
read("../fam/bpsw.gp");
zrank(p,D,Q)={my(m=p-kronecker(D,p), f=factor(m)); for(i=1,#f~, my(r=f[i,1]);
  for(e=1,f[i,2], if(lucUV(m/r,Q,D,p)[1]==0, m/=r, break))); m};
oddp(x)=x>>valuation(x,2);
Uexact(E,Q)=([1,-Q;1,0]^E)[2,1];
M0=default_m0; DLIST=[5,-7,-11,13,-15,17,-19,-23]; EMAX=2^23;
proc(m,j1,j2,kr)={my(E1=(1<<j1)*oddp(m-1), E2=(1<<j2)*oddp(if(kr==-1, m+1, m-1)), G, F);
  if(E1>EMAX || E2>EMAX, gbig++; return);
  G=gcd((1<<E1)-1, Uexact(E2,gQ)); gng++; if(G>gmaxG, gmaxG=G);
  if(G==1, return);
  F=factor(G,10^7);
  for(i=1,#F~, my(q=F[i,1]); if(m%q==0, next);
    if(!ispseudoprime(q), gunf++; print("  UNFACTORED part D=",gD," m=",m," digits ",#digits(q)); next);
    gcand++; if(bpsw(m*q), listput(ghits,[gD,m,q]); print("BPSW HIT D=",gD," m=",m," q=",q)))};
rec(m,i0,j1,j2,A,B,kr)={if(m>1, gnm++; proc(m,j1,j2,kr));
  for(i=i0,#gL, my(p=gL[i], d); if(m*p>M0, break); d=mapget(gdat,p);
    if(m>1 && (d[1]!=j1 || d[2]!=j2), next);
    if(m>1 && (gcd(A,d[4])>1 || gcd(B,d[3])>1 || A%p==0 || B%p==0 || gcd(m, d[3]*d[4])>1), next);
    rec(m*p, i+1, d[1], d[2], lcm(A,d[3]), lcm(B,d[4]), kr*d[5]))};
{foreach(DLIST, D, gD=D; gQ=(1-D)/4; gnm=0; gng=0; gbig=0; gcand=0; gunf=0; gmaxG=1; ghits=List();
  gdat=Map(); foreach(primes([3,M0]), p, if(D%p==0 || gQ%p==0, next); my(o=znorder(Mod(2,p)), z=zrank(p,D,gQ));
     mapput(gdat, p, [valuation(o,2), valuation(z,2), oddp(o), oddp(z), kronecker(D,p)]));
  gL=select(p->mapisdefined(gdat,p), primes([3,M0]));
  rec(1,1,-1,-1,1,1,1);
  printf("D=%d: cofactors m<=%d tested %d, gcds %d, skipped(E>EMAX) %d, unfactored %d, candidate q %d, max log10 G %.1f, hits %d\n",
     D, M0, gnm, gng, gbig, gunf, gcand, log(gmaxG)/log(10), #ghits));}
