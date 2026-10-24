\\ Exhaustive check: every BPSW psp n with Selfridge D* = D whose prime factors all have ord_p 2 | E and
\\ z_D(p) | E (E = lcm(1..k)) is a product of the "usable" primes of G(E) = gcd(2^E - 1, U_E(1,Q)).
\\ Necessary filters (Lemma 1): one 2-adic class (v2 ord, v2 z) for all primes; odd(ord) and odd(z)
\\ coprime across the whole set. Every surviving subset of size >= 2 is tested with the exact bpsw().
default(parisizemax, 5*10^9);
read("../fam/bpsw.gp");
zrank(p,D,Q)={my(m=p-kronecker(D,p), f=factor(m)); for(i=1,#f~, my(r=f[i,1]);
  for(e=1,f[i,2], if(lucUV(m/r,Q,D,p)[1]==0, m/=r, break))); m};
oddp(x)=x>>valuation(x,2);
Uexact(E,Q)=([1,-Q;1,0]^E)[2,1];
run(D,k)={my(Q=(1-D)/4, E=lcm(vector(k,i,i)), G=gcd((1<<E)-1, Uexact(E,Q)), P=select(p->D%p && Q%p && p>2, factor(G)[,1]), U=List(), byc=Map(), tested=0, hits=0);
  foreach(P, p, my(o=znorder(Mod(2,p)), z=zrank(p,D,Q), kr=kronecker(D,p));
    if(kr==-1 || gcd(oddp(o),oddp(z))==1, listput(U,[p,valuation(o,2),valuation(z,2),oddp(o),oddp(z)])));
  foreach(U, u, my(c=[u[2],u[3]], L=List()); if(mapisdefined(byc,c,&L), listput(L,u); mapput(byc,c,L), mapput(byc,c,List([u]))));
  foreach(Mat(byc)[,1]~, c, my(L=Vec(mapget(byc,c)), r=#L);
    forsubset(r, S, if(#S<2, next); my(A=1, B=1, n=1, ok=1);
      for(i=1,#S, my(u=L[S[i]]); A=lcm(A,u[4]); B=lcm(B,u[5]); n*=u[1]);
      if(gcd(A,B)>1, next); tested++;
      if(bpsw(n), hits++; print("BPSW HIT D=",D," n=",n," primes=",apply(i->L[i][1],Vec(S))))));
  printf("D=%d k=%d: usable %d, classes %d, consistent subsets tested %d, BPSW hits %d\n", D, k, #U, #byc, tested, hits)};
foreach([5,-7,-11,13,-15,17,-19,-23], D, run(D,20));
