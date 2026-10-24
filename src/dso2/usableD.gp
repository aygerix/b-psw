\\ Same census for other Selfridge D: G = gcd(2^E - 1, U_E(1,Q)), E = lcm(1..k).
default(parisizemax, 5*10^9);
read("../fam/bpsw.gp");
zrank(p,D,Q)={my(m=p-kronecker(D,p), f=factor(m)); for(i=1,#f~, my(r=f[i,1]);
  for(e=1,f[i,2], if(lucUV(m/r,Q,D,p)[1]==0, m/=r, break))); m};
oddp(x)=x>>valuation(x,2);
Uexact(E,Q)=([1,-Q;1,0]^E)[2,1];
cen(D,k)={my(Q=(1-D)/4, E=lcm(vector(k,i,i)), G=gcd((1<<E)-1, Uexact(E,Q)), P=select(p->D%p && Q%p && p>2, factor(G)[,1]), U=List());
  foreach(P, p, my(o=znorder(Mod(2,p)), z=zrank(p,D,Q), kr=kronecker(D,p));
    if(kr==-1 || gcd(oddp(o),oddp(z))==1, listput(U,[p,kr])));
  U=Vec(U); printf("D=%d k=%d: %d primes (max %d digits), usable %d (inert %d), largest usable %s\n", D, k, #P,
    #digits(vecmax(concat(P,[1]))), #U, #select(u->u[2]==-1,U), if(#U, vecmax(apply(u->u[1],U)), 0))};
foreach([5,-7,-11,13,-15,17], D, foreach([16,20], k, cen(D,k)));
