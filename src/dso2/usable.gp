default(parisizemax, 5*10^9);
read("../fam/bpsw.gp");
zrank(p,D,Q)={my(m=p-kronecker(D,p), f=factor(m)); for(i=1,#f~, my(r=f[i,1]);
  for(e=1,f[i,2], if(lucUV(m/r,Q,D,p)[1]==0, m/=r, break))); m};
oddp(x)=x>>valuation(x,2);
E=lcm(vector(20,i,i)); G=gcd((1<<E)-1, fibonacci(E)); P=factor(G)[,1];
P=select(p->p>5, P); print("primes in G(lcm(1..20)) (p > 5): ", #P, "   max ", vecmax(P));
U=List(); C=Map();
{foreach(P, p, my(o=znorder(Mod(2,p)), z=zrank(p,5,-1), kr=kronecker(5,p), ok);
   ok = (kr==-1) || gcd(oddp(o), oddp(z))==1;
   if(ok, listput(U,[p,kr,valuation(o,2),valuation(z,2),oddp(o),oddp(z)]);
     my(key=[valuation(o,2),valuation(z,2)], c=0); if(mapisdefined(C,key,&c), mapput(C,key,c+1), mapput(C,key,1))));}
U=Vec(U);
print("BPSW-usable (inert, or split with coprime odd orders): ", #U, "  of which inert: ", #select(u->u[2]==-1,U));
print("by 2-adic class [j1,j2]: ", Mat(C));
print("largest usable: ", vecsort(apply(u->u[1],U))[max(1,#U-9)..#U]);
