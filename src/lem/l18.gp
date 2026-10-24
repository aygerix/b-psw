\\ Lemma 18 check: for p inert (resp. split) w.r.t. D, v2(z_D(p)) = v2(p+1) (resp. v2(p-1)) iff (Q/p) = -1,
\\ where Q = (1-D)/4; and v2(ord_p 2) = v2(p-1) iff (2/p) = -1. Exhaustive over primes p < 2e5, 10 Selfridge D.
read("../fam/bpsw.gp");
zrank(p,D,Q)={my(m=p-kronecker(D,p), f=factor(m)); for(i=1,#f~, my(r=f[i,1]);
  for(e=1,f[i,2], if(lucUV(m/r,Q,D,p)[1]==0, m/=r, break))); m};
{foreach([5,-7,-11,13,-15,17,-19,-23,29,-31], D, my(Q=(1-D)/4, bad=0, nin=0, nsp=0);
  forprime(p=3, 2*10^5, if(D%p==0 || Q%p==0, next); my(k=kronecker(D,p), z=zrank(p,D,Q), o=znorder(Mod(2,p)), full);
    full = (valuation(z,2) == valuation(p-k,2));
    if(full != (kronecker(Q,p)==-1), bad++);
    if((valuation(o,2)==valuation(p-1,2)) != (kronecker(2,p)==-1), bad++);
    if(k==-1, nin++, nsp++));
  printf("D=%d: inert %d split %d, violations %d\n", D, nin, nsp, bad));}
