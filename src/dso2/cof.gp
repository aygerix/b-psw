default(parisizemax, 5*10^9);
read("../fam/bpsw.gp");
zrank(p,D,Q)={my(m=p-kronecker(D,p), f=factor(m)); for(i=1,#f~, my(r=f[i,1]);
  for(e=1,f[i,2], if(lucUV(m/r,Q,D,p)[1]==0, m/=r, break))); m};
cof(k)={my(E=lcm(vector(k,i,i)), G=gcd((1<<E)-1, fibonacci(E)), F=factor(G,10^7), c=select(q->q>10^7,F[,1]));
  if(#c==0, return); c=c[1]; printf("k=%d: cofactor %d digits; factoring...\n", k, #digits(c));
  my(t=getabstime(), f=factor(c)); printf("  factored in %d s: %s\n", (getabstime()-t)\1000, f[,1]~);
  foreach(f[,1], p, my(o=znorder(Mod(2,p)), z=zrank(p,5,-1), kr=kronecker(5,p));
    printf("  p=%d (%d digits) (5/p)=%d ord=%s z=%s idxF=%s idxL=%s\n", p, #digits(p), kr, factor(o), factor(z),
      factor((p-1)/o), factor((p-kr)/z)))};
cof(18);
alarm(3000, cof(20));
