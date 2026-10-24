\\ Lucas-Wieferich test: for Wieferich primes p (base 2) and Selfridge D, is z_D(p^2) = z_D(p),
\\ i.e. p^2 | U_{z(p)} for P=1, Q=(1-D)/4 ?  (Needed for p^2 | n in a BPSW psp: Lemma 1, row F3.)
read("../fam/bpsw.gp");
zr(p,D,Q)={my(m=p-kronecker(D,p), f=factor(m)); for(i=1,#f~, my(r=f[i,1]);
  for(e=1,f[i,2], if(lucUV(m/r,Q,D,p)[1]==0, m/=r, break))); m};
{foreach([1093,3511], p,
  if(Mod(2,p^2)^(p-1)!=1, print(p," not Wieferich"); next);
  my(hits=List(), D=5, cnt=0);
  while(abs(D)<=20000,
    if(!issquare(D) && D%p, my(Q=(1-D)/4); if(Q%p,
      cnt++; my(z=zr(p,D,Q)); if(lucUV(z,Q,D,p^2)[1]==0, listput(hits,[D,z,kronecker(D,p)]))));
    D=if(D>0,-D-2,-D+2));
  print("p=",p," ord_p(2)=",znorder(Mod(2,p))," D tested ",cnt," Lucas-Wieferich hits [D, z(p), (D/p)]: ",Vec(hits)));}
