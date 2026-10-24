\\ Positive control for the last-prime divisibilities (each half separately, on real pseudoprimes):
\\  spsp(2) n = m*q:  q | 2^(2^j1 * odd(m-1)) - 1, j1 = v2(ord_q 2)
\\  strong Lucas psp (Selfridge D) n = m*q: q | U_(2^j2 * odd(m+-1)), sign + if q split, - if q inert; j2 = v2(z(q))
read("../fam/bpsw.gp");
zrank(p,D,Q)={my(m=p-kronecker(D,p), f=factor(m)); for(i=1,#f~, my(r=f[i,1]);
  for(e=1,f[i,2], if(lucUV(m/r,Q,D,p)[1]==0, m/=r, break))); m};
oddp(x)=x>>valuation(x,2);
nf=0; bf=0; nl=0; bl=0;
{forstep(n=9,2*10^6,2, if(isprime(n) || issquare(n), next);
  if(sprp2(n) && issquarefree(n), my(f=factor(n)[,1]); nf++;
    foreach(f,q, my(m=n/q, j1=valuation(znorder(Mod(2,q)),2)); if(Mod(2,q)^((1<<j1)*oddp(m-1))!=1, bf++)));
  if(slprp(n) && issquarefree(n), my(D=selD(n), Q=(1-D)/4, f=factor(n)[,1]); nl++;
    foreach(f,q, my(m=n/q, z=zrank(q,D,Q), j2=valuation(z,2), E=(1<<j2)*oddp(if(kronecker(D,q)==1, m+1, m-1)));
      if(lucUV(E,Q,D,q)[1]!=0, bl++))));}
print("spsp(2) squarefree n < 2e6: ",nf,"  violations of the Fermat last-prime divisibility: ",bf);
print("strong Lucas psp squarefree n < 2e6: ",nl,"  violations of the Lucas last-prime divisibility: ",bl);
