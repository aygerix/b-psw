\\ Chernick-type pool at FIXED m: members p_d = 1 + d m (d <= Dmax), all-inert for D = 5, and the Lucas
\\ residual group they generate (odd z, the side NOT automated).  Sharing lemma predicts
\\ log2 lcm(odd z) ~ sum_d (log2 z_d - O(log2 Dmax)), i.e. N/B_L ~ 1/log2 m.  Output: poolNB.out
read("univ_lib.gp");
default(parisizemax, 1000000000);
doit(m,Dmax)={my(N=0,Lz=1,sumlog=0.,sumz=0.,maxshare=0);
  for(d=1,Dmax, my(p=1+d*m); if(!ispseudoprime(p) || kronecker(5,p)!=-1, next);
     my(z=oddp(zrank(p,5))); N++; my(g=gcd(Lz,z)); maxshare=max(maxshare,g); Lz=lcm(Lz,z); sumlog+=log(p)/log(2); sumz+=log(z)/log(2));
  my(B=log(Lz)/log(2));
  print("m=",m," (log2 m=",round(log(m)/log(2)),") Dmax=",Dmax,": pool N=",N," mean log2 p=",if(N,sumlog/N,0),
        " sum log2 z=",round(sumz)," B_L=log2 lcm(odd z)=",round(B)," N/B_L=",if(B,N/B,0)," largest gcd(new z, lcm so far)=",maxshare)};
\\ m must be prime to 5 (else every member is split for D = 5)
doit(2^4*3^2*7*11*13,20000);
doit(2^6*3^3*7^2*11*13*17*19,20000);
doit(2^10*3^5*7^2*11*13*17*19*23*29,20000);
doit(lcm(vector(40,i,i))/5^2,20000);
