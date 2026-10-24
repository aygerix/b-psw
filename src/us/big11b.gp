default(parisizemax, 10^9);
{foreach([[-15,86495844719236350721,0,0],[-23,69444174200193754009,2,1]], t,
  my(t0=getabstime(), r=testm(t[1], t[2], t[3], t[4]));
  print("D=",t[1]," m=",t[2]," class (",t[3],",",t[4],"): testm ", r, "  time ", (getabstime()-t0)\1000, " s"))}
\\ D=-7 small splits (A=1,3,O/3... ) checked directly
{my(D=-7,Q=2,m=2560959248555062987,O=(m-1)>>valuation(m-1,2),F=factor(O)); print("D=-7 O=",O," F=",F);
 forvec(s=vector(#F~,i,[0,1]), my(A=prod(i=1,#F~,if(s[i],F[i,1]^F[i,2],1)),B=O/A,W1=2*A,W2=2*B);
   if(min(W1,W2)<1000, my(g); if(W1<W2, my(N1=(1<<W1)-1); g=gcd(N1,lift(lucUV(W2,Q,D,N1)[1])), my(U=Uexact(W2,Q)); g=gcd(U,lift(Mod(2,U)^W1)-1));
     print("  A=",A," W1=",W1," W2=",W2," g=",g," factors ",if(g>1,factor(g),[]))))}
quit
