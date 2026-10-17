default(parisizemax, 4*10^9);
\\ Psi_d = min poly of 2cos(2pi/d), via x^m Psi(x+1/x) = Phi_d(x) (palindromic reduction)
lpsi(d)={if(d==1,return('x-2)); if(d==2,return('x+2)); my(P=polcyclo(d,'t), m=poldegree(P)/2, Q=0, T);
  while(poldegree(P)>=m && P!=0, my(k=poldegree(P)-m, c=pollead(P)); Q+=c*'x^k; P-= c*'t^(m-k)*('t^2+1)^k; P=simplify(P)); Q};
N=300; PS=vector(N,d,lpsi(d));
chk={my(ok=1); for(d=3,60, my(Y=2*cos(2*Pi/d)); if(abs(subst(PS[d],'x,Y))>1e-20, ok=0)); ok};
print("Psi sanity: ", chk);
E=900; PE=vector(E,d,lpsi(d));
ide(f)={my(g=f/pollead(f)); for(e=1,E, if(poldegree(PE[e])==poldegree(g) && PE[e]==g, return(e))); 0};
{for(d=1,N, foreach([-1,1],s, my(F=factor(PS[d]+s), L=List(), ok=1);
   for(i=1,#F~, my(e=ide(F[i,1])); if(e==0, ok=0; break); listput(L,[e,F[i,2]]));
   if(ok, print("d=",d," Psi",if(s<0,"-1","+1")," = ",Vec(L)))));}
