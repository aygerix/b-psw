read("bpsw.gp");
F1=['t^2-'t-1, 't^3-'t^2-2*'t+1, 't^3-3*'t-1];
F2=concat(F1, ['t^4-3*'t^2-2*'t+1]);
\\ direct BPSW test of n(Y) for the Fermat-identity families, all integer Y in range
scan(F,lo,hi)={my(c=0,h=0); for(Y=lo,hi, my(v=vector(#F,i,subst(F[i],'t,Y))); if(vecmin(v)<3, next);
   my(n=prod(i=1,#v,v[i])); if(n%2==0, next);
   if(sprp2(n), c++; if(slprp(n), h++; print("BPSW PASS Y=",Y," factors ",v)))); print("F size ",#F,": Y in [",lo,",",hi,"]: strong-Fermat passes ",c,", BPSW passes ",h)};
scan(F1,-200000,200000);
scan(F2,-200000,200000);
