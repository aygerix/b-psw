\\ Faithful BPSW (problem statement): strong Fermat base 2 + strong Lucas, Selfridge Method A.
sprp2(n)={my(d=n-1,s=valuation(d,2),x); d>>=s; x=Mod(2,n)^d; if(x==1||x==-1,return(1)); for(r=1,s-1, x=x^2; if(x==-1,return(1))); 0};
selD(n)={if(issquare(n),return(0)); my(D=5); while(kronecker(D,n)!=-1, D=if(D>0,-D-2,-D+2)); D};
\\ U_k,V_k mod n for P=1,Q, via binary method (n odd)
lucUV(k,Q,D,n)={my(U=Mod(0,n),V=Mod(2,n),Qk=Mod(1,n),b=binary(k),i2=Mod(2,n)^-1);
  for(i=1,#b, U=U*V; V=V^2-2*Qk; Qk=Qk^2; if(b[i], my(U2=(U+V)*i2, V2=(D*U+V)*i2); U=U2; V=V2; Qk=Qk*Q)); [U,V,Qk]};
slprp(n)={my(D=selD(n)); if(D==0,return(0)); my(Q=(1-D)/4, e=n+1, t=valuation(e,2)); e>>=t;
  my(r=lucUV(e,Q,D,n), U=r[1], V=r[2], Qk=r[3]); if(U==0||V==0,return(1));
  for(i=1,t-1, V=V^2-2*Qk; Qk=Qk^2; if(V==0,return(1))); 0};
bpsw(n)={n>2 && n%2 && sprp2(n) && slprp(n)};
