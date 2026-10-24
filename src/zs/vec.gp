\\ Unit / log-derivative data of Psi_d at the special points theta_e (roots of Psi_e), e in EA.
\\ Korselt identity for a set S of Psi_d  <=>  for every special point e of S (root of some Psi_d - 1, d in S):
\\   prod_{d in S} Psi_d(theta_e) = 1  and, if some member has multiplicity >= 2 at e,
\\   sum_{d in S} Psi_d'(theta_e)/Psi_d(theta_e) = 0   (multiplicity 3+ excluded here).
default(parisizemax, 2*10^9);
Psi(d)={my(f); if(d==1,return('Y-2)); if(d==2,return('Y+2));
  f=polresultant(polcyclo(d,'w), 'w^2-'Y*'w+1, 'w); f=factor(f)[1,1]; f/pollead(f)};
EA=[1,2,3,4,6,5,10,8,12,7,14,9,18]; NMAX=400;
BNF=vector(#EA,i, if(poldegree(Psi(EA[i]))>=2, bnfinit(subst(Psi(EA[i]),'Y,'x),1), 0));
EM=Map(); {for(e=1,2*NMAX+2, mapput(EM,Psi(e),e));}
{for(d=3,NMAX, my(P=Psi(d), F=factor(P-1), sp=List(), ok=1, out="");
  for(i=1,#F~, my(g=F[i,1]/pollead(F[i,1]), e=0); if(mapisdefined(EM,g,&e) && setsearch(Set(EA),e), listput(sp,[e,F[i,2]]), ok=0));
  if(!ok, next);
  out=Str(d,"|",poldegree(P),"|"); for(i=1,#sp, out=Str(out,sp[i][1],":",sp[i][2],if(i<#sp,",","")));
  for(i=1,#EA, my(e=EA[i], g=Psi(e), r=P%g, dr);
    out=Str(out,"|",e,"=");
    if(r==0, out=Str(out,"Z"); next);
    dr=(deriv(P)*Mod(1,g))/(P*Mod(1,g));
    if(poldegree(g)==1, my(t=-polcoeff(g,0), v=subst(P,'Y,t), dv=subst(deriv(P),'Y,t));
        if(abs(v)!=1, out=Str(out,"N"), out=Str(out,"U",sign(v),";",dv/v)),
      my(b=BNF[i], u=bnfisunit(b, subst(lift(Mod(P,g)),'Y,'x)));
        if(#u==0, out=Str(out,"N"),
          out=Str(out,"U",apply(lift,u~),";",Vecrev(lift(dr),poldegree(g))))));
  print(out));}
