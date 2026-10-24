\\ Both sides: special points of Psi_d - 1 (sign -1 side, n(theta)=+1) and Psi_d + 1 (n(theta)=-1).
\\ Line: d|deg|minus specials e:m,...|plus specials e:m,...|e=STATUS ... for e in EA
\\ STATUS: N (non-unit), Z (zero), or U<unit exps incl. torsion>;<log-derivative coords>
default(parisizemax, 4*10^9);
Psi(d)={my(f); if(d==1,return('Y-2)); if(d==2,return('Y+2));
  f=polresultant(polcyclo(d,'w), 'w^2-'Y*'w+1, 'w); f=factor(f)[1,1]; f/pollead(f)};
NMAX=400; DEGMAX=4;
EA=select(e->(e<=2 || eulerphi(e)/2<=DEGMAX), [1..4*NMAX]);
EA=select(e->e<=120, EA);
PS=vector(120,e,Psi(e)); EM=Map(); for(e=1,#PS, mapput(EM,PS[e],e));
BNF=vector(#EA,i, if(poldegree(PS[EA[i]])>=2, bnfinit(subst(PS[EA[i]],'Y,'x),1), 0));
spec(P)={my(F=factor(P), L=List()); for(i=1,#F~, my(g=F[i,1]/pollead(F[i,1]), e=0);
  if(mapisdefined(EM,g,&e) && setsearch(Set(EA),e), listput(L,[e,F[i,2]]), return(0))); Vec(L)};
fmt(L)={my(s=""); for(i=1,#L, s=Str(s,L[i][1],":",L[i][2],if(i<#L,",",""))); s};
{for(d=3,NMAX, my(P=Psi(d), sm=spec(P-1), sp=spec(P+1), out);
  if(sm==0 && sp==0, next);
  out=Str(d,"|",poldegree(P),"|",if(sm==0,"X",fmt(sm)),"|",if(sp==0,"X",fmt(sp)));
  for(i=1,#EA, my(e=EA[i], g=PS[e], r=P%g, dr, u);
    out=Str(out,"|",e,"=");
    if(r==0, out=Str(out,"Z"); next);
    dr=(deriv(P)*Mod(1,g))/(P*Mod(1,g));
    if(poldegree(g)==1, my(t=-polcoeff(g,0), v=subst(P,'Y,t), dv=subst(deriv(P),'Y,t));
        out=Str(out, if(abs(v)!=1, "N", Str("U",sign(v),";",dv/v))),
      u=bnfisunit(BNF[i], subst(lift(Mod(P,g)),'Y,'x));
      out=Str(out, if(#u==0, "N", Str("U",apply(lift,u~),";",Vecrev(lift(dr),poldegree(g)))))));
  print(out));}
