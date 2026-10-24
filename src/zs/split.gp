\\ For each d, factor Psi_d(Y) - 1 over Q and identify factors as Psi_e (min poly of 2cos(2pi/e)).
\\ Output: d, deg Psi_d, list [e, multiplicity] or "other(deg)" for non-cyclotomic factors.
default(parisizemax, 2*10^9);
Psi(d)={my(f); if(d==1,return('Y-2)); if(d==2,return('Y+2));
  f=polresultant(polcyclo(d,'w), 'w^2-'Y*'w+1, 'w); f=factor(f)[1,1]; f/pollead(f)};
N=400; EM=Map(); PS=vector(2*N+2,e,Psi(e)); for(e=1,#PS, mapput(EM,PS[e],e));
{for(d=3,N, my(F=factor(PS[d]-1), L=List(), full=1, me=0);
  for(i=1,#F~, my(g=F[i,1]/pollead(F[i,1]), e=0);
    if(mapisdefined(EM,g,&e), listput(L,[e,F[i,2]]); me=max(me,e), full=0; listput(L,["x",poldegree(g),F[i,2]])));
  if(full && me<=12, print(d," deg ",poldegree(PS[d])," : ",Vec(L))));}
print("--- d<=",N," with Psi_d-1 fully cyclotomic, any e: count");
{my(c=0,cs=0); for(d=3,N, my(F=factor(PS[d]-1), full=1); for(i=1,#F~, my(g=F[i,1]/pollead(F[i,1]), e=0); if(!mapisdefined(EM,g,&e), full=0)); c+=full); print(c," of ",N-2);}
