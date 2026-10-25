\\ src/cs/ident/cubic_census.gp -- degree-3 escape hatch for the extremal principle (Theorem 10 Steps 2-3 need degree 2).
\\ Members f_j = 1 + A(t) (a_j t + b_j): shared quadratic Fermat factor A, linear cofactors (NOT composed-linear).
\\ Exact Korselt check f_j - 1 | F - 1 in Q[t]; Lucas defect via h_j = gcd(F+1, f_j+1); dead = deg h_j = 0.
\\ Reports families with NO dead member (the only kind that can give infinitely many BPSW psp), then an exact
\\ end-to-end scan (bpsw.gp) of the families with the fewest dead members.
read("univ_lib.gp");
default(parisizemax, 1500000000);
ccensus(Apol,k,Amax,Bmax,Tscan)={
  my(mem=List()); for(a=1,Amax, for(b=-Bmax,Bmax, listput(mem,[a,b])));
  my(M=#mem, pols=vector(M,i, 1+Apol*(mem[i][1]*'t+mem[i][2])), res=List(), nK=0, hist=vector(k+1));
  forsubset([M,k], v,
    my(F=prod(i=1,k,pols[v[i]]), ok=1);
    for(i=1,k, if((F-1)%(pols[v[i]]-1)!=0, ok=0; break));
    if(!ok, next); nK++;
    my(def=vector(k,i, 3-poldegree(gcd(F+1, pols[v[i]]+1))), nd=#select(x->x==3,def));
    hist[nd+1]++; listput(res,[vector(k,i,mem[v[i]]), def, nd]));
  my(best=if(#res, vecmin(vector(#res,i,res[i][3])), -1));
  print("A=",Apol," k=",k," a<=",Amax," |b|<=",Bmax,": Korselt families=",nK,"  #dead histogram (0..k)=",hist,"  min #dead=",best);
  my(nsc=0);
  for(r=1,#res, if(res[r][3]>best, next); nsc++; if(nsc>12, break);
    my(fam=res[r][1], fs=vector(k,i,1+Apol*(fam[i][1]*'t+fam[i][2])), np=0, nf=0, nl=0, nb=0);
    for(tt=-Tscan,Tscan, my(P=vector(k,i,subst(fs[i],'t,tt)), ok=1);
      if(vecmin(P)<3 || #Set(P)<k, next);
      for(i=1,k, if(!ispseudoprime(P[i]), ok=0; break)); if(!ok,next);
      np++; my(n=prod(i=1,k,P[i]), f2=sprp2(n), fl=slprp(n)); nf+=f2; nl+=fl;
      if(f2&&fl, nb++; print("   BPSW HIT fam=",fam," t=",tt," n=",n," factors=",P)));
    print("   fam ",fam," Lucas defects ",res[r][2],": prime tuples |t|<=",Tscan,": ",np,"  sprp2 ",nf,"  slprp ",nl,"  BPSW ",nb));
  res};
{
foreach(['t^2+1, 't^2+'t+1, 't^2+2, 't^2+3, 't^2-'t+1, 't^2+'t+2, 't*('t+1), 't*('t+2), 't^2-2], A,
   ccensus(A,3,8,10,5000));
}
