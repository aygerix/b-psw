\\ src/cs/ident/mixed_census.gp -- shared-factor Korselt families f_j = 1 + A(t) g_j(t), g_j = a_j t + b_j with a_j >= 0
\\ (constant cofactors allowed: mixed degrees). The earlier quad/cubic censuses (a_j >= 1) found ONLY composed-linear
\\ families (all g_j proportional, i.e. f_j = 1 + c_j u with u = A (t - r)), which are Korselt over Q trivially and are
\\ covered by the sign lemma. Here we flag NON-composed families (Lemma 5 is one: A = (t-2)(t+1), g = 1, t, t+1).
read("univ_lib.gp");
default(parisizemax, 1500000000);
composed(fam)={my(k=#fam); for(i=2,k, if(fam[i][1]*fam[1][2]-fam[i][2]*fam[1][1]!=0, return(0))); 1};
mcensus(Apol,k,Amax,Bmax,Tscan)={
  my(mem=List()); for(a=0,Amax, for(b=-Bmax,Bmax, if(a==0&&b<=0,next); listput(mem,[a,b])));
  my(M=#mem, pols=vector(M,i, 1+Apol*(mem[i][1]*'t+mem[i][2])), nK=0, nNC=0, hist=vector(k+1), out=List());
  forsubset([M,k], v,
    my(fam=vector(k,i,mem[v[i]]));
    if(composed(fam), next);
    my(F=prod(i=1,k,pols[v[i]]), ok=1);
    for(i=1,k, if((F-1)%(pols[v[i]]-1)!=0, ok=0; break));
    if(!ok, next); nNC++;
    my(def=vector(k,i, poldegree(pols[v[i]])-poldegree(gcd(F+1, pols[v[i]]+1))),
       dead=vector(k,i, poldegree(gcd(F+1, pols[v[i]]+1))==0), nd=vecsum(dead));
    hist[nd+1]++; listput(out,[fam,def,nd]));
  print("A=",Apol," k=",k," 0<=a<=",Amax," |b|<=",Bmax,": NON-composed Korselt families=",nNC,"  #dead histogram (0..k)=",hist);
  for(r=1,min(#out,15), my(fam=out[r][1], fs=vector(k,i,1+Apol*(fam[i][1]*'t+fam[i][2])), np=0, nf=0, nl=0, nb=0);
    for(tt=-Tscan,Tscan, my(P=vector(k,i,subst(fs[i],'t,tt)), ok=1);
      if(vecmin(P)<3 || #Set(P)<k, next);
      for(i=1,k, if(!ispseudoprime(P[i]), ok=0; break)); if(!ok,next);
      np++; my(n=prod(i=1,k,P[i]), f2=sprp2(n), fl=slprp(n)); nf+=f2; nl+=fl;
      if(f2&&fl, nb++; print("   BPSW HIT fam=",fam," t=",tt," n=",n," factors=",P)));
    print("   fam ",fam," degrees ",vector(k,i,poldegree(fs[i]))," Lucas defects ",out[r][2]," dead ",out[r][3],": prime tuples |t|<=",Tscan,": ",np,"  sprp2 ",nf,"  slprp ",nl,"  BPSW ",nb));
  out};
{
foreach(['t^2-'t-2, 't^2+1, 't^2+'t+1, 't^2-2, 't*('t+1), 't*('t-3), 't^2-'t-1, 't^2+2], A,
   mcensus(A,3,5,7,20000));
}
