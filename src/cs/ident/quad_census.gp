\\ src/cs/ident/quad_census.gp -- "Chernick forms with polynomial cofactors" (R6 first experiment, quadratic case).
\\ Members f_j = 1 + t (a_j t + b_j): shared Fermat factor A(t) = t, linear cofactors g_j = a_j t + b_j (a_j >= 1).
\\ Korselt identity f_j - 1 | F - 1 in Q[t] is checked exactly. Lucas side: h_j = gcd(F + 1, f_j + 1) in Q[t],
\\ defect delta_j = 2 - deg h_j. delta_j = 2 ("dead"): z(f_j(t)) | C_j (Bezout constant), finitely many t [PROVED, advisor].
\\ delta_j = 1: open (needs index >= |linear cofactor|). Williams-good: delta_j = 0. Theorem 10: not all delta_j = 0.
\\ Then an exact end-to-end scan of the best families (fewest dead members, least total defect).
read("univ_lib.gp");
census(k,A,B,Tscan)={
  my(mem=List()); for(a=1,A, for(b=-B,B, listput(mem,[a,b])));
  my(M=#mem, pols=vector(M,i, 1+'t*(mem[i][1]*'t+mem[i][2])), res=List(), nK=0, hist=matrix(2*k+1,k+1));
  forsubset([M,k], v,
    my(F=prod(i=1,k,pols[v[i]]), ok=1);
    for(i=1,k, if((F-1)%(pols[v[i]]-1)!=0, ok=0; break));
    if(!ok, next); nK++;
    my(def=vector(k,i, 2-poldegree(gcd(F+1, pols[v[i]]+1))), td=vecsum(def), nd=#select(x->x==2,def));
    hist[td+1,nd+1]++;
    listput(res,[vector(k,i,mem[v[i]]), def, td, nd]));
  print("k=",k," a<=",A," |b|<=",B,": Korselt quadratic families = ",nK);
  print("  histogram rows = total Lucas defect 0..",2*k,", columns = #dead members 0..",k,":"); print(hist);
  \\ scan every family with the minimum number of dead members
  my(best=vecmin(vector(#res,i,res[i][4])));
  print("  minimum #dead members = ",best);
  my(nsc=0);
  for(r=1,#res, if(res[r][4]>best, next);
    my(fam=res[r][1], fs=vector(k,i,1+'t*(fam[i][1]*'t+fam[i][2])), np=0, nf=0, nl=0, nb=0);
    nsc++; if(nsc>40, break);
    for(tt=-Tscan,Tscan, my(P=vector(k,i,subst(fs[i],'t,tt)), ok=1);
      if(vecmin(P)<3 || #Set(P)<k, next);
      for(i=1,k, if(!ispseudoprime(P[i]), ok=0; break)); if(!ok,next);
      np++; my(n=prod(i=1,k,P[i])); my(f2=sprp2(n), fl=slprp(n)); nf+=f2; nl+=fl;
      if(f2&&fl, nb++; print("   BPSW HIT fam=",fam," t=",tt," n=",n," factors=",P)));
    print("  fam ",fam," defects ",res[r][2],": prime tuples |t|<=",Tscan," : ",np,"  sprp2 ",nf,"  slprp ",nl,"  BPSW ",nb));
  res};
default(parisizemax, 1500000000);
R3=census(3,10,15,20000);
R4=census(4,5,7,20000);
