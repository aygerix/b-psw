\\ src/cs/ident/univ_lib.gp -- universal-form reductions ("Theorem U") and positive controls.
\\ Three one-parameter linear-type universal forms, members indexed by a set S of positive integers d:
\\   F-form (Chernick):        p_d = 1 + d m,        m = lcm(S) t, lcm(S) | sum(S)  -> p_d - 1 | n - 1 (Korselt)
\\   L-form (Lucas-Carmichael): p_d = d m - 1,       m = lcm(S) t, lcm(S) | sum(S), k odd -> p_d + 1 | n + 1
\\   T-form (Chebyshev T_2):    p_d = 2 (d m)^2 - 1,  m = 2 lcm(S) t, lcm(S)^2 | sum d^2, k odd -> p_d + 1 | n + 1
\\     (the T-form is the L-form in u = m^2 with cofactors 2 d^2; recorded as a check, not a new family)
\\ The other side reduces to a FIXED integer (independent of m):
\\   F-form: odd z(p_d)    | n + 1  <=>  odd z(p_d)    | NF(S,d) = d^(k-1) - prod_{e!=d} (d - 2e)
\\   L-form: odd ord_p(2)  | n - 1  <=>  odd ord_p(2)  | NL(S,d) = prod_{e!=d} (2e - d) - d^(k-1)
\\   T-form: odd ord_p(2)  | n - 1  <=>  odd ord_p(2)  | NT(S,d) = prod_{e!=d} (2e^2 - d^2) - d^(2(k-1))
read("../../fam/bpsw.gp");
oddp(x)= if(x==0,0,x>>valuation(x,2));
lucU(k,Q,D,p)= lucUV(k,Q,D,p)[1];
\\ rank of apparition of U(1,Q) mod p (D = 1 - 4Q); 0 if p | D Q
zrank(p,D)={my(Q=(1-D)/4, kr=kronecker(D,p)); if(kr==0 || Q%p==0, return(0));
  my(e=p-kr, f=factor(e));
  for(i=1,#f~, my(q=f[i,1]); for(j=1,f[i,2], if(lucU(e/q,Q,D,p)==0, e/=q, break))); e};
NF(S,i)={my(d=S[i],k=#S); d^(k-1) - prod(j=1,k, if(j==i,1, d-2*S[j]))};
NL(S,i)={my(d=S[i],k=#S); prod(j=1,k, if(j==i,1, 2*S[j]-d)) - d^(k-1)};
NT(S,i)={my(d=S[i],k=#S); prod(j=1,k, if(j==i,1, 2*S[j]^2-d^2)) - d^(2*(k-1))};
member(form,d,m)= if(form=="F", 1+d*m, if(form=="L", d*m-1, 2*(d*m)^2-1));
Nval(form,S,i)= if(form=="F", NF(S,i), if(form=="L", NL(S,i), NT(S,i)));
\\ v2 class of a prime for the selected D: [v2(ord_p 2), v2(z(p))]
cls(p,D)=[valuation(znorder(Mod(2,p)),2), valuation(zrank(p,D),2)];

\\ run(form, Dmax, k, Tmax, maxtuples): enumerate S (gcd 1; additive condition where needed), t <= Tmax.
run(form,Dmax,k,Tmax,maxtuples)={
  my(nS=0,nT=0,nKor=0,nW=0,nF=0,nL=0,nB=0,agree=0,dis=0,inert=0,split=0,defsum=0.,logsum=0.,nbad=0,ngood=0,
     badhist=vector(k+1),minmax=0,allinert=0,onecls=0,ctlF=List(),ctlL=List());
  if(form!="F" && k%2==0, error("L/T forms need k odd"));
  forsubset([Dmax,k], v,
    my(S=Vec(v), L=lcm(S), s=vecsum(S));
    if(gcd(S)!=1, next);
    if(form!="T" && s%L, next);
    if(form=="T" && sum(i=1,k,S[i]^2)%(L^2), next);
    nS++;
    my(bad=vector(k,i, Nval(form,S,i)!=0), nb=vecsum(bad));
    badhist[nb+1]++; if(bad[1] && bad[k], minmax++);
    for(t=1,Tmax,
      my(m=if(form=="T",2*L*t,L*t), P=vector(k,i,member(form,S[i],m)), ok=1);
      if(P[1]<3, next);
      for(i=1,k, if(!ispseudoprime(P[i]), ok=0; break));
      if(!ok, next);
      nT++;
      my(n=prod(i=1,k,P[i]));
      if(form=="F", my(c=1); for(i=1,k, if((n-1)%(P[i]-1), c=0)); nKor+=c,
                    my(c=1); for(i=1,k, if((n+1)%(P[i]+1), c=0)); nW+=c);
      my(f2=sprp2(n), fl=slprp(n)); nF+=f2; nL+=fl;
      if(f2 && #ctlF<3, listput(ctlF,[S,t,n]));
      if(fl && #ctlL<3, listput(ctlL,[S,t,n]));
      if(f2&&fl, nB++; print("BPSW HIT ",form," S=",S," t=",t," n=",n));
      my(D=selD(n), ai=1, cl=0, oc=1);
      for(i=1,k, my(p=P[i], kr=kronecker(D,p));
        if(kr!=-1, split++; ai=0; next);
        inert++;
        my(c=cls(p,D)); if(cl===0, cl=c, if(c!=cl, oc=0));
        my(o = if(form=="F", oddp(zrank(p,D)), oddp(znorder(Mod(2,p)))));
        if(o==0, next);
        my(N=Nval(form,S,i), A=if(form=="F", (n+1)%o==0, (n-1)%o==0), B=(N%o==0));
        if(A==B, agree++, dis++; print("DISAGREE ",form," S=",S," t=",t," i=",i));
        if(N!=0, nbad++; defsum += log(o/gcd(o,if(form=="F",n+1,n-1)))/log(2); logsum += log(p)/log(2), ngood++));
      allinert+=ai; if(ai, onecls+=oc);
      if(nT>=maxtuples, break(2))));
  print("[",form,"] Dmax=",Dmax," k=",k," Tmax=",Tmax,": valid S=",nS," prime tuples=",nT,
        if(form=="F"," Korselt=",", p+1|n+1 for all="), if(form=="F",nKor,nW),
        " sprp2=",nF," slprp=",nL," BPSW=",nB);
  print("   bad-member histogram over S (entry j+1 = #S with j members having N_d!=0): ",badhist,"  S with min AND max bad: ",minmax);
  print("   tuples all-inert for their D: ",allinert," of which one 2-adic class: ",onecls,"; inert members=",inert," split=",split);
  print("   Theorem U equivalence on inert members: agree=",agree," disagree=",dis);
  print("   inert members with N_d!=0: ",nbad,"  mean residual (defect) bits=",if(nbad,defsum/nbad,0)," mean log2 p=",if(nbad,logsum/nbad,0),
        "  defect/log2p=",if(logsum,defsum/logsum,0),";  inert members with N_d=0: ",ngood);
  print("   controls (first passes): sprp2: ",ctlF,"  slprp: ",ctlL);
  [nT,nF,nL,nB]
};
