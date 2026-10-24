\\ Twisted Lucas-forced Korselt families (D = 5, Y = L_M with M odd).
\\ P_d(Y) = i^(-phi(d)/2) Psi_d(iY) for 4 | d: at Y = L_M, iY = u + 1/u with u = i phi^M, so any prime
\\ p | P_d(L_M) has (i phi^M)^d = 1, phi^(Md) = 1, z(p) | M d (forced). Korselt: P_a - 1 | prod P - 1.
\\ Necessary pair condition: Res(P_a - 1, P_b) = +-1 for every a != b in the family (both directions).
\\ Usage: gp -q twist.gp  (DMAX, KMAX below); also re-runs the untwisted Psi case as a positive control.
default(parisizemax, 2*10^9);
Psi(d)={my(f); if(d==1,return('Y-2)); if(d==2,return('Y+2));
  f=polresultant(polcyclo(d,'w), 'w^2-'Y*'w+1, 'w); f=factor(f)[1,1]; f/pollead(f)};
Tw(d)={my(g=subst(Psi(d),'Y,I*'Y)*I^(-poldegree(Psi(d)))); g=simplify(g); if(type(g)!="t_POL" || variable(g)!='Y, error("bad"));
  for(k=0,poldegree(g), if(imag(polcoeff(g,k))!=0, error("non-real P_",d))); real(g)};
ft(S)={for(j=1,#S, my(P=1); for(i=1,#S, if(i!=j, P*=gPS[S[i]])); if((P-1)%(gPS[S[j]]-1)!=0, return(0))); 1};
rec(S,cand)={my(k=#S); if(k>=2, gcnt[k]++; if(ft(S), listput(gF,vector(#S,i,gDL[S[i]])))); if(k==gK, return);
  for(i=1,#cand, my(c=cand[i], nc=List()); for(t=i+1,#cand, if(gG[c,cand[t]]&&gG[cand[t],c], listput(nc,cand[t]))); rec(concat(S,[c]),Vec(nc)))};
fam(PS, DL, KMAX, tag)={my(n=#DL); gPS=PS; gDL=DL; gK=KMAX; gG=matrix(n,n); gcnt=vector(KMAX); gF=List();
  for(a=1,n, for(b=1,n, if(a!=b, gG[a,b]=(abs(polresultant(PS[a]-1,PS[b]))==1))));
  for(a=1,n, my(nc=List()); for(b=a+1,n, if(gG[a,b]&&gG[b,a], listput(nc,b))); rec([a],Vec(nc)));
  print(tag,": candidates ",n,"  cliques by size ",gcnt,"  Korselt families: ",Vec(gF))};
DMAX=600; KMAX=8;
DL=select(d->d%4==0, [1..DMAX]); PT=vector(#DL,i,Tw(DL[i]));
print("P_d samples: ", vector(4,i,[DL[i],PT[i]]));
fam(PT, DL, KMAX, "twisted 4|d (5|d allowed), d<=600");
DL2=select(d->d>=3, [1..60]); fam(vector(#DL2,i,Psi(DL2[i])), DL2, 4, "control: untwisted Psi, d<=60");
