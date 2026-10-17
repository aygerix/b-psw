default(parisizemax, 4*10^9);
\\ All monic integer polys of degree 2..4 with small coeffs, irreducible, positive for large t.
L=List();
forvec(v=vector(2,i,[-8,8]), my(f='t^2+v[1]*'t+v[2]); if(polisirreducible(f), listput(L,f)));
forvec(v=vector(3,i,[-4,4]), my(f='t^3+v[1]*'t^2+v[2]*'t+v[3]); if(polisirreducible(f), listput(L,f)));
forvec(v=vector(4,i,[-3,3]), my(f='t^4+v[1]*'t^3+v[2]*'t^2+v[3]*'t+v[4]); if(polisirreducible(f), listput(L,f)));
S=Vec(L); N=#S; print("polys: ",N);
\\ unit-resultant graph, both signs, both directions
Gm=matrix(N,N); Gp=matrix(N,N);
{for(a=1,N, my(Am=S[a]-1, Ap=S[a]+1); for(b=1,N, if(a!=b,
   Gm[a,b]=(abs(polresultant(Am,S[b]))==1); Gp[a,b]=(abs(polresultant(Ap,S[b]))==1))));}
GW=matrix(N,N,a,b, Gm[a,b]&&Gm[b,a]&&Gp[a,b]&&Gp[b,a]);
GF=matrix(N,N,a,b, Gm[a,b]&&Gm[b,a]);
print("Williams edges ",sum(a=1,N,sum(b=a+1,N,GW[a,b])),"  Fermat-half edges ",sum(a=1,N,sum(b=a+1,N,GF[a,b])));
tst(T,full)={my(D=sum(i=1,#T,poldegree(S[T[i]]))); for(j=1,#T, if(3*poldegree(S[T[j]])>D, return(0)));
  for(j=1,#T, my(P=1); for(i=1,#T, if(i!=j, P*=S[T[i]])); if((P-1)%(S[T[j]]-1)!=0, return(0)); if(full && (P-1)%(S[T[j]]+1)!=0, return(0))); 1};
rec(R,cand,G,full)={my(k=#R); if(k>=3, CNT[k]++; if(tst(R,full), listput(HITS,R)));
  if(k<6, for(i=1,#cand, my(c=cand[i], nc=List()); for(t=i+1,#cand, if(G[c,cand[t]], listput(nc,cand[t]))); rec(concat(R,[c]),Vec(nc),G,full)))};
run(G,full,name)={CNT=vector(6); HITS=List(); for(a=1,N, my(nc=List()); for(b=a+1,N, if(G[a,b], listput(nc,b))); rec([a],Vec(nc),G,full));
  print(name," cliques ",CNT," hits ",#HITS); for(i=1,min(#HITS,30), print("  ",vector(#HITS[i],j,S[HITS[i][j]])))};
run(GW,1,"WILLIAMS");
run(GF,0,"FERMAT-HALF");
print("done");
