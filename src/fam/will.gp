default(parisizemax, 4*10^9);
lpsi(d)={if(d==1,return('x-2)); if(d==2,return('x+2)); my(P=polcyclo(d,'t), m=poldegree(P)/2, Q=0);
  while(poldegree(P)>=m && P!=0, my(k=poldegree(P)-m, c=pollead(P)); Q+=c*'x^k; P-= c*'t^(m-k)*('t^2+1)^k; P=simplify(P)); Q};
N=90; PS=vector(N,d,lpsi(d));
\\ pairwise unit filter for BOTH signs, both directions
U=matrix(N,N); {for(a=1,N, for(b=1,N, if(a!=b, U[a,b]=(abs(polresultant(PS[a]-1,PS[b]))==1 && abs(polresultant(PS[a]+1,PS[b]))==1))));}
G=matrix(N,N,a,b, U[a,b]&&U[b,a]);
print("edges: ", sum(a=1,N,sum(b=a+1,N,G[a,b])));
\\ full Williams test on a set
wt(S)={for(j=1,#S, my(P=1); for(i=1,#S, if(i!=j, P*=PS[S[i]])); if((P-1)%(PS[S[j]]-1)!=0 || (P-1)%(PS[S[j]]+1)!=0, return(0))); 1};
\\ half tests
ft(S)={for(j=1,#S, my(P=1); for(i=1,#S, if(i!=j, P*=PS[S[i]])); if((P-1)%(PS[S[j]]-1)!=0, return(0))); 1};
\\ clique enumeration up to size K
K=6; cnt=vector(K); hits=List();
rec(S,cand)={my(k=#S); if(k>=3, cnt[k]++; if(wt(S), listput(hits,S); print("WILLIAMS ",S)));
  if(k==K, return); for(i=1,#cand, my(c=cand[i], nc=List()); for(t=i+1,#cand, if(G[c,cand[t]], listput(nc,cand[t]))); rec(concat(S,[c]),Vec(nc)))};
{for(a=1,N, my(nc=List()); for(b=a+1,N, if(G[a,b], listput(nc,b))); rec([a],Vec(nc)));}
print("cliques by size: ",cnt); print("Williams sets: ",#hits);
