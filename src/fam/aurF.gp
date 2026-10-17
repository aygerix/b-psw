default(parisizemax, 4*10^9);
lpsi(d)={if(d==1,return('x-2)); if(d==2,return('x+2)); my(P=polcyclo(d,'t), m=poldegree(P)/2, Q=0);
  while(poldegree(P)>=m && P!=0, my(k=poldegree(P)-m, c=pollead(P)); Q+=c*'x^k; P-= c*'t^(m-k)*('t^2+1)^k; P=simplify(P)); Q};
\\ shapes: irreducible factors (up to sign/content) of a list of polynomials in 'z
shapes(L)={my(S=List()); foreach(L,P, my(F=factor(P)); for(i=1,#F~, my(f=F[i,1]); if(poldegree(f)>=1, f=f/content(f); if(pollead(f)<0,f=-f); if(!setsearch(Set(S),f), listput(S,f))))); Vec(S)};
\\ side s=-1: need f-1 | P-1 (Fermat identity, Lucas-exponential shapes)
\\ side s=+1: need f+1 | P-1 ... i.e. (n+1) divisible by f+1 when n=f*P (Lucas identity, Fermat-exponential shapes)
search(S,s,K,name)={my(N=#S,A=matrix(N,N),E=0,hits=List(),cnt=vector(K));
  for(a=1,N, for(b=1,N, if(a!=b, A[a,b]=(abs(polresultant(S[a]+s,S[b]))==1))));
  my(G=matrix(N,N,a,b,A[a,b]&&A[b,a]));
  my(tst=(T)->for(j=1,#T, my(P=1); for(i=1,#T, if(i!=j, P*=S[T[i]])); if((P-1)%(S[T[j]]+s)!=0, return(0))); 1);
  my(rec=(R,cand,self)->my(k=#R); if(k>=2, cnt[k]++; if(tst(R), listput(hits,R))); if(k<K, for(i=1,#cand, my(c=cand[i], nc=List()); for(t=i+1,#cand, if(G[c,cand[t]], listput(nc,cand[t]))); self(concat(R,[c]),Vec(nc),self))));
  for(a=1,N, my(nc=List()); for(b=a+1,N, if(G[a,b], listput(nc,b))); rec([a],Vec(nc),rec));
  print(name,": shapes ",N,", sym edges ",sum(a=1,N,sum(b=a+1,N,G[a,b])),", cliques ",cnt,", HITS ",#hits);
  foreach(hits,h, print("   ",vector(#h,i,S[h[i]]))); hits};
D=60;


SF=shapes(concat(vector(D,d,subst(polcyclo(d),'x,2*'z^2)), vector(D,d,polcyclo(d,'z))));
hf=search(SF,1,5,"Fermat-exp shapes in X (Lucas identity)");
print("done");
