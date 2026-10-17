default(parisizemax, 4*10^9);
lpsi(d)={if(d==1,return('x-2)); if(d==2,return('x+2)); my(P=polcyclo(d,'t), m=poldegree(P)/2, Q=0);
  while(poldegree(P)>=m && P!=0, my(k=poldegree(P)-m, c=pollead(P)); Q+=c*'x^k; P-= c*'t^(m-k)*('t^2+1)^k; P=simplify(P)); Q};
shapes(L)={my(S=List()); foreach(L,P, my(F=factor(P)); for(i=1,#F~, my(f=F[i,1]); if(poldegree(f)>=1, f=f/content(f); if(pollead(f)<0,f=-f); if(!setsearch(Set(S),f), listput(S,f))))); Vec(S)};
okres(r,l)={if(r==0,return(0)); r=abs(r); my(f=factor(abs(l))[,1]); for(i=1,#f, while(r%f[i]==0, r/=f[i])); r==1};
tst(T)={for(j=1,#T, my(P=1); for(i=1,#T, if(i!=j, P*=SS[T[i]])); if((P-1)%(SS[T[j]]+SG)!=0, return(0))); 1};
rec(R,cand)={my(k=#R); if(k>=2, CNT[k]++; if(tst(R), listput(HITS,R))); if(k<KK, for(i=1,#cand, my(c=cand[i], nc=List()); for(t=i+1,#cand, if(GG[c,cand[t]], listput(nc,cand[t]))); rec(concat(R,[c]),Vec(nc))))};
search(S,s,K,name)={SS=S; SG=s; KK=K; CNT=vector(K); HITS=List(); my(N=#S,A=matrix(N,N));
  for(a=1,N, for(b=1,N, if(a!=b, A[a,b]=okres(polresultant(S[a]+s,S[b]),pollead(S[a])))));
  GG=matrix(N,N,a,b,A[a,b]&&A[b,a]);
  for(a=1,N, my(nc=List()); for(b=a+1,N, if(GG[a,b], listput(nc,b))); rec([a],Vec(nc)));
  print(name,": shapes ",N,", sym edges ",sum(a=1,N,sum(b=a+1,N,GG[a,b])),", cliques ",CNT,", HITS ",#HITS);
  foreach(HITS,h, print("   ",vector(#h,i,SS[h[i]])))};
D=60;
\\ Lucas-exp: k odd (Y=5Z^2-2) and k even (Y=5Z^2+2) kept SEPARATE
search(shapes(vector(D,d,subst(lpsi(d),'x,5*'z^2-2))),-1,5,"Lucas-exp, k odd");
search(shapes(vector(D,d,subst(lpsi(d),'x,5*'z^2+2))),-1,5,"Lucas-exp, k even");
search(shapes(concat(vector(D,d,subst(polcyclo(d),'x,2*'z^2)), vector(D,d,polcyclo(d,'z)))),1,5,"Fermat-exp (X=2^m, 2X^2=2^(2m+1))");
print("done");
