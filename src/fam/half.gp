default(parisizemax, 4*10^9);
read("bpsw.gp");
lpsi(d)={if(d==1,return('x-2)); if(d==2,return('x+2)); my(P=polcyclo(d,'t), m=poldegree(P)/2, Q=0);
  while(poldegree(P)>=m && P!=0, my(k=poldegree(P)-m, c=pollead(P)); Q+=c*'x^k; P-= c*'t^(m-k)*('t^2+1)^k; P=simplify(P)); Q};
N=150; PS=vector(N,d,lpsi(d));
H=matrix(N,N); {for(a=3,N, for(b=3,N, if(a!=b, H[a,b]=(abs(polresultant(PS[a]-1,PS[b]))==1))));}
G=matrix(N,N,a,b, H[a,b]&&H[b,a]);
ft(S)={for(j=1,#S, my(P=1); for(i=1,#S, if(i!=j, P*=PS[S[i]])); if((P-1)%(PS[S[j]]-1)!=0, return(0))); 1};
K=5; fams=List(); cnt=vector(K);
rec(S,cand)={my(k=#S); if(k>=3, cnt[k]++; if(ft(S), listput(fams,S)));
  if(k==K, return); for(i=1,#cand, my(c=cand[i], nc=List()); for(t=i+1,#cand, if(G[c,cand[t]], listput(nc,cand[t]))); rec(concat(S,[c]),Vec(nc)))};
{for(a=3,N, my(nc=List()); for(b=a+1,N, if(G[a,b], listput(nc,b))); rec([a],Vec(nc)));}
print("cliques: ",cnt); print("half-families (Fermat identity): ", Vec(fams));
\\ direct BPSW test on n(Y), Y = L_M and Y = -L_M, M <= 400
{foreach(fams,S, my(hit=0); for(M=1,400, foreach([1,-1],sg, my(Y=sg*(fibonacci(M-1)+fibonacci(M+1)), n=prod(i=1,#S,subst(PS[S[i]],'x,Y)));
   if(n>3 && n%2 && bpsw(n), hit++; print("BPSW PASS S=",S," M=",M," sign=",sg," isprime(n)=",isprime(n))))));
   print(S," tested, passes=",hit));}
