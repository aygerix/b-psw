\\ Full identity test on triangles of the unit-resultant graph.
\\ Lucas-exp: p_j = Psi_{d_j}(Y), Y=L_{2^a}: need (Psi_j - 1) | (prod_{i!=j} Psi_i - 1)
\\ Fermat-exp: p_j = Phi_{d_j}(X), X=2^(2^a): need (Phi_j + 1) | (prod_{i!=j} Phi_i - 1)
lpsi(d)={if(d==1,return(x-2)); if(d==2,return(x+2)); my(f=factor(polresultant(polcyclo(d,y), y^2-x*y+1, y))[1,1]); f/pollead(f)};
tst(V,s)={for(j=1,#V, my(P=1); for(i=1,#V, if(i!=j, P*=V[i])); if((P-1)%(V[j]+s)!=0, return(0))); 1};
run(F,s,name)={my(N=#F,A=matrix(N,N),E=0,T=0,H=0);
 for(a=1,N, for(b=1,N, if(a!=b, A[a,b]=(abs(polresultant(F[a]+s,F[b]))==1))));
 for(a=1,N, for(b=a+1,N, if(A[a,b]&&A[b,a], E++;
   for(c=b+1,N, if(A[a,c]&&A[c,a]&&A[b,c]&&A[c,b], T++;
     if(tst([F[a],F[b],F[c]],s), H++; print(name," HIT ",[a,b,c])))))));
 print(name,": sym edges ",E,", triangles ",T,", full-identity triples ",H); A};
N=100;
AL=run(vector(N,d,lpsi(d)),-1,"Lucas");

