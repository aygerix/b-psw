\\ Fermat-exponential family: p_j = Phi_{d_j}(X), X = 2^(2^a).
\\ Lucas side needs Phi_{d_j}(X)+1 | prod_i Phi_{d_i}(X)+1 in Q[X]
\\ => necessary: Res(Phi_a+1, Phi_b) = +-1 for all pairs a != b in the set.
N=400;
cnt=0;
for(a=1,N, A=polcyclo(a)+1; L=List(); for(b=1,N, if(b!=a, r=abs(polresultant(A,polcyclo(b))); if(r==1, listput(L,b)))); if(#L>0 && #L<40, print(a," : ",Vec(L))); if(#L>=40, print(a," : ",#L," partners")); cnt+=#L);
print("total unit pairs ",cnt);
