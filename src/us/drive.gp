\\ Driver: every BPSW psp with Selfridge D* = D and EXACTLY ONE inert prime factor q is n = m*q, where m is a
\\ product of admissible split primes (odd(ord) coprime to odd(z), min(v2 ord, v2 z) <= 1), all in one class,
\\ with odd(ord_p) coprime to odd(z_p') across the set. Enumerate such m (primes <= X from the scan, m <= MMAX)
\\ and test every candidate q from lastq(). Needs gq.gp loaded first; ROWS = [[D,p,ord,z],...].
oddp(x)=x>>valuation(x,2);
adm(r)=min(valuation(r[3],2),valuation(r[4],2))<=1;
runD(D, MMAX, COSTMAX)={my(R=select(r->r[1]==D && adm(r), ROWS), cls=Set(apply(r->[valuation(r[3],2),valuation(r[4],2)],R)),
    nm=0, ncand=0, nskip=0, skipped=List(), hits=List());
  foreach(cls, c, my(L=vecsort(select(r->[valuation(r[3],2),valuation(r[4],2)]==c, R), 2), k=#L);
    \\ DFS over subsets with product <= MMAX and cross-coprimality
    my(stack=List([[0,1,1,1]]));  \\ [last index, m, A=lcm odd ord, B=lcm odd z]
    while(#stack, my(t=stack[#stack]); listpop(stack);
      for(i=t[1]+1, k, my(p=L[i][2], m=t[2]*p, A=lcm(t[3],oddp(L[i][3])), B=lcm(t[4],oddp(L[i][4])));
        if(m>MMAX, break); if(gcd(A,B)>1, next);
        listput(stack,[i,m,A,B]);
        my(O=oddp(m-1), w=omega(O), cost=(2^w)*(1<<max(c[1],c[2]))*sqrtint(O));
        if(cost>COSTMAX, nskip++; listput(skipped,[m,cost]); next);
        my(res=testm(D,m,c[1],c[2])); nm++; ncand+=res[1]; foreach(res[2],q,listput(hits,[m,q])))));
  printf("D=%d: admissible split primes %d, classes %d, cofactors m tested %d, candidate q %d, skipped (cost) %d, hits %d\n",
    D, #R, #cls, nm, ncand, nskip, #hits);
  if(#skipped, print("  skipped: ", Vec(skipped)[1..min(5,#skipped)]));
  if(#hits, print("  HITS: ", Vec(hits)))};
