default(parisizemax,10^9)
\\ Last inert prime for a split cofactor (Lemma 14 + the split-prime lemma).
\\ If n = m*q is a BPSW psp whose ONLY inert prime is q, then m is a product of admissible split primes, all in
\\ class (j1,j2), and q satisfies odd(ord_q 2)*odd(z_q) | O = odd(m-1), v2 ord_q 2 = j1, v2 z_q = j2.
\\ So q | g_A = gcd(2^(2^j1 A) - 1, U_(2^j2 B)) for some unitary divisor A of O (B = O/A).
\\ Each g_A is computed modulo the SMALLER of the two numbers (sizes about 2^j sqrt(O) at worst).
read("../fam/bpsw.gp");
oddp(x)=x>>valuation(x,2);
Uexact(E,Q)=([1,-Q;1,0]^E)[2,1];
lgalpha(D,Q)=if(D>0, log((1+sqrt(D))/2)/log(2), log(sqrt(abs(Q)))/log(2));
\\ returns [list of prime factors q of the union of g_A, number of A tried, unfactored cofactors]
lastq(D,m,j1,j2)={my(Q=(1-D)/4, O=oddp(m-1), F=factor(O), r=#F~, qs=List(), unf=List(), la=lgalpha(D,Q));
  forvec(s=vector(r,i,[0,1]), my(A=prod(i=1,r,if(s[i],F[i,1]^F[i,2],1)), B=O/A, W1=(1<<j1)*A, W2=(1<<j2)*B, g);
    if(W1 <= W2*la,
      my(N1=(1<<W1)-1); g=gcd(N1, lift(lucUV(W2,Q,D,N1)[1])),
      my(U=Uexact(W2,Q)); g=gcd(U, lift(Mod(2,U)^W1)-1));
    if(g>1, my(G=factor(g,10^9)); for(i=1,#G~, my(q=G[i,1]); if(q>10^18 && !ispseudoprime(q), listput(unf,q); next);
      if(q>2 && m%q && D%q && Q%q, listput(qs,q)))));
  [Set(qs), 1<<r, Vec(unf)]};
\\ test all candidates; return [#candidates, hits]
testm(D,m,j1,j2)={my(R=lastq(D,m,j1,j2), hits=List());
  foreach(R[1], q, if(bpsw(m*q), listput(hits,q); print("BPSW HIT D=",D," m=",m," q=",q)));
  if(#R[3], print("  UNFACTORED D=",D," m=",m," parts=",R[3]));
  [#R[1], Vec(hits), R[2]]};
