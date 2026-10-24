\\ DSO census for highly composite exponents: G(E) = gcd(2^E - 1, F_E) (D = 5) contains EVERY prime p with
\\ ord_p 2 | E and z(p) | E. Attempt D measured ord, z <= 1500 only; here E = lcm(1..k) reaches 2.3e8 (smooth orders).
\\ For each E: factor G completely (trial division to 1e7, then ECM/MPQS via factor()), report the largest primes,
\\ and for each prime its (ord, z, Fermat index, Lucas index, (5/p)).
default(parisizemax, 5*10^9);
read("../fam/bpsw.gp");
zrank(p,D,Q)={my(m=p-kronecker(D,p), f=factor(m)); for(i=1,#f~, my(r=f[i,1]);
  for(e=1,f[i,2], if(lucUV(m/r,Q,D,p)[1]==0, m/=r, break))); m};
run(k)={my(E=lcm(vector(k,i,i)), t=getabstime(), G, F, P);
  G=gcd((1<<E)-1, fibonacci(E)); F=factor(G, 10^7); P=F[,1];
  printf("k=%d E=%d: log10 G = %.1f, %d prime-power factors after trial division, gcd+td time %d s\n",
    k, E, log(G)/log(10), #P, (getabstime()-t)\1000);
  my(big=select(q->q>10^7, P)); if(#big, print("  cofactors > 1e7: ", apply(q->[#digits(q), isprime(q)], big)));
  my(L=vecsort(select(q->isprime(q) && q>5, P),,4)); \\ descending
  for(i=1,min(12,#L), my(p=L[i], o=znorder(Mod(2,p)), z=zrank(p,5,-1));
    printf("  p=%d (5/p)=%d ord=%d z=%d  idxF=%d idxL=%d\n", p, kronecker(5,p), o, z, (p-1)/o, (p-kronecker(5,p))/z));
  #select(q->isprime(q), P)};
foreach([12,14,16,18,20], k, run(k));
