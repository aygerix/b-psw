\\ C2: structured families that pass one test automatically; run the exact BPSW (bpsw.gp) on each composite member.
\\  M_q = 2^q - 1 (q prime, composite): spsp(2) automatically.  W_q = (2^q + 1)/3 (q prime, composite): spsp(2).
\\  F_k = 2^(2^k) + 1 (k >= 5): spsp(2).
\\  L_(2^k) composite: the inert part is a strong Lucas psp for D = 5 when n == -1 mod 2^(k+1); here we test L_(2^k)
\\  itself and all products of an odd number of its known prime factors, for k <= 9.
default(parisizemax, 4*10^9);
read("../fam/bpsw.gp");
QMAX = 10000;
{my(cnt=0, hits=0); forprime(q=3, QMAX, my(n=(1<<q)-1); if(ispseudoprime(n), next); cnt++;
   if(bpsw(n), hits++; print("HIT Mersenne q=",q)));
 print("Mersenne composites q <= ",QMAX,": tested ",cnt,", BPSW passes ",hits);}
{my(cnt=0, hits=0); forprime(q=5, QMAX, my(n=((1<<q)+1)/3); if(ispseudoprime(n), next); cnt++;
   if(bpsw(n), hits++; print("HIT Wagstaff q=",q)));
 print("Wagstaff composites q <= ",QMAX,": tested ",cnt,", BPSW passes ",hits);}
{my(hits=0); for(k=5,16, my(n=(1<<(1<<k))+1); if(bpsw(n), hits++; print("HIT Fermat k=",k)));
 print("Fermat numbers F_5..F_16: BPSW passes ",hits);}
{for(k=3,9, my(L=fibonacci(2^k-1)+fibonacci(2^k+1), f=factor(L,10^9)[,1], P=select(p->ispseudoprime(p),f), hit=0, t=0);
   forsubset(#P, s, if(#s<2, next); my(n=prod(i=1,#s,P[s[i]])); t++; if(bpsw(n), hit++; print("HIT Lucas k=",k," n=",n)));
   printf("L_(2^%d): %d prime factors found (%s), subsets tested %d, passes %d\n", k, #P, if(#P==#f,"complete","partial"), t, hit));}
