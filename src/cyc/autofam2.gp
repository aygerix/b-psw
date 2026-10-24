\\ One-side-automatic families against the faithful BPSW (fam/bpsw.gp), over primes p in [P0, P1].
\\ Set P0, P1 (and DOFERMAT) before reading this file. Run from src/.
\\  M_p = 2^p - 1 and W_p = (2^p + 1)/3 are spsp(2) automatically, so only the strong Lucas test decides.
\\  Fib(p), p == +-2 mod 5, is a strong Lucas psp for Selfridge D = 5 automatically, so only spsp(2) decides.
\\ Unlike autofam.gp there is NO ispseudoprime filter: a composite that fooled PARI's own BPSW would have been
\\ skipped there. Every n passing both tests is printed and labelled by whether p is a known prime/PRP index.
read("fam/bpsw.gp");
MERS = [2,3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423,9689,9941,11213,19937];
WAG = [3,5,7,11,13,17,19,23,31,43,61,79,101,127,167,191,199,313,347,701,1709,2617,3539,5807,10501,10691,11279,12391,14479];
FIBP = [3,4,5,7,11,13,17,23,29,43,47,83,131,137,359,431,433,449,509,569,571,2971,4723,5387,9311,9677,14431];
both(name,p,known)={print("BOTH ",name,"_",p, if(known," (known prime/PRP index)"," *** NOT a known prime index: HIT candidate ***"))};
{
my(c=vector(3), l=vector(3), h=vector(3), t0=getabstime());
forprime(p=max(3,P0), P1,
  my(n=2^p-1, k=setsearch(MERS,p)>0);
  if(!k, c[1]++); if(slprp(n), if(!k, l[1]++); if(sprp2(n), if(!k, h[1]++); both("M",p,k)));
  n=(2^p+1)/3; k=setsearch(WAG,p)>0;
  if(!k, c[2]++); if(slprp(n), if(!k, l[2]++); if(sprp2(n), if(!k, h[2]++); both("W",p,k)));
  if(p%5==2 || p%5==3,
    n=fibonacci(p); k=setsearch(FIBP,p)>0;
    if(!k, c[3]++); if(sprp2(n), if(!k, l[3]++); if(slprp(n), if(!k, h[3]++); both("Fib",p,k)))));
printf("range [%d,%d], %d s\n", P0, P1, (getabstime()-t0)\1000);
printf("Mersenne composites %d: strong Lucas %d, BPSW %d\n", c[1], l[1], h[1]);
printf("Wagstaff composites %d: strong Lucas %d, BPSW %d\n", c[2], l[2], h[2]);
printf("Fib(p), p=+-2 mod 5, composites %d: spsp2 %d, BPSW %d\n", c[3], l[3], h[3]);
if(DOFERMAT, for(m=5,17, my(n=2^(2^m)+1); printf("F_%d: spsp2 %d strong Lucas %d\n", m, sprp2(n), slprp(n))));
}
quit
