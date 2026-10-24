\\ Check of Lemma 12(c): for odd composite non-square n, if the Selfridge D* is NOT a prime discriminant
\\ then either D* = -15, or gcd(n, Q*) > 1 (strong Lucas fails), or 5*7 | n (strong base-2 fails: v2(ord_5 2)=2,
\\ v2(ord_7 2)=0). Scan all odd n <= NMAX; count selected D* by type and any violation.
read("../fam/bpsw.gp");
isfund(D)={my(a=abs(D)); isprime(a) && D%4==1};      \\ prime discriminant l* (l >= 5)
NMAX=3*10^6; cnt=Map(); bad=0; c15=0; tot=0;
{forstep(n=9,NMAX,2, if(issquare(n) || isprime(n), next); my(D=selD(n)); if(D==0, next); tot++;
  if(isfund(D), next);
  my(Q=(1-D)/4, why="");
  if(D==-15, c15++; next);
  if(gcd(n,Q)>1, why="gcdQ"); if(why=="" && n%35==0, why="35");
  my(c=0); if(mapisdefined(cnt,[D,why],&c), mapput(cnt,[D,why],c+1), mapput(cnt,[D,why],1));
  if(why=="", bad++; if(bad<10, print("VIOLATION n=",n," D=",D))));}
print("odd composites scanned: ",tot,"  D*=-15: ",c15,"  violations: ",bad);
print("non-prime-discriminant D* (other than -15) by [D, reason]: ", Mat(cnt));
