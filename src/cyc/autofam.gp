\\ Canonical one-side-automatic families, checked against the faithful BPSW (fam/bpsw.gp).
\\  M_p = 2^p-1 (composite, p prime) and W_p = (2^p+1)/3 are spsp(2) automatically; F_m = 2^(2^m)+1 too.
\\  Fib(p), p == +-2 mod 5 prime, composite: strong Lucas psp for Selfridge D = 5 automatically
\\  (z(q) = p for every q | Fib(p), and p | odd(Fib(p)+1) because Fib(p) == -1 mod p).
read("fam/bpsw.gp");
chk(name,p,n)={my(a=sprp2(n),b=slprp(n)); if(a&&b,print("HIT ",name,"_",p)); [a,b]};
P=20000;
{
my(c=vector(3),s=vector(3),l=vector(3),h=vector(3),r);
forprime(p=3,P,
  n=2^p-1;
  if(!ispseudoprime(n), c[1]++; r=chk("M",p,n); s[1]+=r[1]; l[1]+=r[2]; h[1]+=r[1]*r[2]; if(r[2],print("M_",p," strong Lucas")));
  n=(2^p+1)/3;
  if(!ispseudoprime(n), c[2]++; r=chk("W",p,n); s[2]+=r[1]; l[2]+=r[2]; h[2]+=r[1]*r[2]; if(r[2],print("W_",p," strong Lucas")));
  if(p%5==2||p%5==3,
    n=fibonacci(p);
    if(!ispseudoprime(n), c[3]++; r=chk("Fib",p,n); s[3]+=r[1]; l[3]+=r[2]; h[3]+=r[1]*r[2]; if(r[1],print("Fib_",p," spsp(2)")))));
printf("Mersenne composites %d: spsp2 %d, strong Lucas %d, BPSW %d\n",c[1],s[1],l[1],h[1]);
printf("Wagstaff composites %d: spsp2 %d, strong Lucas %d, BPSW %d\n",c[2],s[2],l[2],h[2]);
printf("Fib(p), p=+-2 mod 5, composites %d: spsp2 %d, strong Lucas %d, BPSW %d\n",c[3],s[3],l[3],h[3]);
for(m=5,17, n=2^(2^m)+1; printf("F_%d: spsp2 %d strong Lucas %d\n",m,sprp2(n),slprp(n)));
}
