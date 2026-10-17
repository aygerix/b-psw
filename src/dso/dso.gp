default(parisizemax, 4*10^9);
\\ Doubly-small-order (DSO) primes: p | gcd(2^a-1, F_b), a,b <= A.
A=1500;
M=vector(A,a,2^a-1); F=vector(A,b,fibonacci(b));
P=Map();
addp(g)={my(f=factor(g,2^40)[,1]); for(i=1,#f, my(p=f[i]); if(p>5 && ispseudoprime(p) && !mapisdefined(P,p), mapput(P,p,1)))};
rowa(a)={for(b=3,A, my(g=gcd(M[a],F[b])); if(g>1, addp(g)))};
for(a=1,A, rowa(a); if(a%300==0, print("a=",a," primes so far ",#P)));
zent(p)={my(m=p-kronecker(5,p)); fordiv(m,d, if(fibonacci(d)%p==0, return(d)))};
L=vecsort(apply(p->[p,znorder(Mod(2,p)),zent(p),kronecker(5,p)], Vec(P)),1);
print("DSO primes (ord,z <= ",A,"): ",#L);
write("dso_primes.txt", L);
for(i=max(1,#L-25),#L, print(L[i]));
