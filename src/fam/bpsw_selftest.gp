default(parisizemax,2*10^9);
read("bpsw.gp");
print("spsp2<1e5 (A001262 starts 2047,3277,4033,4681,8321): ", select(n->!isprime(n)&&sprp2(n), [3..100000]));
print("slprp<1e5 (A217255: 5459,5777,10877,16109,18971,22499,24569,25199,40309,58519,75077,97439): ", select(n->!isprime(n)&&n%2&&slprp(n), [3..100000]));
print("primes failing: ", #select(p->!bpsw(p), primes([3,200000])));
print("composite BPSW < 1e6: ", #select(n->!isprime(n)&&bpsw(n), [3..300000]));
