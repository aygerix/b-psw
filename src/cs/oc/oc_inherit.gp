\\ Line OC: (a) Fibonacci primitive parts Phi*_f = prod_{d|f} F_d^mu(f/d) / gcd(., f) that are PRP, f composite <= FMAX,
\\ analysed like the other families (Lucas side forced: z = f); (b) inheritance test: does p -+ 1 of a cyclotomic
\\ primitive-part prime (p = Phi_e(2)/gcd, e composite) or of a Phi*_f prime carry more 1500-smooth content than a
\\ random prime of the same size in the same residue class (p == 1 mod e resp. p == -+1 mod f)?  Raw smooth bits,
\\ no index absorption.  Control: KC random primes per member.
read("oc_fam.gp");
Dl = [5,-7,-11,13,-15,17,-19,-23];
smbits(N) = { my(N2 = N >> valuation(N, 2), S = smoothpart(N2, 1500)[1]); log(S)/log(2); }
phistar(f) = { my(P = prod(i = 1, #divisors(f), my(d = divisors(f)[i]); fibonacci(d)^moebius(f/d))); P / gcd(P, f); }
if(type(FMAX) != "t_INT", FMAX = 1500);
if(type(EMAX) != "t_INT", EMAX = 1500);
if(type(KC) != "t_INT", KC = 20);
ctl(p, modu, res) = { my(b = #binary(p), cs = 0., cm = 0., r); \\ random primes of b bits, == res mod modu
  for(i = 1, KC, until(ispseudoprime(r), r = (random(2^(b-1)) + 2^(b-1)); r = r - (r % modu) + res; if(r < 2^(b-1), r += modu));
    cs += smbits(r + 1); cm += smbits(r - 1)); [cs/KC, cm/KC]; }
{
setrand(12345);
print("# fcyc analysis rows");
for(f = 4, FMAX, if(isprime(f), next); P = phistar(f); if(P > 3 && ispseudoprime(P), analyze(P, "fcy", f, Dl)));
print("# inheritance: fam par bits smooth(p+1) smooth(p-1) ctl_smooth(p+1) ctl_smooth(p-1)");
for(e = 4, EMAX, if(isprime(e), next); P = polcyclo(e, 2); P /= gcd(P, e);
  if(P > 2^20 && ispseudoprime(P), c = ctl(P, e, 1); printf("INH cyc %d %d %.1f %.1f %.1f %.1f\n", e, #binary(P), smbits(P+1), smbits(P-1), c[1], c[2])));
for(f = 4, FMAX, if(isprime(f), next); P = phistar(f);
  if(P > 2^20 && ispseudoprime(P), eps = kronecker(5, P); c = ctl(P, f, (eps + f) % f);
     printf("INH fcy %d %d %.1f %.1f %.1f %.1f\n", f, #binary(P), smbits(P+1), smbits(P-1), c[1], c[2])));
}
