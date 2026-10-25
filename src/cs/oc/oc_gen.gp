\\ Line OC: generate members of each family (PRP via ispseudoprime) up to size caps, and analyse them.
\\ Usage: echo 'FAM="wag"; CAP=1500; read("oc_gen.gp")' | gp -q
read("oc_fam.gp");
Dl = [5,-7,-11,13,-15,17,-19,-23];
{
if(FAM == "mer", forprime(q = 3, CAP, p = 2^q - 1; if(ispseudoprime(p), analyze(p, "mer", q, Dl))));
if(FAM == "wag", forprime(q = 5, CAP, p = (2^q + 1)/3; if(ispseudoprime(p), analyze(p, "wag", q, Dl))));
if(FAM == "fib", for(f = 3, CAP, if(!isprime(f) && f != 4, next); p = fibonacci(f); if(p > 3 && ispseudoprime(p), analyze(p, "fib", f, Dl))));
if(FAM == "luc", for(f = 2, CAP, if(!isprime(f) && f != 2^valuation(f,2), next); p = fibonacci(f-1) + fibonacci(f+1); if(p > 3 && ispseudoprime(p), analyze(p, "luc", f, Dl))));
if(FAM == "gm", forprime(q = 3, CAP, k = (q-1)/2;
     p = 2^q + 2^(k+1) + 1; if(ispseudoprime(p), analyze(p, "gm+", q, Dl));
     p = 2^q - 2^(k+1) + 1; if(p > 5 && ispseudoprime(p), analyze(p, "gm-", q, Dl))));
if(FAM == "cyc", for(e = 3, CAP, if(isprime(e), next); P = polcyclo(e, 2); P /= gcd(P, e);
     if(P > 3 && ispseudoprime(P), analyze(P, "cyc", e, Dl))));
if(FAM == "rie", forstep(c = 1, CMAX, 2, for(m = 2, CAP, p = c*2^m - 1; if(p > 3 && ispseudoprime(p), analyze(p, Str("rie", c), m, Dl)))));
if(FAM == "pro", forstep(c = 1, CMAX, 2, for(m = 2, CAP, p = c*2^m + 1; if(p > 5 && ispseudoprime(p), analyze(p, Str("pro", c), m, Dl)))));
}
