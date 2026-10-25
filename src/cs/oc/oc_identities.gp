\\ Line OC: check the algebraic identities used for the "largest piece" column, and the cyclotomic content of p + 1
\\ for Aurifeuillian (Gaussian-Mersenne) primes.  All exact.
F(n) = fibonacci(n);
L(n) = fibonacci(n-1) + fibonacci(n+1);
ok = 1;
nchk = 0;
{for(k = 1, 60, nchk++; \
  if(F(4*k+1) - 1 != F(2*k)*L(2*k+1), ok = 0; print("F4k+1 - 1 fails ", k)); \
  if(F(4*k+3) - 1 != F(2*k+2)*L(2*k+1), ok = 0; print("F4k+3 - 1 fails ", k)); \
  if(F(4*k+1) + 1 != F(2*k+1)*L(2*k), ok = 0; print("F4k+1 + 1 fails ", k)); \
  if(F(4*k+3) + 1 != F(2*k+1)*L(2*k+2), ok = 0; print("F4k+3 + 1 fails ", k)); \
  if(L(4*k+1) - 1 != 5*F(2*k)*F(2*k+1), ok = 0; print("L4k+1 - 1 fails ", k)); \
  if(L(4*k+3) - 1 != L(2*k+1)*L(2*k+2), ok = 0; print("L4k+3 - 1 fails ", k)); \
  if(L(4*k+1) + 1 != L(2*k)*L(2*k+1), ok = 0; print("L4k+1 + 1 fails ", k)); \
  if(L(4*k+3) + 1 != 5*F(2*k+1)*F(2*k+2), ok = 0; print("L4k+3 + 1 fails ", k)));}
print("Fibonacci/Lucas -+1 identities, ", nchk, " values of k: ", if(ok && nchk == 60, "ALL HOLD", "FAIL"));
\\ Aurifeuillian: for odd q, k = (q-1)/2, p+- = 2^q +- 2^(k+1) + 1, (p+- + 1)/2 = x^(2k) +- x^k + 1 at x = 2.
\\ largest cyclotomic index in x^(2k) + x^k + 1 = Phi_3(x^k) and x^(2k) - x^k + 1 = Phi_6(x^k)
cycidx(P) = { my(f = factor(P)[,1], best = 0, d); for(i = 1, #f, d = 1; while(polcyclo(d) != f[i], d++; if(d > 20000, return(-1))); best = max(best, d)); best; }
{ for(q = 3, 61, if(!isprime(q), next); k = (q-1)/2;
    if(2^q + 2^(k+1) + 1 + 1 != 2*(2^(2*k) + 2^k + 1) || 2^q - 2^(k+1) + 1 + 1 != 2*(2^(2*k) - 2^k + 1), print("GM identity fails ", q));
    d3 = cycidx(x^(2*k) + x^k + 1); d6 = cycidx(x^(2*k) - x^k + 1);
    if(q < 40, print("q=", q, " k=", k, "  p+ : largest cyclotomic index ", d3, " (phi ", eulerphi(d3), ")   p- : ", d6, " (phi ", eulerphi(d6), ")"))); }
\\ Wagstaff: (p+1)/4 = (2^(q-2)+1)/3, largest index 2(q-2), phi(q-2)
{ for(q = 5, 40, if(!isprime(q), next); P = (2^q+1)/3; if((P+1) != 4*(2^(q-2)+1)/3, print("Wag fails ", q))); print("Wagstaff p+1 = 4(2^(q-2)+1)/3 checked for q < 40"); }
