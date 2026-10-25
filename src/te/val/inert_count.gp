\\ independent count of inert primes p <= X (kronecker(D,p) = -1, p not dividing 2QD), Q = (1-D)/4
ic(D, X) = { my(Q = (1 - D)/4, c = 0); forprime(p = 3, X, if((2*Q*D) % p && kronecker(D, p) == -1, c++)); c; }
{
  foreach([[5, 10^6], [5, 10^7], [5, 10^8], [-7, 10^6], [-7, 10^7], [13, 10^6]], v,
    print("INERT D=", v[1], " X=", v[2], " count=", ic(v[1], v[2])));
}
quit
