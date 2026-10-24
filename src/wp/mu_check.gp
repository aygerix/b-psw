\\ Lemma 20 check: for odd k, the only common root of h(z) = (1+z)^k - 1 - z^k and h(-z) is z = 0.
\\ (S(X) = (X+mu)^k + nu with S(0)=0, S(1)=1, S(-1)=-1 forces h(mu) = h(-mu) = 0.)
{bad=0;
forstep(k=3,301,2, my(h=(1+z)^k-1-z^k, g=gcd(h, subst(h,z,-z)));
  if(poldegree(g)!=1 || g!=z*pollead(g), bad++; print("k=",k," gcd=",g)));
print("odd k in [3,301]: exceptions ", bad);}
quit
