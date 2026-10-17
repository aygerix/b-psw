// k4q.sing translated to Calyx (Magma language) by sing2calyx.py. Run: calyx < k4q.m
// Singular: ring R=0,(c1,a2,a3,a4,b2,b3,b4,c2,c3,c4,u2,u3,u4),dp; option(redSB); std(I); dim; size.
SetColumns(0);
P<c1,a2,a3,a4,b2,b3,b4,c2,c3,c4,u2,u3,u4> := PolynomialRing(Rationals(), 13, "grevlex");
I := ideal<P |
2*a2*a3*a4*c1^3 - 2*a2*a3*a4*c1 - a2*a3*c1^2*c4 + a2*a3*c4 - a2*a4*c1^2*c3 + a2*a4*c3 - a2*b3*b4*c1^2 + a2*b3*b4 - a3*a4*c1^2*c2 + a3*a4*c2 - a3*b2*b4*c1^2 + a3*b2*b4 - a4*b2*b3*c1^2 + a4*b2*b3 + c2*c3*c4 - 1,
-a2*a3*b4*c1^2 + a2*a3*b4 - a2*a4*b3*c1^2 + a2*a4*b3 - a3*a4*b2*c1^2 + a3*a4*b2 + b2*c3*c4 + b3*c2*c4 + b4*c2*c3,
3*a2*a3*a4*c1^2 + a2*a3*a4 - 2*a2*a3*c1*c4 - 2*a2*a4*c1*c3 - 2*a2*b3*b4*c1 + a2*c3*c4 - 2*a3*a4*c1*c2 - 2*a3*b2*b4*c1 + a3*c2*c4 - 2*a4*b2*b3*c1 + a4*c2*c3 + b2*b3*c4 + b2*b4*c3 + b3*b4*c2,
-2*a2*a3*b4*c1 - 2*a2*a4*b3*c1 + a2*b3*c4 + a2*b4*c3 - 2*a3*a4*b2*c1 + a3*b2*c4 + a3*b4*c2 + a4*b2*c3 + a4*b3*c2 + b2*b3*b4,
-4*a2^2*a3*a4*b2^2*c2^2*u2^6 + 4*a2^2*a3*a4*b2^2*u2^6 + 2*a2*a3*a4*c2^3*u2^4 - 2*a2*a3*a4*c2*u2^4 + 2*a2*a3*b2*b4*c2^2*u2^4 - 2*a2*a3*b2*b4*u2^4 + 2*a2*a4*b2*b3*c2^2*u2^4 - 2*a2*a4*b2*b3*u2^4 + a3*a4*b2^2*c2^2*u2^4 - a3*a4*b2^2*u2^4 - a3*a4*c1*c2^2*u2^2 + a3*a4*c1*u2^2 - a3*c2^2*c4*u2^2 + a3*c4*u2^2 - a4*c2^2*c3*u2^2 + a4*c3*u2^2 - b3*b4*c2^2*u2^2 + b3*b4*u2^2 + c1*c3*c4 - 1,
-8*a2^2*a3*a4*b2^3*c2*u2^6 + 6*a2*a3*a4*b2*c2^2*u2^4 - 2*a2*a3*a4*b2*u2^4 + 4*a2*a3*b2^2*b4*c2*u2^4 + 4*a2*a4*b2^2*b3*c2*u2^4 + 2*a3*a4*b2^3*c2*u2^4 - 2*a3*a4*b2*c1*c2*u2^2 - 2*a3*b2*c2*c4*u2^2 - a3*b4*c2^2*u2^2 + a3*b4*u2^2 - 2*a4*b2*c2*c3*u2^2 - a4*b3*c2^2*u2^2 + a4*b3*u2^2 - 2*b2*b3*b4*c2*u2^2 + b3*c1*c4 + b4*c1*c3,
-8*a2^3*a3*a4*b2^2*c2*u2^6 - 4*a2^2*a3*a4*b2^4*u2^6 + 4*a2^2*a3*a4*c2^2*u2^4 + 4*a2^2*a3*b2*b4*c2*u2^4 + 4*a2^2*a4*b2*b3*c2*u2^4 + 8*a2*a3*a4*b2^2*c2*u2^4 - 2*a2*a3*a4*c1*c2*u2^2 + 2*a2*a3*b2^3*b4*u2^4 - 2*a2*a3*c2*c4*u2^2 + 2*a2*a4*b2^3*b3*u2^4 - 2*a2*a4*c2*c3*u2^2 - 2*a2*b3*b4*c2*u2^2 + a3*a4*b2^4*u2^4 - a3*a4*b2^2*c1*u2^2 - a3*a4*c2^2*u2^2 + a3*a4*u2^2 - a3*b2^2*c4*u2^2 - 2*a3*b2*b4*c2*u2^2 + a3*c1*c4 - a4*b2^2*c3*u2^2 - 2*a4*b2*b3*c2*u2^2 + a4*c1*c3 - b2^2*b3*b4*u2^2 + b3*b4*c1 + c3*c4,
-8*a2^3*a3*a4*b2^3*u2^6 + 8*a2^2*a3*a4*b2*c2*u2^4 + 4*a2^2*a3*b2^2*b4*u2^4 + 4*a2^2*a4*b2^2*b3*u2^4 + 4*a2*a3*a4*b2^3*u2^4 - 2*a2*a3*a4*b2*c1*u2^2 - 2*a2*a3*b2*c4*u2^2 - 2*a2*a3*b4*c2*u2^2 - 2*a2*a4*b2*c3*u2^2 - 2*a2*a4*b3*c2*u2^2 - 2*a2*b2*b3*b4*u2^2 - 2*a3*a4*b2*c2*u2^2 - a3*b2^2*b4*u2^2 + a3*b4*c1 - a4*b2^2*b3*u2^2 + a4*b3*c1 + b3*c4 + b4*c3,
-4*a2*a3^2*a4*b3^2*c3^2*u3^6 + 4*a2*a3^2*a4*b3^2*u3^6 + 2*a2*a3*a4*c3^3*u3^4 - 2*a2*a3*a4*c3*u3^4 + 2*a2*a3*b3*b4*c3^2*u3^4 - 2*a2*a3*b3*b4*u3^4 + a2*a4*b3^2*c3^2*u3^4 - a2*a4*b3^2*u3^4 - a2*a4*c1*c3^2*u3^2 + a2*a4*c1*u3^2 - a2*c3^2*c4*u3^2 + a2*c4*u3^2 + 2*a3*a4*b2*b3*c3^2*u3^4 - 2*a3*a4*b2*b3*u3^4 - a4*c2*c3^2*u3^2 + a4*c2*u3^2 - b2*b4*c3^2*u3^2 + b2*b4*u3^2 + c1*c2*c4 - 1,
-8*a2*a3^2*a4*b3^3*c3*u3^6 + 6*a2*a3*a4*b3*c3^2*u3^4 - 2*a2*a3*a4*b3*u3^4 + 4*a2*a3*b3^2*b4*c3*u3^4 + 2*a2*a4*b3^3*c3*u3^4 - 2*a2*a4*b3*c1*c3*u3^2 - 2*a2*b3*c3*c4*u3^2 - a2*b4*c3^2*u3^2 + a2*b4*u3^2 + 4*a3*a4*b2*b3^2*c3*u3^4 - a4*b2*c3^2*u3^2 + a4*b2*u3^2 - 2*a4*b3*c2*c3*u3^2 - 2*b2*b3*b4*c3*u3^2 + b2*c1*c4 + b4*c1*c2,
-8*a2*a3^3*a4*b3^2*c3*u3^6 - 4*a2*a3^2*a4*b3^4*u3^6 + 4*a2*a3^2*a4*c3^2*u3^4 + 4*a2*a3^2*b3*b4*c3*u3^4 + 8*a2*a3*a4*b3^2*c3*u3^4 - 2*a2*a3*a4*c1*c3*u3^2 + 2*a2*a3*b3^3*b4*u3^4 - 2*a2*a3*c3*c4*u3^2 + a2*a4*b3^4*u3^4 - a2*a4*b3^2*c1*u3^2 - a2*a4*c3^2*u3^2 + a2*a4*u3^2 - a2*b3^2*c4*u3^2 - 2*a2*b3*b4*c3*u3^2 + a2*c1*c4 + 4*a3^2*a4*b2*b3*c3*u3^4 + 2*a3*a4*b2*b3^3*u3^4 - 2*a3*a4*c2*c3*u3^2 - 2*a3*b2*b4*c3*u3^2 - 2*a4*b2*b3*c3*u3^2 - a4*b3^2*c2*u3^2 + a4*c1*c2 - b2*b3^2*b4*u3^2 + b2*b4*c1 + c2*c4,
-8*a2*a3^3*a4*b3^3*u3^6 + 8*a2*a3^2*a4*b3*c3*u3^4 + 4*a2*a3^2*b3^2*b4*u3^4 + 4*a2*a3*a4*b3^3*u3^4 - 2*a2*a3*a4*b3*c1*u3^2 - 2*a2*a3*b3*c4*u3^2 - 2*a2*a3*b4*c3*u3^2 - 2*a2*a4*b3*c3*u3^2 - a2*b3^2*b4*u3^2 + a2*b4*c1 + 4*a3^2*a4*b2*b3^2*u3^4 - 2*a3*a4*b2*c3*u3^2 - 2*a3*a4*b3*c2*u3^2 - 2*a3*b2*b3*b4*u3^2 - a4*b2*b3^2*u3^2 + a4*b2*c1 + b2*c4 + b4*c2,
-4*a2*a3*a4^2*b4^2*c4^2*u4^6 + 4*a2*a3*a4^2*b4^2*u4^6 + 2*a2*a3*a4*c4^3*u4^4 - 2*a2*a3*a4*c4*u4^4 + a2*a3*b4^2*c4^2*u4^4 - a2*a3*b4^2*u4^4 - a2*a3*c1*c4^2*u4^2 + a2*a3*c1*u4^2 + 2*a2*a4*b3*b4*c4^2*u4^4 - 2*a2*a4*b3*b4*u4^4 - a2*c3*c4^2*u4^2 + a2*c3*u4^2 + 2*a3*a4*b2*b4*c4^2*u4^4 - 2*a3*a4*b2*b4*u4^4 - a3*c2*c4^2*u4^2 + a3*c2*u4^2 - b2*b3*c4^2*u4^2 + b2*b3*u4^2 + c1*c2*c3 - 1,
-8*a2*a3*a4^2*b4^3*c4*u4^6 + 6*a2*a3*a4*b4*c4^2*u4^4 - 2*a2*a3*a4*b4*u4^4 + 2*a2*a3*b4^3*c4*u4^4 - 2*a2*a3*b4*c1*c4*u4^2 + 4*a2*a4*b3*b4^2*c4*u4^4 - a2*b3*c4^2*u4^2 + a2*b3*u4^2 - 2*a2*b4*c3*c4*u4^2 + 4*a3*a4*b2*b4^2*c4*u4^4 - a3*b2*c4^2*u4^2 + a3*b2*u4^2 - 2*a3*b4*c2*c4*u4^2 - 2*b2*b3*b4*c4*u4^2 + b2*c1*c3 + b3*c1*c2,
-8*a2*a3*a4^3*b4^2*c4*u4^6 - 4*a2*a3*a4^2*b4^4*u4^6 + 4*a2*a3*a4^2*c4^2*u4^4 + 8*a2*a3*a4*b4^2*c4*u4^4 - 2*a2*a3*a4*c1*c4*u4^2 + a2*a3*b4^4*u4^4 - a2*a3*b4^2*c1*u4^2 - a2*a3*c4^2*u4^2 + a2*a3*u4^2 + 4*a2*a4^2*b3*b4*c4*u4^4 + 2*a2*a4*b3*b4^3*u4^4 - 2*a2*a4*c3*c4*u4^2 - 2*a2*b3*b4*c4*u4^2 - a2*b4^2*c3*u4^2 + a2*c1*c3 + 4*a3*a4^2*b2*b4*c4*u4^4 + 2*a3*a4*b2*b4^3*u4^4 - 2*a3*a4*c2*c4*u4^2 - 2*a3*b2*b4*c4*u4^2 - a3*b4^2*c2*u4^2 + a3*c1*c2 - 2*a4*b2*b3*c4*u4^2 - b2*b3*b4^2*u4^2 + b2*b3*c1 + c2*c3,
-8*a2*a3*a4^3*b4^3*u4^6 + 8*a2*a3*a4^2*b4*c4*u4^4 + 4*a2*a3*a4*b4^3*u4^4 - 2*a2*a3*a4*b4*c1*u4^2 - 2*a2*a3*b4*c4*u4^2 + 4*a2*a4^2*b3*b4^2*u4^4 - 2*a2*a4*b3*c4*u4^2 - 2*a2*a4*b4*c3*u4^2 - a2*b3*b4^2*u4^2 + a2*b3*c1 + 4*a3*a4^2*b2*b4^2*u4^4 - 2*a3*a4*b2*c4*u4^2 - 2*a3*a4*b4*c2*u4^2 - a3*b2*b4^2*u4^2 + a3*b2*c1 - 2*a4*b2*b3*b4*u4^2 + b2*c3 + b3*c2,
a2*u2 - 1,
a3*u3 - 1,
a4*u4 - 1
>;
printf "generators: %o, variables: %o\n", #Basis(I), Rank(P);
t0 := Cputime(); r0 := Realtime();
G := GroebnerBasis(I);
printf "GroebnerBasis: %o s cpu, %o s real\n", Cputime(t0), Realtime(r0);
t1 := Cputime();
d := Dimension(I);
printf "Dimension: %o s cpu\n", Cputime(t1);
print "dim:"; print d; print "size:"; print #G;
F := "/Users/jmercedes/Desktop/Projects/Research/b-psw/src/alg/k4q_calyx.out";
PrintFile(F, "// reduced Groebner basis of k4q (grevlex, " cat Sprint(#G) cat " elements), Krull dimension " cat Sprint(d) : Overwrite := true);
for g in G do PrintFile(F, g); end for;
quit;
