// k3_elim.sing translated to Calyx (Magma language). Run: calyx < k3_elim.m
// k=3 polynomial families: (x,y,z)=(f1,f2,f3)(t) must lie on
// Q1: yz-1 = c1(x^2-1), Q2: xz-1 = c2(y^2-1), Q3: xy-1 = c3(z^2-1).
// Singular's ring R = 0,(y,z,c1,c2,c3,x),dp is grevlex with y > z > c1 > c2 > c3 > x.
SetColumns(0);
P<y,z,c1,c2,c3,x> := PolynomialRing(Rationals(), 6, "grevlex");
I := ideal<P | y*z-1-c1*(x^2-1), x*z-1-c2*(y^2-1), x*y-1-c3*(z^2-1)>;
t0 := Cputime();
GI := GroebnerBasis(I);
print "== GB(I) grevlex, size:"; print #GI;
for g in GI do print g; end for;
// eliminate(I, y*z): keep the variables c1, c2, c3, x.
K := EliminationIdeal(I, {c1, c2, c3, x});
GK := GroebnerBasis(K);
print "== GB(K) grevlex, size:"; print #GK;
for g in GK do print g; end for;
// Coefficients of each generator of K with respect to x, zeros dropped.
C := [P | ];
for f in Basis(K) do
    C cat:= [g : g in Coefficients(f, x) | g ne 0];
end for;
S<c1,c2,c3> := PolynomialRing(Rationals(), 3, "grevlex");
h := hom<P -> S | [0, 0, S.1, S.2, S.3, 0]>;
J := ideal<S | [h(g) : g in C]>;
G := GroebnerBasis(J);
print "== GB(C) in S grevlex, size:"; print #G;
for g in G do print g; end for;
d := Dimension(J);
print "== dim:"; print d;
// Calyx has no minAssGTZ / RadicalDecomposition. Check Singular's single minimal
// prime Q instead: J is in Q, and Q is in the radical of J, so V(J) = V(Q).
Q := ideal<S | c1+c2+c3-3, c2^2*c3+c2*c3^2-3*c2*c3+1>;
GQ := GroebnerBasis(Q);
print "== GB(Q) grevlex, size:"; print #GQ;
for g in GQ do print g; end for;
print "== J in Q (normal forms of GB(J) mod Q):"; print [NormalForm(g, Q) : g in G];
print "== Q in rad(J):"; print [IsInRadical(q, J) : q in Basis(Q)];
print "== Q proper, dim:"; print IsProper(Q), Dimension(Q);
printf "== cputime %o s\n", Cputime(t0);
quit;
