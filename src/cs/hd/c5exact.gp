\\ EXACT certificate of the k = 5 two-dead family (all M): members f = 1 + a x with a in {1, y1, y2, y3, M}, where y1, y2, y3
\\ are the roots of P(y, M) = c y^3 + c2 y^2 + c1 y + c0 (below).  Claim: every root y of P is a Lucas-good member, i.e.
\\ prod over the other four slopes s of (y - 2 s) = y^4.  With P(z) = c (z - y)(z - y')(z - y''), (y - 2y')(y - 2y'') = -8 P(y/2)/(c y),
\\ so the claim is the polynomial congruence  8 (y - 2)(y - 2M) P(y/2) + c y^5 == 0  mod P(y)  in Q(M)[y].
P(y, M) = (7*M^3 + 57*M^2 + 57*M + 7)*y^3 + (-98*M^3 - 308*M^2 - 98*M)*y^2 + (392*M^3 + 392*M^2)*y - 448*M^3;
{
  my(c = 7*M^3 + 57*M^2 + 57*M + 7, Py = P('y, 'M), R = 8*('y - 2)*('y - 2*'M)*subst(Py, 'y, 'y/2) + c*'y^5);
  print("remainder of 8(y-2)(y-2M)P(y/2) + c y^5 modulo P(y), in Q(M)[y]:  ", R % Py);
  print("disc_y P = ", factor(poldisc(Py, 'y)));
  \\ dead members: a = 1 and a = M; check they are dead generically: residual prod_{s != 1}(1 - 2s) - 1 and the same for M
  my(r = polroots(P('y, 135)));
  print("M = 135: roots of P = ", r);
  my(v = concat([1], concat(vector(3, i, real(r[i])), [135])));
  print("  prod_(i!=j)(1 - 2a_i/a_j), j = 1..5: ", vector(5, j, prod(i = 1, 5, if (i == j, 1, 1 - 2*v[i]/v[j]))));
  print("  Sturm: number of real roots of P(y,135) in (0, 135): ", polsturm(P('y, 135), [0, 135]));
  print("  exact rational mixed-sign point M = -119/54: P factors as ", factor(P('y, -119/54)));
}
