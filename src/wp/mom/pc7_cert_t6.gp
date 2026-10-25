read("pc7_cert_shape.gp");
\\ reality interval at (-3/10, 4/5): roots of (w lam - e1)^2 - 4 (e2 + w mu) and w^2 lam^2 - 4 w (mu + lam e1) + 4
{
  my(R1 = -3/10, R2 = 4/5, E1 = R1 + R2, E2 = R1*R2, Dl = (1 - R1^2)*(1 - R2^2), LA = 2*E1/Dl, MU = 2*(1 - E1^2 + E2)/Dl);
  my(q1 = (x*LA - E1)^2 - 4*(E2 + x*MU), q2 = x^2*LA^2 - 4*x*(MU + LA*E1) + 4);
  print("disc_- roots: ", polroots(q1), "\ndisc_+ roots: ", polroots(q2));
}
dump("pc7_cert_t6.dat", [[-3/10, 4/5]]);
