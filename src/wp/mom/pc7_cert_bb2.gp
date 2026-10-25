\\ reality interval and exact B&B at a second shape (-57/80, 19/20) (a grid shape of q22_k7_lp.py where constant alpha fails
\\ without reality: Farkas witness pc7_cert_obs2.gp)
read("pc7_cert_bb.gp");
{
  my(R1 = -57/80, R2 = 19/20, E1 = R1 + R2, E2 = R1*R2, Dl = (1 - R1^2)*(1 - R2^2), LA = 2*E1/Dl, MU = 2*(1 - E1^2 + E2)/Dl);
  print("disc_- roots: ", polroots((x*LA - E1)^2 - 4*(E2 + x*MU))~, "\ndisc_+ roots: ", polroots(x^2*LA^2 - 4*x*(MU + LA*E1) + 4)~);
}
read("pc7_cert_alpha2.gp");
SHAPE = [-57/80, 19/20]; LO = 2/25; HI = 11/25;
foreach([4, 8, 12, 16], n, N = n; print("== N = ", N); if(bb() > 0, break));
