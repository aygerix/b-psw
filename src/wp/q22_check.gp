\\ Exact check of the (2,2) quartic pencil reduction (gen_quartic22.py): with rational u1, u2, c_i,
\\ F == prod_i (1 + (1 - c_i/c_j) X) mod K-_j  and  F == prod_i (-1 + (1 - c_i/c_j) Y) mod K+_j.
ok = 1;
{
for(trial = 1, 30,
  u1 = random(1000)/1001 - 1/2; u2 = random(1000)/1001 + 1/3; if(u1 == u2, next);
  G = t^2 - 1; H = (t - u1)*(t - u2);
  \\ f0 = 1 + G (al t + be), f0(u) = -1 at u1, u2
  M = [(u1^2-1)*u1, (u1^2-1); (u2^2-1)*u2, (u2^2-1)]; s = matsolve(M, [-2, -2]~); al = s[1]; be = s[2];
  f0 = 1 + G*(al*t + be); P = G*H; X = f0 - 1; Y = f0 + 1;
  r0 = Y / H; if(type(r0) != "t_POL", ok = 0; print("H does not divide f0+1"));
  k = 5; c = vector(k, i, (random(97) + 1)/(random(13) + 1)); fs = vector(k, i, f0 + c[i]*P); F = prod(i = 1, k, fs[i]);
  j = random(k) + 1;
  Km = al*t + be + c[j]*H; Kp = r0 + c[j]*G;
  pm = prod(i = 1, k, 1 + (1 - c[i]/c[j])*X); pp = prod(i = 1, k, -1 + (1 - c[i]/c[j])*Y);
  if((F - pm) % Km != 0 || (F - pp) % Kp != 0, ok = 0; print("FAIL at trial ", trial)));
}
print(if(ok, "exact identity verified on 30 random rational pencils", "FAILED"));
