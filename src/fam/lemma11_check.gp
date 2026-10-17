\\ Checks for Lemma 11 (Lemma 5 family, Y >= 3): gcd(Y(Y-1), n+1) | 2, gcd(Y(Y^2-3), n+1) | 8, (Y^2-2) | n+1;
\\ odd(z(Psi10)), odd(z(Psi18)) for the selected D; inert/split pattern. Y < 0: |Psi14|+1, |Psi18|+1 | n-1.
read("bpsw.gp");
zrank(p,D,Q)={my(m=p-kronecker(D,p), f=factor(m)); for(i=1,#f~, my(r=f[i,1]);
  for(e=1,f[i,2], if(lucUV(m/r,Q,D,p)[1]==0, m/=r, break))); m};
oddp(x)=x>>valuation(x,2);
{my(rows=0, g1=0, g3=0, g2=0, pw=0, allin=0, Dh=Map(), negok=0, negrows=0);
for(Y=-200000,200000, my(a=Y^2-Y-1, b=Y^3-Y^2-2*Y+1, d=Y^3-3*Y-1, n=abs(a*b*d));
  if(vecmin([abs(a),abs(b),abs(d)])<3 || !isprime(abs(a)) || !isprime(abs(b)) || !isprime(abs(d)), next);
  if(Y<0, negrows++; if((n-1)%(abs(b)+1)==0 && (n-1)%(abs(d)+1)==0, negok++); next);
  my(D=selD(n), Q=(1-D)/4); rows++;
  if(2%gcd(Y*(Y-1),n+1)==0, g1++); if(8%gcd(Y*(Y^2-3),n+1)==0, g3++); if((n+1)%(Y^2-2)==0, g2++);
  if(oddp(zrank(a,D,Q))==1 || oddp(zrank(d,D,Q))==1, pw++);
  if(kronecker(D,a)==-1 && kronecker(D,b)==-1 && kronecker(D,d)==-1, allin++);
  my(c=0); if(mapisdefined(Dh,D,&c), mapput(Dh,D,c+1), mapput(Dh,D,1)));
print("Y>0 rows ",rows,": gcd(Y(Y-1),n+1)|2: ",g1,"  gcd(Y(Y^2-3),n+1)|8: ",g3,"  (Y^2-2)|n+1: ",g2);
print("Y>0 rows with odd(z(Psi10))=1 or odd(z(Psi18))=1: ",pw,"   all three inert for selected D: ",allin);
print("Y>0 Selfridge D histogram: ",Mat(Dh));
print("Y<0 rows ",negrows,": |Psi14|+1 and |Psi18|+1 both divide n-1: ",negok);}
