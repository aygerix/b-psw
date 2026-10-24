\\ Lemma 15 checks. (1) the identity prod_cyc(p_j p_k - 1) - prod(p_i^2 - 1) = 1/2 sum_cyc (p_i^2 - 1)(p_j - p_k)^2.
\\ (2) the bound on N for x1 < x2 < x3 on random triples.
\\ (3) live test on 3-factor odd squarefree n <= 2e7 with all p inert for D = 5: the c_i = (n/p_i - 1)/O_i are
\\ not integers in general, but whenever all three are integers the inequality prod I_i >= prod(p_i^2-1)/N must hold.
read("../fam/bpsw.gp");
zrank(p,D,Q)={my(m=p-kronecker(D,p), f=factor(m)); for(i=1,#f~, my(r=f[i,1]);
  for(e=1,f[i,2], if(lucUV(m/r,Q,D,p)[1]==0, m/=r, break))); m};
oddp(x)=x>>valuation(x,2);
L=(a*b-1)*(a*c-1)*(b*c-1)-(a^2-1)*(b^2-1)*(c^2-1);
R=((a^2-1)*(b-c)^2+(b^2-1)*(a-c)^2+(c^2-1)*(a-b)^2)/2;
print("identity holds: ", L==R);
{my(bad=0); for(t=1,10^5, my(x=vecsort([random(10^9)+2,random(10^9)+2,random(10^9)+2])); if(#Set(x)<3, next);
  my(N=subst(subst(subst(R,a,x[1]),b,x[2]),c,x[3])); if(N>=3/2*x[2]^2*x[3]^2, bad++));
 print("random triples violating N < 3/2 p2^2 p3^2: ", bad);}
\\ (3): all products of three inert primes (D=5) below 2e7 for which every c_i is an integer
{my(P=select(p->kronecker(5,p)==-1, primes([3,30000])), O=vector(#P), cnt=0, ok=0);
 for(i=1,#P, my(p=P[i]); O[i]=oddp(znorder(Mod(2,p)))*oddp(zrank(p,5,-1)));
 for(i=1,#P, for(j=i+1,#P, if(P[i]*P[j]^2>2*10^7, break); for(k=j+1,#P, my(n=P[i]*P[j]*P[k]); if(n>2*10^7, break);
   if((P[j]*P[k]-1)%O[i] || (P[i]*P[k]-1)%O[j] || (P[i]*P[j]-1)%O[k], next); cnt++;
   my(I=prod(t=1,3,(P[[i,j,k][t]]^2-1)/O[[i,j,k][t]]), x=[P[i],P[j],P[k]],
      N=((x[1]^2-1)*(x[2]-x[3])^2+(x[2]^2-1)*(x[1]-x[3])^2+(x[3]^2-1)*(x[1]-x[2])^2)/2);
   if(I*N>=prod(t=1,3,x[t]^2-1), ok++); print("  all c_i integral: ",x," prod I = ",I," bound ",prod(t=1,3,x[t]^2-1)/N*1.))));
 print("triples with all c_i integral: ",cnt,", satisfying the Lemma 15 bound: ",ok);}
\\ (4) non-vacuous test of the integrality argument: for ANY distinct odd x_i > 1 put O_i = gcd(x_i^2-1, x_j x_k - 1)
\\ (largest admissible order), I_i = (x_i^2-1)/O_i; the lemma's inequality must hold for every triple.
{my(bad=0, tight=0, tot=0); forstep(a=3,301,2, forstep(b=a+2,301,2, forstep(c=b+2,301,2, tot++;
  my(x=[a,b,c], I=1); for(t=1,3, my(u=x[t], v=prod(s=1,3,if(s==t,1,x[s]))); I*=(u^2-1)/gcd(u^2-1, v-1));
  my(N=((a^2-1)*(b-c)^2+(b^2-1)*(a-c)^2+(c^2-1)*(a-b)^2)/2, B=prod(t=1,3,x[t]^2-1));
  if(I*N<B, bad++); if(I*N<2*B, tight++))));
 print("odd triples <= 301: ",tot,"  violations of prod I >= prod(x^2-1)/N: ",bad,"  within factor 2 of the bound: ",tight);}
