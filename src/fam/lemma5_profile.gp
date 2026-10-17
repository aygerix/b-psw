\\ G7 property profile of the Lemma 5 family n = Psi10(Y) Psi14(Y) Psi18(Y), all |values| prime.
\\ Per n: Selfridge D, Korselt, strong-Fermat / strong-Lucas (direct, bpsw.gp), the Lemma 1 prediction
\\ from per-prime orders (must agree with the direct tests), A_F, B_F, Fermat/Lucas defects, indices.
\\ Usage (from src/fam): gp -q lemma5_profile.gp < /dev/null > lemma5_profile.out
default(parisizemax, 2*10^9);
read("bpsw.gp");
\\ rank of apparition of U (P=1, Q) mod p, p odd, p not | QD
zrank(p,D,Q)={my(m=p-kronecker(D,p), f=factor(m)); for(i=1,#f~, my(r=f[i,1]);
  for(e=1,f[i,2], if(lucUV(m/r,Q,D,p)[1]==0, m/=r, break))); m};
oddp(x)=x>>valuation(x,2);
lg(x)=if(x<=1,0.,log(x)/log(2));
profile(Y)={my(v=[Y^2-Y-1, Y^3-Y^2-2*Y+1, Y^3-3*Y-1], P=apply(abs,v), n=abs(prod(i=1,3,v[i])));
  if(#Set(P)<3 || vecmin(P)<3, return(0)); for(i=1,3, if(!isprime(P[i]), return(0)));
  my(D=selD(n), Q=(1-D)/4, s=valuation(n-1,2), t=valuation(n+1,2), A=1, B=1, jf=Set(), jl=Set(), kor=1, ix=[], ixf=[], kr=[], bad=0);
  for(i=1,3, my(p=P[i]); if(D%p==0 || Q%p==0, bad=1; next);
    my(o=znorder(Mod(2,p)), z=zrank(p,D,Q));
    A=lcm(A,oddp(o)); B=lcm(B,oddp(z)); jf=setunion(jf,[valuation(o,2)]); jl=setunion(jl,[valuation(z,2)]);
    if((n-1)%(p-1), kor=0); kr=concat(kr,kronecker(D,p)); ix=concat(ix,(p-kronecker(D,p))/z); ixf=concat(ixf,(p-1)/o));
  if(bad, return([Y,n,D,"p|QD"]));
  my(dF=A/gcd(A,n-1), dL=B/gcd(B,n+1), fc=(#jf==1 && jf[1]<=s), lc=(#jl==1 && jl[1]<=t));
  my(f2=sprp2(n), sl=slprp(n));
  [Y, n, D, kr, kor, f2, sl, fc, lc, lg(A), lg(B), lg(dF), lg(dL), lg(n), ix, ixf, (f2==(dF==1&&fc)) && (sl==(dL==1&&lc)), jf, jl]};

YMAX=200000; R=List();
{for(Y=-YMAX,YMAX, my(r=profile(Y)); if(r!=0, if(#r==4, print("skip (p | QD): ",r), listput(R,r))));}
R=Vec(R);
print("rows (|Y| <= ",YMAX,", all three |values| prime, distinct, >= 3): ",#R);
{my(Dh=Map()); for(i=1,#R, my(D=R[i][3], c=0); if(mapisdefined(Dh,D,&c), mapput(Dh,D,c+1), mapput(Dh,D,1)));
 print("Selfridge D histogram [D, count]: ", Mat(Dh));}
{my(pos=select(r->r[1]>0,R), neg=select(r->r[1]<0,R));
 foreach([[pos,"Y>0"],[neg,"Y<0"]], g, my(S=g[1]); if(#S==0, next);
  print(g[2],": rows ",#S,
   "  Korselt ",#select(r->r[5],S),
   "  strong-Fermat(2) ",#select(r->r[6],S),
   "  one Fermat 2-adic class ",#select(r->r[8],S),
   "  Fermat defect 1 ",#select(r->r[12]==0,S),
   "  one Lucas 2-adic class ",#select(r->r[9],S),
   "  Lucas defect 1 ",#select(r->r[13]==0,S),
   "  strong-Lucas ",#select(r->r[7],S),
   "  BPSW ",#select(r->r[6]&&r[7],S),
   "  Lemma-1 prediction mismatches ",#select(r->!r[17],S)));}
\\ size of the Lucas defect relative to n, and A_F relative to n
{my(rl=vecsort(apply(r->r[13]/r[14],R)), ra=vecsort(apply(r->r[10]/r[14],R)), m=#R);
 printf("log2(delta_L)/log2(n): min %.3f  median %.3f  max %.3f\n", rl[1], rl[(m+1)\2], rl[m]);
 printf("log2(A_F)/log2(n):     min %.3f  median %.3f  max %.3f\n", ra[1], ra[(m+1)\2], ra[m]);}
{my(S=vecsort(R,(a,b)->sign(a[13]-b[13]))); print("10 smallest Lucas defects [Y, D, log2 delta_L, log2 n]:");
 for(i=1,min(10,#S), printf("  Y=%d D=%d  %.1f  %.1f\n", S[i][1], S[i][3], S[i][13], S[i][14]));}
\\ index distribution: Lucas index (p-(D/p))/z and Fermat index (p-1)/ord_p 2 over all primes
{my(I=concat(apply(r->r[15],R)), J=concat(apply(r->r[16],R)));
 print("Lucas index (p-(D/p))/z(p): <=2: ",#select(x->x<=2,I),"  <=10: ",#select(x->x<=10,I),"  <=100: ",#select(x->x<=100,I),"  of ",#I);
 print("Fermat index (p-1)/ord_p 2:  <=2: ",#select(x->x<=2,J),"  <=10: ",#select(x->x<=10,J),"  <=100: ",#select(x->x<=100,J),"  of ",#J);}
print("--- rows with |Y| < 20000: [Y, D, (D/p_j), Korselt, sF2, sL, Fclass, Lclass, log2 A, log2 B, log2 dF, log2 dL, log2 n, Lucas idx, Fermat idx]");
{for(i=1,#R, my(r=R[i]); if(abs(r[1])<20000,
  printf("%d %d %s %d %d %d %d %d %.1f %.1f %.1f %.1f %.1f %s %s\n", r[1], r[3], r[4], r[5], r[6], r[7], r[8], r[9], r[10], r[11], r[12], r[13], r[14], r[15], r[16])));}
