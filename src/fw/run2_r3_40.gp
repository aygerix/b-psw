\\ Forced + one free Williams member (D = 5). Members at Y = L_2K: forced Psi_{d_i}(Y) (Korselt
\\ Psi_d - 1 | n - 1; z | 2Kd forced) and one free f in Q[Y] with f - 1 | n - 1, f + 1 | n + 1.
\\ With P = prod Psi_{d_i}: f - 1 = c D1, f + 1 = c D2 with D1, D2 monic divisors of P - 1,
\\ D2 - D1 = kappa constant, c = 2/kappa, and L = lcm(Psi_{d_i} - 1) | f P - 1. Exact over Q.
\\ Usage: gp -q fw1.gp  (edit R, DMAX below)
default(parisizemax, 2*10^9);
Psi(d)={my(f); if(d==1,return('Y-2)); if(d==2,return('Y+2));
  f=polresultant(polcyclo(d,'w), 'w^2-'Y*'w+1, 'w); f=factor(f)[1,1]; f/pollead(f)};
divs(F)={my(L=List([1])); for(i=1,#F~, my(g=F[i,1]/pollead(F[i,1]), e=F[i,2], M=List());
   foreach(L, D, my(t=D); for(k=0,e, listput(M,t); t*=g)); L=M); Vec(L)};
search(S)={my(P=prod(i=1,#S,PS[S[i]]), L=1, F, DV, n2, res=List());
  for(i=1,#S, L=lcm(L, PS[S[i]]-1));
  F=factor(P-1); n2=poldegree(P-1)\2;
  DV=select(D->poldegree(D)>=1 && poldegree(D)<=n2, divs(F));
  for(i=1,#DV, my(D1=DV[i]); for(j=1,#DV, my(D2=DV[j], k); if(j==i || poldegree(D2)!=poldegree(D1), next);
     k=D2-D1; if(k==0 || poldegree(k)>0, next); k=polcoeff(k,0); if(k<0, next);
     if(((P-1)%(D1*D2))!=0, next);
     my(f=1+(2/k)*D1);
     pre++; if(!mapisdefined(PSM,f), npsi++; if(!mapisdefined(FM,f), mapput(FM,f,S); print("NEWF ",S," f=",f," D1=",D1," D2=",D2)));
     if(((f*P-1)%L)==0, listput(res,[S,f,D1,D2,k]))));
  Vec(res)};
DMAX=40; R=3;
PS=vector(DMAX,d,Psi(d));
DL=select(d->d>=3 && d%5, [1..DMAX]);
cnt=0; hits=0; pre=0; npsi=0; FM=Map(); PSM=Map(); for(e=1,4*DMAX, mapput(PSM,Psi(e),e));
{forsubset([#DL,R], s, my(S=vector(R,i,DL[s[i]]), h=search(S)); cnt++;
  for(t=1,#h, hits++; print("HIT ",h[t])));}
print("R=",R," DMAX=",DMAX," sets tested: ",cnt,"  (D1,D2) pairs before Korselt check: ",pre,"  of which f not a Psi_e: ",npsi,"  distinct non-Psi f: ",#FM,"  hits: ",hits);
