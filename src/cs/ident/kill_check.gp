\\ src/cs/ident/kill_check.gp -- (1) numeric check of the two "universal kill" corollaries of Theorem U:
\\  (a) n = (6t+1)(12t+1)(18t+1) (all prime) is never a Lucas psp U_{n+1}(1,Q) == 0 mod n, for ANY Q (D = 1-4Q):
\\      member 12t+1 has odd(N_2) = 1 and v2(12t+2) = 1, so z | 2, impossible (U_1 = U_2 = 1).
\\  (b) n = (6t-1)(12t-1)(18t-1) (all prime) is never a base-2 Fermat psp: member 12t-1 has odd(N''_2) = 1 and
\\      v2(12t-2) = 1, so ord_p 2 | 2, i.e. p | 3.
\\ (2) census: fraction of Korselt-type sets S (F-form, lcm | sum) having a member with odd(N_d) = 1 and 4 | d*lcm(S)
\\     (then that member forces z | 2: the whole family has no Lucas psp with P = 1, for any Q).
read("univ_lib.gp");
lucasU(k,Q,n)={my(D=1-4*Q); lucUV(k,Q,D,n)[1]};
na=0; bada=0; nb=0; badb=0;
{
for(t=1,20000, my(P=[6*t+1,12*t+1,18*t+1]); if(ispseudoprime(P[1])&&ispseudoprime(P[2])&&ispseudoprime(P[3]),
   my(n=P[1]*P[2]*P[3]); na++; for(Q=-60,60, if(Q==0||gcd(Q*(1-4*Q),n)>1||issquare(1-4*Q), next); if(lucasU(n+1,Q,n)==0, bada++; print("(a) counterexample t=",t," Q=",Q)))));
for(t=1,20000, my(P=[6*t-1,12*t-1,18*t-1]); if(ispseudoprime(P[1])&&ispseudoprime(P[2])&&ispseudoprime(P[3]),
   my(n=P[1]*P[2]*P[3]); nb++; if(Mod(2,n)^(n-1)==1, badb++; print("(b) counterexample t=",t))));
}
print("(a) Chernick tuples t<=20000: ",na,", Lucas psp for some Q in [-60,60]: ",bada);
print("(b) Lucas-Carmichael tuples t<=20000: ",nb,", base-2 Fermat psp: ",badb);
census(Dmax,k)={my(nS=0,nk=0); forsubset([Dmax,k],v, my(S=Vec(v),L=lcm(S)); if(gcd(S)!=1||vecsum(S)%L,next); nS++;
   my(kill=0); for(i=1,k, if(oddp(abs(NF(S,i)))==1 && (S[i]*L)%4==0, kill=1)); nk+=kill); [nS,nk]};
for(k=3,6, my(r=census(if(k<=4,120,if(k==5,60,40)),k)); print("k=",k,": Korselt sets=",r[1]," with a universal-kill member=",r[2]));
