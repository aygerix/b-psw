\\ positive control: n = 5*13*17 style numbers are not BPSW; instead verify that lastq recovers q for a planted case:
\\ pick an inert q and a split p with O_q | p-1 artificially is hard, so test consistency: every q returned must satisfy
\\ odd(ord_q 2)*odd(z_q) | odd(m-1) and class (j1,j2)-compatible orders dividing 2^j1 A, 2^j2 B.

zrank(p,D,Q)={my(mm=p-kronecker(D,p), f=factor(mm)); for(i=1,#f~, my(r=f[i,1]);
  for(e=1,f[i,2], if(lucUV(mm/r,Q,D,p)[1]==0, mm/=r, break))); mm};
chk(D,m,j1,j2)={my(Q=(1-D)/4, R=lastq(D,m,j1,j2), bad=0);
  foreach(R[1], q, my(o=znorder(Mod(2,q)), z=zrank(q,D,Q)); if((oddp(m-1)%(oddp(o)*oddp(z)) && kronecker(D,q)==-1) || valuation(o,2)>j1 || valuation(z,2)>j2, bad++));
  printf("D=%d m=%d: %d candidates, %d inconsistent, A-subsets %d\n", D, m, #R[1], bad, R[2]); R[1]};
\\ m = 61681 (D=5, class (3,1)); m = 127 (D=13, class (0,0)); m = 5581 (D=-11, class (2,0))
print(chk(5,61681,3,1)); print(chk(13,127,0,0)); print(chk(-11,5581,2,0));
