default(parisizemax, 3*10^9);
\\ The three split cofactors skipped by the cost filter in run11 (Theorem B at 1e11). testm computes each unitary-split
\\ gcd modulo the smaller number, so only the actual min(W1, W2*log2|alpha|) matters.
{foreach([[-7,2560959248555062987,1,1],[-15,86495844719236350721,0,0],[-23,69444174200193754009,2,1]], t,
  my(t0=getabstime(), r=testm(t[1], t[2], t[3], t[4]));
  print("D=",t[1]," m=",t[2]," class (",t[3],",",t[4],"): testm ", r, "  time ", (getabstime()-t0)\1000, " s"))}
