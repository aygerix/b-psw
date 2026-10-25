/* Recompute the cs_rank per-cell formula E(d,e) = 1/(phi(de) ln(de)), d in [2,X], e in [3,X]: sum E and sum E^2/2.
 * phi(de) = phi(d) phi(e) g / phi(g), g = gcd(d,e).  usage: ./aud_taskf X */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
static long gcdl(long a,long b){ while(b){ long t=a%b; a=b; b=t; } return a; }
int main(int c,char**v){ long X=atol(v[1]); long *ph=malloc((X+1)*sizeof(long));
  for(long i=0;i<=X;i++) ph[i]=i; for(long i=2;i<=X;i++) if(ph[i]==i) for(long j=i;j<=X;j+=i) ph[j]-=ph[j]/i;
  double S=0,S2=0; for(long d=2;d<=X;d++) for(long e=3;e<=X;e++){ long g=gcdl(d,e); double f=(double)ph[d]*ph[e]*g/ph[g];
    double E=1.0/(f*log((double)d*e)); S+=E; S2+=E*E/2; }
  printf("X=%ld: sum E = %.4f, sum E^2/2 = %.5f\n",X,S,S2); return 0; }
