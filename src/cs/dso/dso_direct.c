/* Independent check of the E2 census (D = 5) for primes p <= P0: compute ord_p(2) and the Fibonacci entry point
 * z(p) directly for every prime p <= P0 (sieve + trial-division factorisation of p-1, p+1), print the primes with
 * ord <= X and z <= X as "p ord z (5/p)".
 * Usage: ./dso_direct P0 X > direct.txt
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
typedef unsigned long long u64; typedef unsigned __int128 u128;
static u64 mulm(u64 a, u64 b, u64 m){ return (u64)((u128)a*b % m); }
static u64 powm(u64 b, u64 e, u64 m){ u64 r=1%m; b%=m; while(e){ if(e&1) r=mulm(r,b,m); b=mulm(b,b,m); e>>=1; } return r; }
/* F_k mod m by fast doubling: (F_k, F_{k+1}) */
static void fib(u64 k, u64 m, u64 *f0, u64 *f1){
  u64 a=0, b=1%m; /* F_0, F_1 */
  for(int i=63;i>=0;i--){
    /* (a,b)=(F_j,F_{j+1}) -> (F_2j, F_2j+1) */
    u64 t = (2*b % m + m - a) % m;        /* 2F_{j+1} - F_j */
    u64 c = mulm(a, t, m);                 /* F_2j */
    u64 d = (mulm(a,a,m) + mulm(b,b,m)) % m; /* F_2j+1 */
    a=c; b=d;
    if((k>>i)&1){ u64 e2=(a+b)%m; a=b; b=e2; }
  }
  *f0=a; *f1=b;
}
static u64 fibm(u64 k, u64 m){ u64 a,b; fib(k,m,&a,&b); return a; }
static int *smallp; static int nsp;
static int factor(u64 n, u64 *q){ int k=0; for(int i=0;i<nsp && (u64)smallp[i]*smallp[i]<=n;i++){ u64 p=smallp[i]; if(n%p==0){ q[k++]=p; while(n%p==0) n/=p; } } if(n>1) q[k++]=n; return k; }
int main(int argc, char **argv){
  u64 P0 = strtoull(argv[1],0,10); u64 X = strtoull(argv[2],0,10);
  /* odd sieve */
  u64 nb = P0/2+1; unsigned char *comp = calloc(nb,1);
  for(u64 i=3;i*i<=P0;i+=2) if(!comp[i/2]) for(u64 j=i*i;j<=P0;j+=2*i) comp[j/2]=1;
  u64 lim=1; while(lim*lim<=P0+1) lim++;
  smallp = malloc(sizeof(int)*(lim+10)); nsp=0; smallp[nsp++]=2;
  for(u64 i=3;i<=lim;i+=2) if(!comp[i/2]) smallp[nsp++]=(int)i;
  u64 q[64]; long cnt=0, nprimes=0;
  for(u64 p=3;p<=P0;p+=2){
    if(comp[p/2]) continue; nprimes++;
    if(p==5){ printf("5 4 5 0\n"); cnt++; continue; }
    /* ord_p(2) */
    u64 o=p-1; int k=factor(p-1,q);
    for(int i=0;i<k;i++) while(o%q[i]==0 && powm(2,o/q[i],p)==1) o/=q[i];
    if(o>X) continue;
    int r=(int)(p%5); int kr=(r==1||r==4)?1:-1;
    u64 z = (kr==1)? p-1 : p+1; k=factor(z,q);
    if(fibm(z,p)!=0){ fprintf(stderr,"fib check failed p=%llu\n",p); return 1; }
    for(int i=0;i<k;i++) while(z%q[i]==0 && fibm(z/q[i],p)==0) z/=q[i];
    if(z>X) continue;
    printf("%llu %llu %llu %d\n", p, o, z, kr); cnt++;
  }
  fprintf(stderr,"primes 3..P0: %ld, DSO (ord, z <= %llu): %ld\n", nprimes, X, cnt);
  return 0;
}
