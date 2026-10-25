/* Independent audit of the E2 DSO census completeness (written separately from dso_census*.gp).
 * Uses Mobius inversion instead of polcyclo:
 *   A = prod_{d<=X} Phi_d(2)          = prod_{k<=X} (2^k - 1)^{M(floor(X/k))}
 *   B = prod_{e<=X} Phi_e(alpha,beta) = prod_{k<=X} |U_k(1,Q)|^{M(floor(X/k))}
 * (M = Mertens function; sum_{d<=X} sum_{k|d} mu(d/k) f(k) = sum_k f(k) M(X/k)).
 * G = gcd(A, B). Then every prime of the given list is divided out of G (all powers); prints whether every listed
 * prime divides G and whether the cofactor is 1 (=> prime support of G == list).
 * For Q = -1 (D = 5) U_k = F_k from mpz_fib_ui; otherwise the recurrence U_k = U_{k-1} - Q U_{k-2}.
 * usage: ./aud_census X Q primes.txt
 */
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>
static double now(void){ struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return t.tv_sec+1e-9*t.tv_nsec; }
static void prodtree(mpz_t *v, long n, mpz_t out){
  if(n==0){ mpz_set_ui(out,1); return; }
  while(n>1){ long m=0; for(long i=0;i+1<n;i+=2){ mpz_mul(v[m],v[i],v[i+1]); m++; }
    if(n&1){ mpz_swap(v[m],v[n-1]); m++; } for(long i=m;i<n;i++) mpz_set_ui(v[i],1); n=m; }
  mpz_swap(out,v[0]);
}
int main(int argc,char**argv){
  long X=atol(argv[1]); long Q=atol(argv[2]); const char*fn=argv[3];
  double t0=now();
  int *mu=malloc((X+1)*sizeof(int)); char *cp=calloc(X+1,1); long *M=malloc((X+1)*sizeof(long));
  for(long i=0;i<=X;i++) mu[i]=1; mu[0]=0;
  for(long i=2;i<=X;i++) if(!cp[i]){ for(long j=i;j<=X;j+=i){ if(j>i) cp[j]=1; mu[j]=-mu[j]; }
    for(long j=i*i;j<=X;j+=i*i) mu[j]=0; }
  M[0]=0; for(long i=1;i<=X;i++) M[i]=M[i-1]+mu[i];
  mpz_t *num=malloc(X*sizeof(mpz_t)), *den=malloc(X*sizeof(mpz_t)); long nn=0,nd=0;
  mpz_t f,A,B,T,Um1,U,Un; mpz_inits(f,A,B,T,Um1,U,Un,NULL);
  /* A */
  for(long k=1;k<=X;k++){ long m=M[X/k]; if(!m) continue;
    mpz_set_ui(f,0); mpz_setbit(f,k); mpz_sub_ui(f,f,1);
    if(m>0){ mpz_init(num[nn]); mpz_pow_ui(num[nn],f,m); nn++; } else { mpz_init(den[nd]); mpz_pow_ui(den[nd],f,-m); nd++; } }
  prodtree(num,nn,A); prodtree(den,nd,T);
  if(!mpz_divisible_p(A,T)){ printf("ERROR: A numerator not divisible by denominator\n"); return 1; }
  mpz_divexact(A,A,T);
  for(long i=0;i<nn;i++) mpz_clear(num[i]); for(long i=0;i<nd;i++) mpz_clear(den[i]); nn=nd=0;
  printf("A: %zu bits (%.1f s)\n", mpz_sizeinbase(A,2), now()-t0); fflush(stdout);
  /* B */
  mpz_set_ui(Um1,0); mpz_set_ui(U,1);      /* U_0, U_1 */
  for(long k=1;k<=X;k++){
    if(Q==-1) mpz_fib_ui(f,k); else mpz_set(f,U);
    if(mpz_sgn(f)==0){ printf("ERROR: degenerate U_%ld = 0\n",k); return 1; }
    long m=M[X/k];
    if(m){ mpz_abs(T,f); if(m>0){ mpz_init(num[nn]); mpz_pow_ui(num[nn],T,m); nn++; } else { mpz_init(den[nd]); mpz_pow_ui(den[nd],T,-m); nd++; } }
    if(Q!=-1){ mpz_mul_si(Un,Um1,-Q); mpz_add(Un,Un,U); mpz_swap(Um1,U); mpz_swap(U,Un); } /* U_{k+1} = U_k - Q U_{k-1} */
  }
  prodtree(num,nn,B); prodtree(den,nd,T);
  if(!mpz_divisible_p(B,T)){ printf("ERROR: B numerator not divisible by denominator\n"); return 1; }
  mpz_divexact(B,B,T);
  for(long i=0;i<nn;i++) mpz_clear(num[i]); for(long i=0;i<nd;i++) mpz_clear(den[i]);
  printf("B: %zu bits (%.1f s)\n", mpz_sizeinbase(B,2), now()-t0); fflush(stdout);
  mpz_t G; mpz_init(G); mpz_gcd(G,A,B); mpz_clear(A); mpz_clear(B);
  printf("G = gcd(A,B): %zu bits (%.1f s)\n", mpz_sizeinbase(G,2), now()-t0); fflush(stdout);
  FILE*fp=fopen(fn,"r"); unsigned long p; long cnt=0, notdiv=0, maxe=0; unsigned long pmax=0;
  while(fscanf(fp,"%lu",&p)==1){ cnt++; long e=0; while(mpz_divisible_ui_p(G,p)){ mpz_divexact_ui(G,G,p); e++; }
    if(e==0){ notdiv++; printf("listed prime %lu does NOT divide G\n",p); } if(e>maxe) maxe=e; if(p>pmax) pmax=p; }
  fclose(fp);
  printf("listed primes %ld (max %lu), not dividing G: %ld, max exponent %ld\n",cnt,pmax,notdiv,maxe);
  if(mpz_cmp_ui(G,1)==0) printf("COFACTOR = 1: prime support of G equals the list\n");
  else { printf("COFACTOR != 1: %zu bits left", mpz_sizeinbase(G,2)); if(mpz_sizeinbase(G,2)<200) gmp_printf(" = %Zd", G); printf("\n"); }
  printf("total %.1f s\n", now()-t0);
  return 0;
}
