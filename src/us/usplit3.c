/* usplit.c - scan primes p < X and report "usable split" primes for Selfridge discriminants D:
 *   (D/p) = +1, p does not divide D*Q, and gcd(odd(ord_p 2), odd(z_D(p))) = 1,
 * where z_D(p) is the rank of apparition of U_k(1, Q), Q = (1-D)/4.  For split p both orders divide
 * p - 1, so usability forces odd(ord)*odd(z) | odd(p - 1): every split prime of a BPSW psp is of this kind.
 * Method: segmented sieve gives the largest prime factor R (exponent e) of p - 1. Filter: R must fail to
 * divide ord or z, i.e. 2^((p-1)/R^e) == 1 or U_((p-1)/R^e) == 0 mod p. Survivors get the exact test.
 * Usage: ./usplit X D1 D2 ...  or  ./usplit X0:X D1 D2 ... (scan [X0, X); output: D p ord z)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
typedef unsigned __int128 u128;
typedef uint64_t u64;

static u64 mulmod(u64 a, u64 b, u64 m){ return (u64)((u128)a*b % m); }
static u64 powmod(u64 a, u64 e, u64 m){ u64 r=1%m; a%=m; while(e){ if(e&1) r=mulmod(r,a,m); a=mulmod(a,a,m); e>>=1;} return r; }
static int jacobi(int64_t a, u64 n){ /* n odd > 0 */
  int64_t aa = a % (int64_t)n; if(aa<0) aa += n; u64 x=(u64)aa, y=n; int r=1;
  while(x){ while(!(x&1)){ x>>=1; u64 t=y&7; if(t==3||t==5) r=-r; } u64 t=x; x=y; y=t; if((x&3)==3 && (y&3)==3) r=-r; x%=y; }
  return y==1 ? r : 0; }
/* U_k(1,Q) mod p via (U,V) doubling; returns U_k */
static u64 lucasU(u64 k, int64_t Q, u64 p){
  int64_t D = 1 - 4*Q; u64 Dm = (u64)(((D % (int64_t)p) + (int64_t)p) % (int64_t)p);
  u64 Qm = (u64)(((Q % (int64_t)p) + (int64_t)p) % (int64_t)p);
  u64 U=0, V=2, Qk=1, inv2=(p+1)/2;
  if(k==0) return 0;
  int b=63; while(!((k>>b)&1)) b--;
  for(; b>=0; b--){
    U = mulmod(U,V,p); V = (mulmod(V,V,p) + p - mulmod(2,Qk,p)) % p; Qk = mulmod(Qk,Qk,p);
    if((k>>b)&1){ u64 U1 = mulmod((U+V)%p, inv2, p); u64 V1 = mulmod((mulmod(Dm,U,p)+V)%p, inv2, p); U=U1; V=V1; Qk=mulmod(Qk,Qm,p); }
  }
  return U; }

/* Montgomery arithmetic mod odd p < 2^63 */
typedef struct { u64 p, pinv, r1, r2; } mctx;
static inline u64 redc(u128 t, const mctx *c){ u64 m = (u64)t * c->pinv; u128 u = (t + (u128)m * c->p) >> 64; u64 r=(u64)u; return r >= c->p ? r - c->p : r; }
static inline u64 mmul(u64 a, u64 b, const mctx *c){ return redc((u128)a*b, c); }
static void minit(mctx *c, u64 p){ c->p=p; u64 inv=p; for(int i=0;i<6;i++) inv*=2-p*inv; c->pinv=-inv;
  c->r1 = (u64)(((u128)1<<64) % p); c->r2 = (u64)((u128)c->r1 * c->r1 % p); }
static inline u64 tom(u64 a, const mctx *c){ return mmul(a % c->p, c->r2, c); }
static inline u64 madd(u64 a, u64 b, u64 p){ u64 s=a+b; return (s>=p||s<a)? s-p : s; }
static inline u64 msub(u64 a, u64 b, u64 p){ return a>=b ? a-b : a+p-b; }
static int mpow_is1(u64 base, u64 e, const mctx *c){ u64 x=tom(base,c), r=c->r1; while(e){ if(e&1) r=mmul(r,x,c); x=mmul(x,x,c); e>>=1; } return r==c->r1; }
static int mlucas_is0(u64 k, int64_t Q, const mctx *c){ u64 p=c->p;
  int64_t D = 1 - 4*Q; u64 Dm = tom((u64)(((D % (int64_t)p) + (int64_t)p) % (int64_t)p), c);
  u64 Qm = tom((u64)(((Q % (int64_t)p) + (int64_t)p) % (int64_t)p), c), inv2 = tom((p+1)/2, c);
  u64 U=0, V=madd(c->r1,c->r1,p), Qk=c->r1; if(k==0) return 1;
  int b=63; while(!((k>>b)&1)) b--;
  for(; b>=0; b--){ U = mmul(U,V,c); u64 t=mmul(V,V,c); V = msub(t, madd(Qk,Qk,p), p); Qk = mmul(Qk,Qk,c);
    if((k>>b)&1){ u64 U1 = mmul(madd(U,V,p), inv2, c); u64 V1 = mmul(madd(mmul(Dm,U,c),V,p), inv2, c); U=U1; V=V1; Qk=mmul(Qk,Qm,c); } }
  return U==0; }

static u64 oddpart(u64 x){ while(x && !(x&1)) x>>=1; return x; }
static u64 gcd64(u64 a, u64 b){ while(b){ u64 t=a%b; a=b; b=t;} return a; }
/* exact order of 2 and rank z, given factorization of m = p-1 (primes f[], exps g[]) */
static u64 ord2(u64 p, u64 *f, int *g, int nf){ u64 o=p-1; for(int i=0;i<nf;i++) for(int j=0;j<g[i];j++){ if(powmod(2,o/f[i],p)==1) o/=f[i]; else break; } return o; }
static u64 zrank(u64 p, int64_t Q, u64 *f, int *g, int nf){ u64 o=p-1; for(int i=0;i<nf;i++) for(int j=0;j<g[i];j++){ if(lucasU(o/f[i],Q,p)==0) o/=f[i]; else break; } return o; }

int main(int argc, char **argv){
  if(argc<3){ fprintf(stderr,"usage: %s X D1 [D2 ...]\n",argv[0]); return 1; }
  u64 X0 = 3, X; char *colon = strchr(argv[1], ':');
  if(colon){ X0 = strtoull(argv[1],0,10); X = strtoull(colon+1,0,10); if(X0 < 3) X0 = 3; } else X = strtoull(argv[1],0,10);
  int nD = argc-2; int64_t Ds[64], Qs[64]; long cnt[64]; memset(cnt,0,sizeof cnt);
  for(int i=0;i<nD;i++){ Ds[i]=atoll(argv[i+2]); Qs[i]=(1-Ds[i])/4; }
  u64 sq = (u64)sqrtl((long double)X) + 2;
  /* base primes up to sq */
  char *bs = calloc(sq+1,1); u64 *bp = malloc(sizeof(u64)*(sq/2+10)); int nb=0;
  for(u64 i=2;i<=sq;i++){ if(!bs[i]){ bp[nb++]=i; for(u64 j=i*i;j<=sq;j+=i) bs[j]=1; } }
  const u64 S = 1u<<22;
  u64 *rem = malloc(sizeof(u64)*S); u64 *lpf = malloc(sizeof(u64)*S); unsigned char *le = malloc(S); char *comp = malloc(S);
  long nprimes=0, nsurv=0;
  for(u64 lo=X0; lo<X; lo+=S){
    u64 hi = lo+S; if(hi>X) hi=X; u64 len = hi-lo;
    memset(comp,0,len);
    for(int i=0;i<nb;i++){ u64 r=bp[i]; if(r*r>=hi) break; u64 st = ((lo + r - 1)/r)*r; if(st < r*r) st = r*r; for(u64 j=st;j<hi;j+=r) comp[j-lo]=1; }
    /* factor v = n - 1 for n in [lo,hi): v in [lo-1, hi-1) */
    for(u64 i=0;i<len;i++){ rem[i]=lo+i-1; lpf[i]=1; le[i]=0; }
    for(int i=0;i<nb;i++){ u64 r=bp[i]; if(r>hi) break; u64 v0 = lo-1; u64 st = ((v0 + r - 1)/r)*r;
      for(u64 v=st; v<hi-1; v+=r){ u64 k=v-v0; if(comp[k+1-1]) continue; unsigned char c=0; while(rem[k]%r==0){ rem[k]/=r; c++; } lpf[k]=r; le[k]=c; } }
    for(u64 i=0;i<len;i++){
      u64 p = lo+i; if(comp[i] || p<3 || !(p&1)) continue; nprimes++;
      u64 R, e; if(rem[i]>1){ R=rem[i]; e=1; } else { R=lpf[i]; e=le[i]; }
      u64 Re=1; for(u64 t=0;t<e;t++) Re*=R;
      mctx mc; minit(&mc,p); int f2 = (R==2) ? 1 : mpow_is1(2,(p-1)/Re,&mc);
      for(int d=0; d<nD; d++){
        int64_t D=Ds[d], Q=Qs[d];
        if((u64)llabs(D)%p==0 || (u64)llabs(Q)%p==0) continue;
        if(jacobi(D,p)!=1) continue;
        if(R!=2 && !f2 && !mlucas_is0((p-1)/Re,Q,&mc)) continue;
        nsurv++;
        /* exact check: factor p-1 by trial division with base primes */
        u64 m=p-1, f[20]; int g[20], nf=0;
        for(int j=0;j<nb && bp[j]*bp[j]<=m;j++){ u64 r=bp[j]; if(m%r==0){ int c=0; while(m%r==0){m/=r;c++;} f[nf]=r; g[nf]=c; nf++; } }
        if(m>1){ f[nf]=m; g[nf]=1; nf++; }
        u64 o = ord2(p,f,g,nf), z = zrank(p,Q,f,g,nf);
        if(gcd64(oddpart(o),oddpart(z))==1){ cnt[d]++; printf("%lld %llu %llu %llu\n",(long long)D,(unsigned long long)p,(unsigned long long)o,(unsigned long long)z); fflush(stdout); }
      }
    }
  }
  fprintf(stderr,"X=%llu primes %ld survivors %ld usable:",(unsigned long long)X,nprimes,nsurv);
  for(int d=0;d<nD;d++) fprintf(stderr," D=%lld:%ld",(long long)Ds[d],cnt[d]); fprintf(stderr,"\n");
  return 0; }
