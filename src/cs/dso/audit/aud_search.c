/* Independent re-implementation (C) of the E2 exact search: every subset S (|S| = k <= KMAX) of usable census
 * primes inside one 2-adic class (v2 ord, v2 z), with an odd number of inert members, satisfying
 *     n == 1 mod odd(ord_q),  n == -1 mod odd(z_q)  for every q in S   (n = prod S).
 * Enumeration: anchor = member with the largest (M, p), M = odd(ord) odd(z); the product m of the others must be
 * == t (t == 1 mod odd ord_A, t == -s_A mod odd z_A). k = 5 uses a lower-pair / upper-pair split (each 4-set
 * {a<b<c<d} is met exactly once as (a,b) + (c,d)), not the dedup set of dso_analyze.py.
 * Control mode (CTRL=1): Lucas side dropped (z := 1), no parity: hits are strong base-2 psp.
 * usage: ./aud_search rows.txt D KMAX CTRL FILTER12B   (rows: "p o z s")
 * Prints, for k = 2..KMAX: subsets in scope, anchor survivors, full solutions; every full solution as "SOL ...".
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned long long u64; typedef unsigned __int128 u128;
typedef struct { u64 p, o, z; int s, a, b; u64 oo, oz, M, t; } P;
static u64 gcdu(u64 a,u64 b){ while(b){ u64 t=a%b; a=b; b=t; } return a; }
static int v2(u64 x){ return __builtin_ctzll(x); }
static u64 mulm(u64 a,u64 b,u64 m){ return (u64)((u128)a*b%m); }
static u64 inv(u64 a,u64 m){ if(m==1) return 0; long long t=0,nt=1; u64 r=m,nr=a%m;
  while(nr){ u64 q=r/nr; long long tt=t-(long long)q*nt; t=nt; nt=tt; u64 rr=r-q*nr; r=nr; nr=rr; }
  if(r!=1) return 0; if(t<0) t+=m; return (u64)t; }
static int cmpP(const void*x,const void*y){ const P*a=x,*b=y; if(a->M!=b->M) return a->M<b->M?-1:1; return a->p<b->p?-1:(a->p>b->p); }
static int CTRL; static double scope[6]; static long surv[6], full[6];
static void check(P*mem,int*idx,int k){   /* idx: k indices (last = anchor) */
  surv[k]++;
  int ni=0; for(int i=0;i<k;i++) if(mem[idx[i]].s==-1) ni++;
  if(!CTRL && ni%2==0) return;
  for(int i=0;i<k;i++){ P*q=&mem[idx[i]];
    if(q->oo>1){ u64 r=1; for(int j=0;j<k;j++) r=mulm(r,mem[idx[j]].p%q->oo,q->oo); if(r!=1) return; }
    if(!CTRL && q->oz>1){ u64 r=1; for(int j=0;j<k;j++) r=mulm(r,mem[idx[j]].p%q->oz,q->oz); if(r!=q->oz-1) return; } }
  full[k]++; printf("SOL k=%d", k); for(int i=0;i<k;i++) printf(" %llu", mem[idx[i]].p); printf("\n");
}
/* chained hash: key -> list of ints */
typedef struct { int *head, *next; u64 *key; int *val; int cap, n, mask; } H;
static void hinit(H*h,int nmax){ int c=1; while(c<2*nmax+2) c<<=1; h->cap=c; h->mask=c-1; h->head=malloc(c*sizeof(int)); memset(h->head,-1,c*sizeof(int));
  h->next=malloc((nmax+1)*sizeof(int)); h->key=malloc((nmax+1)*sizeof(u64)); h->val=malloc((nmax+1)*sizeof(int)); h->n=0; }
static void hfree(H*h){ free(h->head); free(h->next); free(h->key); free(h->val); }
static inline int hslot(H*h,u64 k){ return (int)((k*0x9E3779B97F4A7C15ULL)>>17 & h->mask); }
static void hput(H*h,u64 k,int v){ int s=hslot(h,k); h->key[h->n]=k; h->val[h->n]=v; h->next[h->n]=h->head[s]; h->head[s]=h->n; h->n++; }
int main(int argc,char**argv){
  FILE*f=fopen(argv[1],"r"); int D=atoi(argv[2]); int KMAX=atoi(argv[3]); CTRL=atoi(argv[4]); int F12=atoi(argv[5]);
  int cap=1<<16, n=0; P*all=malloc(cap*sizeof(P)); u64 p,o,z; int s;
  while(fscanf(f,"%llu %llu %llu %d",&p,&o,&z,&s)==4){
    if(s==0 || p==2) continue;
    if(CTRL){ z=1; s=-1; }
    else { u64 oo=o>>v2(o), oz=z>>v2(z); if(s==1 && gcdu(oo,oz)!=1) continue;
      if(F12 && D==5 && (o%5==0 || z%5==0)) continue; }
    P q; q.p=p; q.o=o; q.z=z; q.s=s; q.a=v2(o); q.b=v2(z); q.oo=o>>q.a; q.oz=z>>q.b; q.M=q.oo*q.oz;
    if(q.oo>1 && q.oz>1 && gcdu(q.oo,q.oz)!=1){ fprintf(stderr,"internal: not coprime\n"); return 1; }
    /* t == 1 mod oo, t == -s mod oz */
    u64 ms = (q.oz==1)?0:((s==1)? q.oz-1 : 1%q.oz);
    if(q.M==1) q.t=0; else { u64 x=1%q.oo; /* x + oo*y == ms mod oz */
      u64 y = (q.oz==1)?0: mulm((ms + q.oz - x%q.oz)%q.oz, inv(q.oo%q.oz,q.oz), q.oz); q.t=(x + (u64)((u128)q.oo*y % q.M))%q.M; 
      if(q.t%q.oo != 1%q.oo || q.t%q.oz != ms){ fprintf(stderr,"crt fail\n"); return 1; } }
    if(n==cap){ cap*=2; all=realloc(all,cap*sizeof(P)); } all[n++]=q; }
  fclose(f);
  printf("# usable %d (ctrl %d, filter12b %d)\n", n, CTRL, F12);
  /* group by class */
  int nc=0; int ca[4096], cb[4096];
  for(int i=0;i<n;i++){ int j; for(j=0;j<nc;j++) if(ca[j]==all[i].a && cb[j]==all[i].b) break; if(j==nc){ ca[nc]=all[i].a; cb[nc]=all[i].b; nc++; } }
  for(int c=0;c<nc;c++){
    int m=0; for(int i=0;i<n;i++) if(all[i].a==ca[c] && all[i].b==cb[c]) m++;
    P*mem=malloc(m*sizeof(P)); m=0; for(int i=0;i<n;i++) if(all[i].a==ca[c] && all[i].b==cb[c]) mem[m++]=all[i];
    qsort(mem,m,sizeof(P),cmpP);
    u64 *res=malloc(m*sizeof(u64)), *iv=malloc(m*sizeof(u64)); char*ok=malloc(m);
    for(int i=1;i<m;i++){
      P*A=&mem[i]; u64 M=A->M, t=A->t; int idx[6];
      for(int j=0;j<i;j++){ res[j]=mem[j].p%M; ok[j]=(gcdu(res[j]==0?M:res[j],M)==1) || M==1; iv[j]= ok[j]?inv(res[j],M):0; }
      double di=i; scope[2]+=di; scope[3]+=di*(di-1)/2; scope[4]+=di*(di-1)*(di-2)/6; scope[5]+=di*(di-1)*(di-2)*(di-3)/24;
      H h; hinit(&h,i); for(int j=0;j<i;j++) if(ok[j]) hput(&h,res[j],j);
      /* k = 2 */
      if(KMAX>=2) for(int e=h.head[hslot(&h,t)];e>=0;e=h.next[e]) if(h.key[e]==t){ idx[0]=h.val[e]; idx[1]=i; check(mem,idx,2); }
      /* k = 3 */
      if(KMAX>=3) for(int j=0;j<i;j++) if(ok[j]){ u64 need=mulm(t,iv[j],M)%(M?M:1); if(M==1) need=0;
        for(int e=h.head[hslot(&h,need)];e>=0;e=h.next[e]) if(h.key[e]==need && h.val[e]>j){ idx[0]=j; idx[1]=h.val[e]; idx[2]=i; check(mem,idx,3); } }
      /* k = 4 */
      if(KMAX>=4) for(int j=0;j<i;j++) if(ok[j]){ u64 b=mulm(t,iv[j],M); for(int l=j+1;l<i;l++) if(ok[l]){ u64 need=(M==1)?0:mulm(b,iv[l],M);
        for(int e=h.head[hslot(&h,need)];e>=0;e=h.next[e]) if(h.key[e]==need && h.val[e]>l){ idx[0]=j; idx[1]=l; idx[2]=h.val[e]; idx[3]=i; check(mem,idx,4); } } }
      /* k = 5 */
      if(KMAX>=5 && i>=4){ int np=0; for(int j=0;j<i;j++) if(ok[j]) for(int l=j+1;l<i;l++) if(ok[l]) np++;
        H hp; hinit(&hp,np); int *pa=malloc((np+1)*sizeof(int)), *pb=malloc((np+1)*sizeof(int)); int q=0;
        for(int j=0;j<i;j++) if(ok[j]) for(int l=j+1;l<i;l++) if(ok[l]){ pa[q]=j; pb[q]=l; hput(&hp,(M==1)?0:mulm(res[j],res[l],M),q); q++; }
        for(int x=0;x<np;x++){ int a=pa[x], b=pb[x]; u64 need=(M==1)?0:mulm(mulm(t,iv[a],M),iv[b],M);
          for(int e=hp.head[hslot(&hp,need)];e>=0;e=hp.next[e]) if(hp.key[e]==need){ int y=hp.val[e]; if(pa[y]>b){ idx[0]=a; idx[1]=b; idx[2]=pa[y]; idx[3]=pb[y]; idx[4]=i; check(mem,idx,5); } } }
        hfree(&hp); free(pa); free(pb); }
      hfree(&h);
    }
    free(res); free(iv); free(ok); free(mem);
  }
  for(int k=2;k<=KMAX;k++) printf("k = %d: subsets in scope %.4g, anchor survivors %ld, full solutions %ld\n", k, scope[k], surv[k], full[k]);
  return 0;
}
