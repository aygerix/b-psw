// Smooth neighbors (x, x+1 both z-smooth) via Conrey-Holmstrom-McLaughlin closure:
// rule1: (s-r) | r(r+1)  -> t = r + r(r+1)/(s-r)        [ (r/(r+1))*((s+1)/s) = t/(t+1) ]
// rule2: (r+s+1) | r*s   -> t = r*s/(r+s+1)             [ ((r+1)/r)*((s+1)/s) = (t+1)/t ]
// Pairwise over current set until closure. u128 arithmetic with overflow guard (values < 2^63).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned __int128 u128;
static uint64_t *S; static size_t nS, capS;
static uint64_t *H; static size_t Hcap;   // hash set, 0 = empty
static int overflow = 0;
static size_t hsh(uint64_t x){ x ^= x>>33; x*=0xff51afd7ed558ccdULL; x^=x>>33; x*=0xc4ceb9fe1a85ec53ULL; x^=x>>33; return x & (Hcap-1);}
static int hhas(uint64_t x){ size_t i=hsh(x); while(H[i]){ if(H[i]==x) return 1; i=(i+1)&(Hcap-1);} return 0;}
static void rehash(void){ size_t oc=Hcap; uint64_t*oh=H; Hcap*=2; H=calloc(Hcap,8); for(size_t i=0;i<oc;i++) if(oh[i]){ size_t j=hsh(oh[i]); while(H[j]) j=(j+1)&(Hcap-1); H[j]=oh[i]; } free(oh);}
static int hadd(uint64_t x){ if(hhas(x)) return 0; if(nS*2>=Hcap) rehash(); size_t i=hsh(x); while(H[i]) i=(i+1)&(Hcap-1); H[i]=x;
  if(nS==capS){capS*=2; S=realloc(S,capS*8);} S[nS++]=x; return 1;}
static int primes[2000], np;
static int smooth(u128 x){ for(int i=0;i<np && x>1;i++){ while(x%primes[i]==0) x/=primes[i]; } return x==1; }
int main(int argc,char**argv){
  int z=atoi(argv[1]);
  for(int p=2;p<=z;p++){int ok=1; for(int q=2;q*q<=p;q++) if(p%q==0){ok=0;break;} if(ok) primes[np++]=p;}
  capS=1<<16; S=malloc(capS*8); Hcap=1<<17; H=calloc(Hcap,8);
  for(uint64_t r=1;r<(uint64_t)z;r++) hadd(r);
  size_t done=0; int pass=0;
  while(1){
    size_t before=nS; pass++;
    for(size_t i=0;i<before;i++){
      for(size_t j=(i<done?done:0); j<before; j++){ // pairs with at least one new element
        if(j<=i && j>=done) { if(j==i) continue; }
        if(i<done && j<done) continue;
        uint64_t a=S[i], b=S[j]; if(a==b) continue;
        uint64_t r = a<b?a:b, s = a<b?b:a;
        u128 rr = (u128)r*(r+1);
        uint64_t d = s-r;
        if(rr % d == 0){ u128 t = (u128)r + rr/d; if(t >> 63) overflow++; else hadd((uint64_t)t); }
        u128 u = (u128)r + s + 1;
        if(((u128)r*s) % u == 0){ u128 t=((u128)r*s)/u; if(t>0) hadd((uint64_t)t); }
      }
    }
    done=before;
    fprintf(stderr,"z=%d pass %d: |S|=%zu new=%zu overflow=%d\n",z,pass,nS,nS-before,overflow);
    if(nS==before) break;
  }
  // verify all smooth & output
  size_t bad=0; for(size_t i=0;i<nS;i++){ if(!smooth(S[i])||!smooth((u128)S[i]+1)) bad++; printf("%llu\n",(unsigned long long)S[i]); }
  fprintf(stderr,"total %zu, nonsmooth %zu, overflow %d\n",nS,bad,overflow);
  return 0;
}
