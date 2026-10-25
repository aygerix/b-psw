// Streaming version of sunit_pairs.c (constant memory): DFS over all S-smooth x <= X, test x+1 by trial
// division over S. Complete for the given X. Usage: ./sunit_stream X q1 ... qk  (prints x; stats on stderr)
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef unsigned __int128 u128;
static uint64_t X, P[256]; static int k; static uint64_t cnt, pairs;
static int smooth(uint64_t v){ for(int i=k-1;i>=0 && v>1;i--){ while(v%P[i]==0) v/=P[i]; } return v==1; }
static void dfs(int i, uint64_t v){
  if(i==k){ cnt++; if(smooth(v+1)){ printf("%llu\n",(unsigned long long)v); pairs++; } return; }
  for(u128 w=v; w<=X; w*=P[i]) dfs(i+1,(uint64_t)w);
}
int main(int argc,char**argv){
  X=strtoull(argv[1],0,10); k=argc-2;
  for(int i=0;i<k;i++) P[i]=strtoull(argv[i+2],0,10);
  for(int i=0;i<k;i++) for(int j=i+1;j<k;j++) if(P[j]>P[i]){uint64_t t=P[i];P[i]=P[j];P[j]=t;}  // descending
  dfs(0,1);
  fprintf(stderr,"S-smooth <= %llu: %llu ; pairs: %llu\n",(unsigned long long)X,(unsigned long long)cnt,(unsigned long long)pairs);
  return 0;
}
