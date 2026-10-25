/* AUDIT COPY of dso_heur.c with a mixed-cell accumulator added (lines marked MIXED); the model itself is unchanged.
 * usage as dso_heur; a 7th argument stops after the mixed sum. */
/* E2 heuristic: expected number E(d,e,s) of primes p with ord_p(2) = d, z(p) = e, (5/p) = s, for D = 5.
 *
 * Model (independent random elements, exact 2-adic rules):
 *   P(ord_p 2 = d | p) = 2 phi(d)/(p-1) if  [2 is a QR (p = +-1 mod 8) and d | (p-1)/2]
 *                                        or [2 is a QNR and d | p-1 and v2(d) = v2(p-1)];  else 0.
 *   gamma = alpha/beta = -alpha^2 has order z(p).
 *     inert (s=-1): gamma in the norm-1 torus T (order p+1); gamma is a square in T iff p = 1 mod 4.
 *     split (s=+1): gamma in F_p^*; gamma is a square iff -1 is, iff p = 1 mod 4.
 *     So with N = p - s:  P(z = e | p) = 2 phi(e)/N if [p = 1 mod 4 and e | N/2] or [p = 3 mod 4 and e | N and
 *     v2(e) = v2(N)];  else 0.
 *   p ranges over x = c + kL (L = lcm(d,e), c = CRT(1 mod d, s mod e)), Pr(x prime) = (L'/phi(L'))/ln x with
 *   L' = lcm(L, 2) (0 for even x), times the (5/x) = s factor (1/2 if 5 does not divide L, else 0/1).
 *   Exact terms for k < K, tail = naive integral times the measured refined/naive ratio.
 * Output: aggregates of E, E^2/2, E^3/6 by category, and E for each occupied cell listed in CELLS file.
 * Usage: ./dso_heur X cells.txt [K] [c]  > heur.out   (cells.txt: lines "d e s"; 'c' = only the listed cells)
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
typedef long long ll; typedef __int128 i128;
static int *phi;
static ll gcdl(ll a, ll b){ while(b){ ll t=a%b; a=b; b=t;} return a<0?-a:a; }
static ll inv_mod(ll a, ll m){ ll g=m, x=0, x1=1, a1=a%m; if(a1<0)a1+=m; ll b=a1; ll old_r=b, r=m, old_s=1, s=0;
  while(r){ ll q=old_r/r, t=old_r-q*r; old_r=r; r=t; t=old_s-q*s; old_s=s; s=t; }
  (void)g;(void)x;(void)x1; if(old_r!=1) return -1; old_s%=m; if(old_s<0) old_s+=m; return old_s; }
static int v2l(ll x){ return __builtin_ctzll(x); }
static __attribute__((unused)) double phil(ll n){ /* phi for n up to X^2 via factoring with small trial division */
  double r=n; ll m=n; for(ll q=2;q*q<=m;q++) if(m%q==0){ r=r/q*(q-1); while(m%q==0) m/=q; } if(m>1) r=r/m*(m-1); return r; }
/* refined weight A(x) B(x) at a specific x (x odd, x = 1 mod d, x = s mod e) */
static double wref(ll x, ll d, ll e, int s){
  double A, B;
  int qr2 = (x%8==1 || x%8==7);
  ll xm = x-1; int vx = v2l(xm), vd = v2l(d);
  if(qr2) A = (vd <= vx-1) ? 2.0*phi[d]/(double)xm : 0.0;
  else    A = (vd == vx) ? 2.0*phi[d]/(double)xm : 0.0;
  ll N = (s==1) ? x-1 : x+1; int vN = v2l(N), ve = v2l(e);
  if(x%4==1) B = (ve <= vN-1) ? 2.0*phi[e]/(double)N : 0.0;
  else       B = (ve == vN) ? 2.0*phi[e]/(double)N : 0.0;
  return A*B;
}
static double wnaive(ll x, ll d, ll e, int s){ ll N=(s==1)?x-1:x+1; return (double)phi[d]*phi[e]/((double)(x-1)*(double)N); }
static int kron5(ll x){ int r=(int)(x%5); return (r==1||r==4)?1:(r==0?0:-1); }
/* E(d,e,s); returns -1 if incompatible */
static double G_small, G_Q, G_C; static ll G_X;   /* per-cell side outputs */
static double Ecell(ll d, ll e, int s, int K){
  G_small = 0; G_Q = 0; G_C = 0;
  ll g = gcdl(d,e);
  if(s==-1 && (2 % g)!=0) return -1;
  if(s==-1 && (d%5==0 || e%5==0)) return 0;   /* inert p cannot be +-1 mod 5 */
  ll L = d/g*e;
  /* c = 1 mod d, c = s mod e: c = 1 + d t, d t = s-1 mod e -> (d/g) t = (s-1)/g mod e/g */
  ll c;
  if(s==1) c = 1 % L; else { ll eg=e/g, dg=d/g; ll rhs=((s-1)/g)%eg; if(rhs<0) rhs+=eg; ll t=0;
    if(eg>1){ ll iv=inv_mod(dg%eg,eg); t=(ll)((i128)rhs*iv%eg); } c=(1+d*t)%L; }
  int five = (L%5==0);
  double phiL = (double)phi[d]*phi[e]/phi[g];               /* phi(lcm) phi(gcd) = phi(d) phi(e) */
  double dens = (L%2==0) ? (double)L/phiL : 2.0*L/phiL;      /* L' = lcm(L,2) */
  /* p | Phi_d(2) and p | Phi_e(alpha,beta): p < 2^phi(d)+1, p < 1.62^phi(e)+2 (crude) */
  double lnU = fmin(phi[d]*log(2.0)+1.0, phi[e]*log(1.6180339887)+1.0);
  double sum=0, sref=0, snai=0; ll x = c; int k;
  /* only the random part: x > X (every prime <= X is in the census deterministically) */
  if(x <= G_X) x += ((G_X - x)/L + 1)*L;
  double dens5 = five ? dens : dens*1.25;          /* condition on 5 not dividing x when 5 does not divide L */
  for(k=0;k<K;k++, x+=L){
    if(log((double)x) > lnU) break;
    if(x%2==0) continue;
    int k5 = kron5(x);
    if(k5==0 || k5!=s) continue;                   /* exact sign of (5/x) */
    double pr = dens5/log((double)x);
    double wr = wref(x,d,e,s), wn = wnaive(x,d,e,s);
    double pi_k = pr*wr; if(pi_k > 1) pi_k = 1;
    sum += pi_k; sref += wr; snai += wn;
    G_Q += pi_k*pi_k; G_C += pi_k*pi_k*pi_k;
  }
  /* tail factor: average of refined/naive over the 2-adic classes of x in the progression.
     refined = naive * 4 * alpha(x) beta(x) with alpha, beta in {0,1} determined by x mod 2^(m+3), m = v2(L). */
  int m2 = (L%2==0) ? v2l(L) : 0; ll M2 = 1LL << (m2+3); int ncls = 0, ngood = 0;
  for(ll r = c % (1LL<<m2); r < M2; r += (1LL<<m2)){
    if(r%2==0) continue;
    ll xr = r + 8*M2;   /* representative with the same residue mod 2^(m+3), large enough */
    ncls++; if(wref(xr,d,e,s) > 0) ngood++;
  }
  double ratio = ncls ? 4.0*ngood/ncls : 0.0; (void)sref; (void)snai;
  /* naive tail: primes in the class, density 1/(phi(L) ln x) per unit x; phi(d)phi(e)/x^2 weight */
  double w5avg = five ? ((kron5(c)==s)?1.0:0.0) : 0.5;
  if(five && kron5(c)==0) w5avg = 0;
  double xt = (double)x, tail = 0;
  if(log(xt) < lnU){ tail = 1.0/(xt*log(xt)); if(lnU < 700) tail -= 1.0/(exp(lnU)*lnU); if(tail<0) tail=0;
    tail *= w5avg*(double)phi[d]*phi[e]/phiL; }
  return sum + ratio*tail;
}
int main(int argc, char **argv){
  int X = atoi(argv[1]); int K = argc>3 ? atoi(argv[3]) : 16;
  phi = malloc((X+1)*sizeof(int));
  for(int i=0;i<=X;i++) phi[i]=i;
  for(int i=2;i<=X;i++) if(phi[i]==i) for(int j=i;j<=X;j+=i) phi[j]-=phi[j]/i;
  /* categories: s (0=inert,1=split) x shared-odd (gcd(odd d, odd e) > 1) x ratio class */
  /* ratio class: (d/g, e/g) with both <= 4 -> 1 + (d/g-1)*4 + (e/g-1), else 0 */
  static double S1[2][2][17], S2[2][2][17], S3[2][2][17];
  /* also by size bin of L: floor(log2 L) */
  static double B1[2][64], B2[2][64], SM[2], EB[2][64];
  G_X = X;
  int cellsonly = (argc>4 && argv[4][0]=='c');
  /* band mode: argv[5] = X0 < X: prime cut x > X0, and separate sums over cells with max(d,e) > X0
     (= expected new primes / 2-cells / 3-cells in the band, beyond a complete census at X0) */
  ll X0 = (argc>5) ? atoll(argv[5]) : X; if(X0 < X) G_X = X0;
  static double BA1[2], BA2[2], BA3[2];
  double MIX=0, MIXB=0;
  for(ll d=2; d<=X; d++) for(ll e=3; e<=X; e++){ ll od=d>>v2l(d), oe=e>>v2l(e);
    if(gcdl(od,oe)!=1 || d%5==0 || e%5==0) continue;
    double Ei=Ecell(d,e,-1,K); if(Ei<=0) continue; double Es=Ecell(d,e,1,K); if(Es<=0) continue;
    MIX+=Ei*Es; if(d > X0 || e > X0) MIXB+=Ei*Es; }
  printf("MIXED cells (1 inert + 1 usable split, same (d,e)), d,e <= %d, primes > %lld: expected %.4e; band max(d,e) > %lld: %.4e\n", X, (long long)G_X, MIX, (long long)X0, MIXB);
  if(argc>6) return 0;
  for(ll d=2; d<=(cellsonly?1:X); d++) for(ll e=3; e<=X; e++) for(int si=0; si<2; si++){
    int s = si ? 1 : -1;
    double E = Ecell(d,e,s,K);
    if(E<=0) continue;
    ll g = gcdl(d,e); ll od=d>>v2l(d), oe=e>>v2l(e); int sh = gcdl(od,oe)>1;
    int rc = (d/g<=4 && e/g<=4) ? 1+(int)((d/g-1)*4+(e/g-1)) : 0;
    /* P(>=2) and P(>=3) for a sum of independent Bernoulli(pi_k): e2 = (S^2 - Q)/2, e3 = (S^3 - 3 S Q + 2 C)/6 */
    double e2 = (E*E - G_Q)/2, e3 = (E*E*E - 3*E*G_Q + 2*G_C)/6; if(e2<0) e2=0; if(e3<0) e3=0;
    S1[si][sh][rc]+=E; S2[si][sh][rc]+=e2; S3[si][sh][rc]+=e3; SM[si]+=G_small;
    if(d > X0 || e > X0){ BA1[si]+=E; BA2[si]+=e2; BA3[si]+=e3; }
    int lb = (int)floor(log2((double)(d/g*e))); B1[si][lb]+=E; B2[si][lb]+=e2;
    { int eb=(int)floor(log10(E))+40; if(eb<0) eb=0; if(eb>63) eb=63; EB[si][eb]+=E; }
  }
  const char *nm[2]={"inert","split"};
  double T1=0,T2=0,T3=0;
  for(int si=0;si<2;si++) for(int sh=0;sh<2;sh++) for(int rc=0;rc<17;rc++){
    if(S1[si][sh][rc]==0) continue;
    T1+=S1[si][sh][rc]; T2+=S2[si][sh][rc]; T3+=S3[si][sh][rc];
    if(rc==0) printf("CAT %s shared=%d ratio=generic  sumE=%.4f  E2cells=%.5f  E3cells=%.3e\n", nm[si], sh, S1[si][sh][rc], S2[si][sh][rc], S3[si][sh][rc]);
    else printf("CAT %s shared=%d ratio=%d:%d  sumE=%.4f  E2cells=%.5f  E3cells=%.3e\n", nm[si], sh, (rc-1)/4+1, (rc-1)%4+1, S1[si][sh][rc], S2[si][sh][rc], S3[si][sh][rc]);
  }
  printf("TOTAL sumE=%.3f E2cells=%.4f E3cells=%.4e\n", T1, T2, T3);
  if(X0 < X) for(int si=0;si<2;si++) printf("BAND %s cells with max(d,e) in (%lld, %d]: E=%.4f E2cells=%.4e E3cells=%.4e\n", nm[si], X0, X, BA1[si], BA2[si], BA3[si]);
  for(int si=0;si<2;si++){ double t=0; for(int sh=0;sh<2;sh++) for(int rc=0;rc<17;rc++) t+=S1[si][sh][rc];
    printf("SIGN %s expected primes p>X: %.3f\n", nm[si], t); }
  for(int si=0;si<2;si++) for(int b=0;b<64;b++) if(B1[si][b]>0) printf("LBIN %s log2L=%d sumE=%.4f E2cells=%.5f\n", nm[si], b, B1[si][b], B2[si][b]);
  for(int si=0;si<2;si++) for(int b=0;b<64;b++) if(EB[si][b]>0) printf("EBIN %s log10E=%d sumE=%.5f\n", nm[si], b-40, EB[si][b]);
  if(argc>2){ FILE *f=fopen(argv[2],"r"); ll d,e; int s;
    while(f && fscanf(f,"%lld %lld %d",&d,&e,&s)==3){ double E=Ecell(d,e,s,K); printf("CELL %lld %lld %d E=%.6g Q=%.6g\n",d,e,s,E,G_Q); }
    if(f) fclose(f); }
  return 0;
}
