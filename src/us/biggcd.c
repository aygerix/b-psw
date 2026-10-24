/* biggcd.c - g = gcd(2^W1 - 1, U_W2(1,Q)) for huge W1, W2 (last-prime step of Theorem B, Lemma 17(c)).
 * Mode "N": work modulo N = 2^W1 - 1 (fast reduction by folding), compute U_W2 mod N, then gcd(N, U).
 * Mode "U": compute U_W2 exactly by the doubling chain, then x = 2^W1 mod U (mpz_powm), g = gcd(U, x - 1).
 * Usage: ./biggcd N|U D W1 W2      prints g (decimal if < 400 digits) and its bit length.  */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
static mpz_t N, mask; static unsigned long W;
static void redN(mpz_t x){ /* x mod 2^W - 1, x >= 0 */
  static mpz_t hi; static int init = 0; if(!init){ mpz_init(hi); init = 1; }
  while(mpz_sizeinbase(x, 2) > W){ mpz_tdiv_q_2exp(hi, x, W); mpz_and(x, x, mask); mpz_add(x, x, hi); }
  if(mpz_cmp(x, N) >= 0) mpz_sub(x, x, N); }
static void red(mpz_t x, int modN, const mpz_t M){ if(modN){ if(mpz_sgn(x) < 0) mpz_mod(x, x, N); else redN(x); } }
/* Lucas U_k, V_k (P=1, Q) by binary doubling; if modN reduce mod N each step, else exact integers */
static void lucas(mpz_t U, mpz_t V, long Q, long D, const char *kstr, int modN){
  mpz_t k, Qk, t, u2, v2; mpz_inits(k, Qk, t, u2, v2, NULL); mpz_set_str(k, kstr, 10);
  mpz_set_ui(U, 0); mpz_set_ui(V, 2); mpz_set_ui(Qk, 1);
  long nb = mpz_sizeinbase(k, 2);
  for(long b = nb - 1; b >= 0; b--){
    mpz_mul(U, U, V); red(U, modN, N);                       /* U_2m = U_m V_m */
    mpz_mul(V, V, V); mpz_submul_ui(V, Qk, 2); mpz_addmul_ui(V, Qk, 0); /* V_2m = V^2 - 2 Q^m */
    red(V, modN, N);
    mpz_mul(Qk, Qk, Qk); red(Qk, modN, N);
    if(mpz_tstbit(k, b)){                                     /* (U,V) -> ((U+V)/2, (D U + V)/2) */
      mpz_add(u2, U, V); mpz_mul_si(v2, U, D); mpz_add(v2, v2, V);
      if(modN){ if(mpz_odd_p(u2)) mpz_add(u2, u2, N); if(mpz_odd_p(v2)) mpz_add(v2, v2, N); }
      mpz_tdiv_q_2exp(U, u2, 1); mpz_tdiv_q_2exp(V, v2, 1);
      if(modN){ red(U, 1, N); if(mpz_sgn(V) < 0) mpz_mod(V, V, N); else red(V, 1, N); }
      mpz_mul_si(Qk, Qk, Q); if(modN){ if(mpz_sgn(Qk) < 0) mpz_mod(Qk, Qk, N); else red(Qk, 1, N); }
    }
    if(modN && mpz_sgn(V) < 0) mpz_mod(V, V, N);
  }
  mpz_clears(k, Qk, t, u2, v2, NULL); }
int main(int argc, char **argv){
  if(argc < 5){ fprintf(stderr, "usage: %s N|U D W1 W2\n", argv[0]); return 1; }
  int modN = argv[1][0] == 'N'; long D = atol(argv[2]); long Q = (1 - D) / 4;
  mpz_t U, V, g, x, e; mpz_inits(U, V, g, x, e, N, mask, NULL);
  if(modN){
    W = strtoul(argv[3], 0, 10);
    mpz_set_ui(N, 0); mpz_setbit(N, W); mpz_sub_ui(N, N, 1); mpz_set(mask, N);
    lucas(U, V, Q, D, argv[4], 1);
    mpz_gcd(g, N, U);
  } else {
    lucas(U, V, Q, D, argv[4], 0); mpz_abs(U, U);
    fprintf(stderr, "U exact: %lu bits\n", (unsigned long)mpz_sizeinbase(U, 2));
    mpz_set_str(e, argv[3], 10); mpz_set_ui(x, 2); mpz_powm(x, x, e, U); mpz_sub_ui(x, x, 1);
    mpz_gcd(g, U, x);
  }
  size_t bits = mpz_sizeinbase(g, 2);
  if(mpz_sizeinbase(g, 10) < 400){ printf("g = "); mpz_out_str(stdout, 10, g); printf("\n"); }
  else printf("g has %zu bits (too long to print)\n", bits);
  printf("bits %zu\n", bits); return 0; }
