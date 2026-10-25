# bbc_diag: the no-deletion variant, instrumented to count at every branching node how many alive primes lie OUTSIDE the
# maximal optimal closure (the primes the deletion lemma would remove).  If this is > 0 while the B&B tree (node count,
# M_ub) is identical to the deleting version, the lemma's prediction (descendant optimal closures unaffected) is confirmed.
sed -e 's/static int SB = 0;/static int SB = 0; static long long ndel_tot = 0, nbr_tot = 0, nbr_del = 0;/' \
    -e 's|      // branch on a row used on both sides by the maximal closure; weight = selected primes per side|      { long long d = 0; for (int i = 0; i < nP; i++) if (nd.alive[i] \&\& !selmax[i]) d++; ndel_tot += d; nbr_tot++; if (d) nbr_del++; }\n      // branch on a row|' \
    -e 's|    sbcl = 0;|    sbcl = 0; printf("  DIAG branching nodes %lld, with deletable primes %lld, total deletable %lld\\n", nbr_tot, nbr_del, ndel_tot); ndel_tot = nbr_tot = nbr_del = 0;|' \
    bbc_nodel.c > bbc_diag.c
grep -c DIAG bbc_diag.c; grep -c ndel_tot bbc_diag.c
gcc -O2 -o bbc_diag bbc_diag.c -lm && echo built
