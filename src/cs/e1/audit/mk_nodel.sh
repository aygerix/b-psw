# bbc variant WITHOUT the maximal-closure deletion lemma: children inherit the parent's alive set (minus the primes
# excluded by the side fix) instead of the maximal optimal closure.  Everything else identical.
sed 's/ch.alive = malloc(nP); memcpy(ch.alive, selmax, nP);/ch.alive = malloc(nP); memcpy(ch.alive, nd.alive, nP);/' ../bbc3.c > bbc_nodel.c
grep -c "memcpy(ch.alive, nd.alive, nP)" bbc_nodel.c
cp ../bbc3.c bbc_audit.c
gcc -O2 -o bbc_audit bbc_audit.c -lm && gcc -O2 -o bbc_nodel bbc_nodel.c -lm && echo built
