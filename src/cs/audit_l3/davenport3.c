/* Brute-force Davenport constant of (Z/3)^r for r = 1, 2, 3: max length of a zero-sum-free sequence.
   Enumerates multisets (nondecreasing sequences of element codes) of length L and checks all nonempty
   subsets via a reachable-sums bitmask.  Prints for each r: whether a zero-sum-free sequence of length L exists. */
#include <stdio.h>
#include <string.h>
#include <stdint.h>
static int r, q, L;            /* q = 3^r */
static int addt[27][27];
static long long count_free;
static int seq[16];
/* reach[d] : bitmask over group elements of sums of nonempty subsets of seq[0..d-1] */
static uint32_t reach[17];
static uint32_t shift_by(uint32_t m, int g) { uint32_t o = 0; for (int x = 0; x < q; x++) if (m >> x & 1) o |= 1u << addt[x][g]; return o; }
static void rec(int d, int lo) {
    if (d == L) { count_free++; return; }
    for (int g = lo; g < q; g++) {
        if (g == 0) continue;                    /* zero element is itself a zero-sum */
        uint32_t nr = reach[d] | (1u << g) | shift_by(reach[d], g);
        if (nr & 1u) continue;                   /* a zero-sum appears */
        reach[d + 1] = nr; seq[d] = g;
        rec(d + 1, g);
    }
}
int main(void) {
    for (r = 1; r <= 3; r++) {
        q = 1; for (int i = 0; i < r; i++) q *= 3;
        for (int a = 0; a < q; a++) for (int b = 0; b < q; b++) {
            int s = 0, pw = 1, x = a, y = b;
            for (int i = 0; i < r; i++) { s += ((x % 3 + y % 3) % 3) * pw; x /= 3; y /= 3; pw *= 3; }
            addt[a][b] = s;
        }
        for (L = 1; L <= 2 * r + 2; L++) {
            count_free = 0; reach[0] = 0; rec(0, 1);
            printf("r=%d L=%d zero-sum-free multisets: %lld\n", r, L, count_free);
            if (!count_free) { printf("r=%d  Davenport constant D = %d (2r+1 = %d)\n", r, L, 2 * r + 1); break; }
        }
    }
    return 0;
}
