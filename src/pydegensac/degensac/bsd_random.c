#include "bsd_random.h"

#include <stdint.h>

/* TYPE_3 additive-feedback generator: x[i] = (x[i-3] + x[i-31]) mod 2^32,
   output is the top 31 bits.  Matches FreeBSD/macOS libc random(). */

#define RAND_DEG 31
#define RAND_SEP 3
#define WARMUP  (10 * RAND_DEG)

static uint32_t rand_tbl[RAND_DEG];
static int fptr_idx = RAND_SEP;
static int rptr_idx = 0;

/* libc srandom() warms the table with 310 discarded draws.  The update is
   linear mod 2^32, so the warmed table equals a fixed 31x31 coefficient
   matrix times the seeded words -- 961 independent multiply-adds the
   compiler can vectorize, instead of a 310-step serial dependency chain. */
static uint32_t warm_mat[RAND_DEG][RAND_DEG];
static uint32_t pow16807[RAND_DEG];  /* 16807^i mod 2^31-1 */
static int tables_ready = 0;

static void init_tables(void)
{
    static uint32_t v[RAND_DEG + WARMUP][RAND_DEG];
    uint64_t p;
    int n, k;

    for (k = 0; k < RAND_DEG; k++)
        v[k][k] = 1;
    for (n = RAND_DEG; n < RAND_DEG + WARMUP; n++)
        for (k = 0; k < RAND_DEG; k++)
            v[n][k] = v[n - RAND_SEP][k] + v[n - RAND_DEG][k];
    /* table slot (RAND_SEP+n)%RAND_DEG receives its last write at warm-up
       step n, for n in the final RAND_DEG steps */
    for (n = WARMUP - RAND_DEG; n < WARMUP; n++)
        for (k = 0; k < RAND_DEG; k++)
            warm_mat[(RAND_SEP + n) % RAND_DEG][k] = v[RAND_DEG + n][k];

    p = 1;
    for (n = 0; n < RAND_DEG; n++) {
        pow16807[n] = (uint32_t) p;
        p = (p * 16807) % 2147483647u;
    }
    tables_ready = 1;
}

/* Park-Miller step via Schrage's trick; FreeBSD substitutes 123459876 for a
   zero input so the LCG cannot stick at zero. */
static int32_t good_rand(int32_t x)
{
    int32_t hi, lo;

    if (x == 0)
        x = 123459876;
    hi = x / 127773;
    lo = x % 127773;
    x = 16807 * lo - 2836 * hi;
    if (x < 0)
        x += 0x7fffffff;
    return x;
}

long degensac_random(void)
{
    uint32_t val;

    rand_tbl[fptr_idx] += rand_tbl[rptr_idx];
    val = rand_tbl[fptr_idx] >> 1;
    if (++fptr_idx >= RAND_DEG)
        fptr_idx = 0;
    if (++rptr_idx >= RAND_DEG)
        rptr_idx = 0;
    return (long) val;
}

void degensac_srandom(unsigned seed)
{
    uint32_t init[RAND_DEG], x[RAND_DEG];
    int i, j, k;

    if (!tables_ready)
        init_tables();

    if ((int32_t) seed > 0 && seed < 2147483647u) {
        /* seeding LCG is Park-Miller, so init[i] = 16807^i * seed mod
           2^31-1: independent per word (vectorizable) instead of a 30-step
           serial chain; reduction by Mersenne-prime folding */
        for (i = 0; i < RAND_DEG; i++) {
            uint64_t m = (uint64_t) seed * pow16807[i];
            uint64_t v = (m & 0x7fffffff) + (m >> 31);
            v = (v & 0x7fffffff) + (v >> 31);
            if (v == 2147483647u)  /* wrap; v stays nonzero for seed > 0 */
                v = 0;
            init[i] = (uint32_t) v;
        }
    } else {
        /* seed 0 or with the top bit set: libc runs the LCG on the raw
           int32 (with the good_rand zero substitution) -- follow it */
        init[0] = (uint32_t) seed;
        for (i = 1; i < RAND_DEG; i++)
            init[i] = (uint32_t) good_rand((int32_t) init[i - 1]);
    }
    /* the warm-up recurrence visits the seeded slots rotated by RAND_SEP */
    for (k = 0; k < RAND_DEG; k++)
        x[k] = init[(k + RAND_SEP) % RAND_DEG];
    for (j = 0; j < RAND_DEG; j++) {
        uint32_t acc = 0;
        for (k = 0; k < RAND_DEG; k++)
            acc += warm_mat[j][k] * x[k];
        rand_tbl[j] = acc;
    }
    fptr_idx = RAND_SEP;
    rptr_idx = 0;
}
