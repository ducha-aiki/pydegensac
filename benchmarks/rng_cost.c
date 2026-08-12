/* Measure the per-reseed cost this branch removes, on whatever platform you
   run it on.

   The benchmark showed a 1.19x/1.16x speed-up on Linux/glibc against 4.2x/2.2x
   on the macOS golden pairs. Two explanations are in play: the libc being
   replaced (macOS srandom() re-derives the TYPE_3 state behind a lock, glibc
   does not) and the ISA (arm64 vs x86-64 vectorise the 31x31 warm-up product
   differently). This separates them: it times exactly the call pair the RANSAC
   loops used to make per iteration, libc against the local implementation, so
   a run on each machine says which one moved.

   Build and run from benchmarks/:

     cc -O3 -I../src/pydegensac/degensac rng_cost.c \
        ../src/pydegensac/degensac/bsd_random.c -o rng_cost && ./rng_cost
*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "bsd_random.h"

#define N_RESEED 200000

static double now_ns(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec * 1e9 + ts.tv_nsec;
}

/* Draws per RANSAC iteration: the minimal sample size, 7 for F, 4 for H. */
static void report(const char *problem, int n_draw)
{
    volatile long sink = 0;
    double t0, before_ns, after_ns;
    unsigned i;
    int k;

    /* Warm up both paths so first-call effects don't land in the timing. */
    for (i = 1; i <= 1000; i++) {
        srandom(i);
        sink += random();
        degensac_srandom(i);
        sink += degensac_random();
    }

    /* Before: reseed via libc, then draw the minimal sample via libc. */
    t0 = now_ns();
    for (i = 1; i <= N_RESEED; i++) {
        srandom(i);
        for (k = 0; k < n_draw; k++)
            sink += random();
    }
    before_ns = (now_ns() - t0) / N_RESEED;

    /* After: no reseed at all (the generator is keyed once per estimator
       call), just the draws, from the local generator. */
    degensac_srandom(1);
    t0 = now_ns();
    for (i = 1; i <= N_RESEED; i++) {
        for (k = 0; k < n_draw; k++)
            sink += degensac_random();
    }
    after_ns = (now_ns() - t0) / N_RESEED;

    printf("%s (%d draws/iteration)\n", problem, n_draw);
    printf("  before: libc srandom + %d x random   %8.1f ns/iteration\n",
           n_draw, before_ns);
    printf("  after:  %d x degensac_random         %8.1f ns/iteration\n",
           n_draw, after_ns);
    printf("  saved                                %8.1f ns  (%.1fx)\n",
           before_ns - after_ns, before_ns / after_ns);
    printf("  predicted wall-clock saving: %.1f ms at 10k iterations, "
           "%.1f ms at 50k\n\n",
           (before_ns - after_ns) * 10000 / 1e6,
           (before_ns - after_ns) * 50000 / 1e6);
    (void)sink;
}

int main(void)
{
    printf("Per-iteration RNG cost removed by this branch\n"
           "(the loops used to re-key the generator every iteration; now it is\n"
           " keyed once per estimator call, and the generator is local)\n\n");
    report("F", 7);
    report("H", 4);
    printf("Compare across machines: if the saving per iteration is much\n"
           "larger on macOS than on glibc Linux, the libc explains the gap;\n"
           "if it is similar, the difference is problem size or ISA.\n");
    return 0;
}
