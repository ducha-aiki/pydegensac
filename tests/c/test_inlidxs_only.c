/* Exact regression gate for the index-only LO scorer. */
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "rtools.h"

#define SENTINEL (-1515870811)

static uint64_t state = UINT64_C(0x7d7a4b39f1e2d3c5);

static double uniform(void)
{
    state = state * UINT64_C(6364136223846793005) + 1;
    return ((state >> 11) * 0x1.0p-53) * 20.0 - 10.0;
}

static void check_case(const char *name, const double *err, int len, double th)
{
    Score full;
    unsigned only;
    int full_inl[259], only_inl[259];
    double saved[257];
    int i, expected = 0;

    assert(len >= 0 && len <= 257);
    memcpy(saved, err, (size_t)len * sizeof(*err));
    for (i = 0; i < len + 2; ++i)
        full_inl[i] = only_inl[i] = SENTINEL;

    full = inlidxs(err, len, th, full_inl);
    only = inlidxs_only(err, len, th, only_inl);
    if (full.I != only) {
        fprintf(stderr, "%s: count differs: %u != %u\n", name, full.I, only);
        assert(0);
    }
    assert(memcmp(full_inl, only_inl, (size_t)len * sizeof(*full_inl)) == 0);
    assert(full_inl[len] == SENTINEL && full_inl[len + 1] == SENTINEL);
    assert(only_inl[len] == SENTINEL && only_inl[len + 1] == SENTINEL);
    assert(memcmp(err, saved, (size_t)len * sizeof(*err)) == 0);

    for (i = 0; i < len; ++i) {
        if (err[i] <= th) {
            assert(expected < (int)only);
            assert(only_inl[expected++] == i);
        }
    }
    assert(expected == (int)only);
}

int main(void)
{
    const int lengths[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
                           13, 14, 15, 16, 17, 65, 257};
    double err[257];
    int i, trial, li;

    for (li = 0; li < (int)(sizeof(lengths) / sizeof(lengths[0])); ++li) {
        int len = lengths[li];
        for (i = 0; i < len; ++i)
            err[i] = (i % 3 == 0) ? -1.0 : (i % 3 == 1) ? 0.0 : 1.0;
        check_case("mixed", err, len, 0.0);
        check_case("negative-threshold", err, len, -0.0);
        check_case("infinite-threshold", err, len, INFINITY);
        check_case("nan-threshold", err, len, NAN);
        for (trial = 0; trial < 8; ++trial) {
            for (i = 0; i < len; ++i)
                err[i] = uniform();
            check_case("random", err, len, uniform());
        }
    }

    for (i = 0; i < 257; ++i)
        err[i] = (i & 1) ? nextafter(1.0, -INFINITY) : nextafter(1.0, INFINITY);
    err[0] = -0.0;
    err[1] = 0.0;
    err[2] = NAN;
    err[3] = INFINITY;
    err[4] = -INFINITY;
    err[5] = 1.0;
    check_case("boundaries", err, 257, 1.0);
    check_case("negative", err, 257, -1.0);
    check_case("all", err, 257, INFINITY);
    for (i = 0; i < 257; ++i)
        err[i] = -1.0;
    check_case("every-index", err, 257, 0.0);
    puts("inlidxs_only matches inlidxs index output");
    return 0;
}
