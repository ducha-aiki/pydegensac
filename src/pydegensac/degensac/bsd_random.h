#ifndef BSD_RANDOM_H
#define BSD_RANDOM_H

/* Thread-local RNGs preserving the estimator's platform-specific streams.
   On Unix, random()/srandom() use the historical BSD TYPE_3 sequence with
   precomputed warm-up matrices. The separate rand()/srand() stream matches
   Darwin's Park-Miller generator or glibc's TYPE_3 seed semantics.
   On Windows, both pairs share the same CRT-compatible 15-bit LCG stream. */

void degensac_srandom(unsigned seed);
long degensac_random(void);

void degensac_srand(unsigned seed);
int degensac_rand(void);

#if defined(_MSC_VER)
#  define DEGENSAC_THREAD_LOCAL __declspec(thread)
#elif defined(__STDC_VERSION__) && __STDC_VERSION__ >= 201112L
#  define DEGENSAC_THREAD_LOCAL _Thread_local
#elif defined(__GNUC__) || defined(__clang__)
#  define DEGENSAC_THREAD_LOCAL __thread
#else
#  define DEGENSAC_THREAD_LOCAL
#endif

#endif /* BSD_RANDOM_H */
