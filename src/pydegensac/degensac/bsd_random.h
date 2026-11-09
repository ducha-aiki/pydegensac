#ifndef BSD_RANDOM_H
#define BSD_RANDOM_H

/* Lock-free reimplementation of the 4.3BSD TYPE_3 random()/srandom() pair,
   bit-compatible with the macOS/FreeBSD libc sequence for every seed and
   with glibc for every nonzero seed.  The RANSAC loops reseed the generator
   on every iteration (deterministic seed chain), which makes libc srandom()
   -- 10*31 warm-up spins behind a lock -- the dominant cost of the whole
   estimator; this local version removes the calling/locking overhead while
   producing the identical stream. */

void degensac_srandom(unsigned seed);
long degensac_random(void);

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
