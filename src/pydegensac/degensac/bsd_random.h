#ifndef BSD_RANDOM_H
#define BSD_RANDOM_H

/* Thread-safe reimplementation of the 4.3BSD TYPE_3 random()/srandom() pair,
   bit-compatible with the macOS/FreeBSD libc sequence for every seed and
   with glibc for every nonzero seed. The RANSAC loops reseed the generator
   on every iteration (deterministic seed chain), which makes libc srandom()
   a major bottleneck; this local version removes the calling/locking overhead
   via precomputed warm-up transformation matrices and thread-local state.

   Also provides degensac_rand()/degensac_srand() for degeneracy checks:
   uses thread-local state matching macOS Park-Miller rand() on Darwin and
   glibc TYPE_3 rand() on Linux, guaranteeing thread-safety without inter-thread
   state corruption. */

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
