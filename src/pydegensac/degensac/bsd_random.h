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

#endif /* BSD_RANDOM_H */
