/* Library of tools wrapping LAPACK utilities and making their usage a bit more comfortable. */
#include <stdlib.h>
#include <stddef.h>

#include "lapwrap.h"

/* Workspace handling.

   Both wrappers below used to do a workspace query (lwork = -1), malloc the
   result, factorise, and free -- on every call.  These are RANSAC inner-loop
   calls on 3x3 and 9x9 matrices, where that bookkeeping costs more than the
   factorisation itself.

   Two changes, both output-preserving:

   - the workspace lives on the stack (LAP_WORK_STACK doubles, 4 KB) whenever
     it fits, so the malloc/free pair disappears;
   - the queried lwork is memoised per matrix shape, so the query call happens
     once per shape instead of once per call.

   LAPACK's optimal lwork is a pure function of the shape and the LAPACK build,
   so replaying a memoised value is exactly what the query would have returned.
   That matters: lwork selects the blocked-vs-unblocked path inside LAPACK, so
   passing a *different* value would change the arithmetic. Passing the same
   value with a larger backing buffer does not -- only work[0..lwork-1] is
   touched.

   The cache is thread-local: no locking, no sharing, and being a plain array
   it cannot leak when a thread exits. */
#define LAP_WORK_STACK 512
#define LAP_CACHE_SLOTS 8

#if defined(_MSC_VER)
#  define LAP_THREAD_LOCAL __declspec(thread)
#elif defined(__STDC_VERSION__) && __STDC_VERSION__ >= 201112L && !defined(__STDC_NO_THREADS__)
#  define LAP_THREAD_LOCAL _Thread_local
#elif defined(__GNUC__)
#  define LAP_THREAD_LOCAL __thread
#else
#  define LAP_THREAD_LOCAL /* single-threaded fallback: re-queries per call */
#endif

typedef struct {
    lapack_int m, n, lwork;
} lap_work_entry;

/* Look up the memoised lwork for a shape, or -1 if this shape is new. */
static lapack_int lap_cache_get(const lap_work_entry *cache, lapack_int m, lapack_int n) {
    int i;
    for (i = 0; i < LAP_CACHE_SLOTS; i++) {
        if (cache[i].lwork > 0 && cache[i].m == m && cache[i].n == n) {
            return cache[i].lwork;
        }
    }
    return -1;
}

/* Record an lwork for a shape. The table is tiny and the call sites use only a
   couple of shapes; if it ever fills, later shapes simply keep re-querying. */
static void lap_cache_put(lap_work_entry *cache, lapack_int m, lapack_int n, lapack_int lwork) {
    int i;
    for (i = 0; i < LAP_CACHE_SLOTS; i++) {
        if (cache[i].lwork <= 0) {
            cache[i].m = m;
            cache[i].n = n;
            cache[i].lwork = lwork;
            return;
        }
    }
}

/* DGESVD prototype (LAPACK).

   Declared unconditionally. It used to be declared, and called, only under
   _WIN32 or __linux__ -- so on macOS, which defines neither, every LAPACK call
   in this file was preprocessed away: lap_SVD/lap_eig returned without
   computing anything. See the note on lap_eig below. */
extern void dgesvd_( char* jobu, char* jobvt, lapack_int* m, lapack_int* n, double* a,
                    lapack_int* lda, double* s, double* u, lapack_int* ldu, double* vt, lapack_int* ldvt,
                    double* work, lapack_int* lwork, lapack_int* info );

/* Standard (=FULL) SVD */
int lap_SVD (double *d, double *a, double *u, lapack_int m, double *vt, lapack_int n) {
  static LAP_THREAD_LOCAL lap_work_entry cache[LAP_CACHE_SLOTS];
  double stack_work[LAP_WORK_STACK];
  lapack_int lda = m, ldu = m, ldvt = n, info = 1, lwork;
  double wkopt;
  double *work;
  double *heap_work = NULL;

  /* Query the optimal workspace, unless this shape has been seen before */
  lwork = lap_cache_get(cache, m, n);
  if (lwork <= 0) {
    lwork = -1;
    dgesvd_( "All", "All", &m, &n, a, &lda, d, u, &ldu, vt, &ldvt, &wkopt, &lwork, &info );
    lwork = (lapack_int) wkopt;
    lap_cache_put(cache, m, n, lwork);
  }

  if (lwork <= LAP_WORK_STACK) {
    work = stack_work;
  } else {
    heap_work = (double *) malloc ( lwork*sizeof(double) );
    if (heap_work == NULL) {
        return 1;
      }
    work = heap_work;
  }

  /* Compute SVD */
  dgesvd_( "All", "All", &m, &n, a, &lda, d, u, &ldu, vt, &ldvt, work, &lwork, &info );

  free(heap_work);
  if (info != 0) {
      return 1;
    } else {
      return 0;
    }
}


/* DSYEV prototype. Declared unconditionally -- see dgesvd_ above. */
extern void dsyev_( char* jobz, char* uplo, lapack_int* n, double* a, lapack_int* lda,
		   double* w, double* work, lapack_int* lwork, lapack_int* info );


/* Eigen-decomposition */
int lap_eig(double *a, double *ev, lapack_int n) {
  static LAP_THREAD_LOCAL lap_work_entry cache[LAP_CACHE_SLOTS];
  double stack_work[LAP_WORK_STACK];
  lapack_int lda = n, info = 1, lwork;
  double wkopt;
  double *work;
  double *heap_work = NULL;

  /* Query the optimal workspace, unless this shape has been seen before */
  lwork = lap_cache_get(cache, n, n);
  if (lwork <= 0) {
    lwork = -1;
    dsyev_( "Vectors", "Upper", &n, a, &lda, ev, &wkopt, &lwork, &info );
    lwork = (lapack_int) wkopt;
    lap_cache_put(cache, n, n, lwork);
  }

  if (lwork <= LAP_WORK_STACK) {
    work = stack_work;
  } else {
    heap_work = (double *) malloc ( lwork*sizeof(double) );
    if (heap_work == NULL) {
        return 1;
      }
    work = heap_work;
  }

  /* Solve eigenproblem */
  dsyev_( "Vectors", "Upper", &n, a, &lda, ev, work, &lwork, &info );

  free(heap_work);
  if( info != 0 ) {
      return 1;
    } else {
      return 0;
    }
}
