/* Library of tools wrapping LAPACK utilities and making their usage a bit more comfortable. */
#include <stdlib.h>
#include <stddef.h>

#include "lapwrap.h"

/* DGESVD prototype (LAPACK) */
#ifdef _WIN32
extern void dgesvd_( char* jobu, char* jobvt, lapack_int* m, lapack_int* n, double* a,
                    lapack_int* lda, double* s, double* u, lapack_int* ldu, double* vt, lapack_int* ldvt,
                    double* work, lapack_int* lwork, lapack_int* info );
#endif

#ifdef __linux__
extern void dgesvd_( char* jobu, char* jobvt, lapack_int* m, lapack_int* n, double* a,
                    lapack_int* lda, double* s, double* u, lapack_int* ldu, double* vt, lapack_int* ldvt,
                    double* work, lapack_int* lwork, lapack_int* info );
#endif

/* Standard (=FULL) SVD */
int lap_SVD (double *d, double *a, double *u, lapack_int m, double *vt, lapack_int n) {
  lapack_int lda = m, ldu = m, ldvt = n, info = 1, lwork;
  double wkopt;
  double *work;
  /* Query and allocate the optimal workspace */
  lwork = -1;
#ifdef _WIN32
  dgesvd_( "All", "All", &m, &n, a, &lda, d, u, &ldu, vt, &ldvt, &wkopt, &lwork, &info );
#endif
#ifdef __linux__
  dgesvd_( "All", "All", &m, &n, a, &lda, d, u, &ldu, vt, &ldvt, &wkopt, &lwork, &info );
#endif

  lwork = (lapack_int) wkopt;
  work = (double *) malloc ( lwork*sizeof(double) );
  /* Compute SVD */
#ifdef _WIN32
  dgesvd_( "All", "All", &m, &n, a, &lda, d, u, &ldu, vt, &ldvt, work, &lwork, &info );
#endif

#ifdef __linux__
  dgesvd_( "All", "All", &m, &n, a, &lda, d, u, &ldu, vt, &ldvt, work, &lwork, &info );
#endif
  free(work);
  if (info != 0) {
      return 1;
    } else {
      return 0;
    }
}


/* DSYEV prototype */
#ifdef _WIN32
extern void dsyev_( char* jobz, char* uplo, lapack_int* n, double* a, lapack_int* lda,
		   double* w, double* work, lapack_int* lwork, lapack_int* info );

#endif

#ifdef __linux__
extern void dsyev_( char* jobz, char* uplo, lapack_int* n, double* a, lapack_int* lda,
		   double* w, double* work, lapack_int* lwork, lapack_int* info );
#endif


/* Eigen-decomposition */
int lap_eig(double *a, double *ev, lapack_int n) {
  lapack_int lda = n, info = 1, lwork;
  double wkopt;
  double *work;
  /* Query and allocate the optimal workspace */
  lwork = -1;
#ifdef _WIN32
  dsyev_( "Vectors", "Upper", &n, a, &lda, ev, &wkopt, &lwork, &info );
#endif

#ifdef __linux__
  dsyev_( "Vectors", "Upper", &n, a, &lda, ev, &wkopt, &lwork, &info );
#endif

  lwork = (lapack_int) wkopt;
  work = (double *) malloc ( lwork*sizeof(double) );
  /* Solve eigenproblem */
#ifdef _WIN32
  dsyev_( "Vectors", "Upper", &n, a, &lda, ev, work, &lwork, &info );
#endif

#ifdef __linux__
  dsyev_( "Vectors", "Upper", &n, a, &lda, ev, work, &lwork, &info );
#endif
  free(work);
  if( info != 0 ) {
      return 1;
    } else {
      return 0;
    }
}
