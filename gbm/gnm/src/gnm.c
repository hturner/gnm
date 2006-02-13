/* vector * matrix */

# include <Rinternals.h> /* for length */
# include <R_ext/Applic.h> /* for dgemm */
# include <R_ext/Lapack.h> /* for dgelsy */

/* for dgelsy */
# ifndef max
#    define max(a, b) ((a > b)? a:b)
# endif

/* copied from src/main/array.c */
static void matprod(double *x, int nrx, int ncx,
		    double *y, int nry, int ncy, double *z)
{
    char *transa = "N", *transb = "N";
    int i,  j, k;
    double one = 1.0, zero = 0.0, sum;
    Rboolean have_na = FALSE;

    if (nrx > 0 && ncx > 0 && nry > 0 && ncy > 0) {
	/* Don't trust the BLAS to handle NA/NaNs correctly: PR#4582
	 * The test is only O(n) here
	 */
	for (i = 0; i < nrx*ncx; i++)
	    if (ISNAN(x[i])) {have_na = TRUE; break;}
	if (!have_na)
	    for (i = 0; i < nry*ncy; i++)
		if (ISNAN(y[i])) {have_na = TRUE; break;}
	if (have_na) {
	    for (i = 0; i < nrx; i++)
		for (k = 0; k < ncy; k++) {
		    sum = 0.0;
		    for (j = 0; j < ncx; j++)
			sum += x[i + j * nrx] * y[j + k * nry];
		    z[i + k * nrx] = sum;
		}
	} else
	    F77_CALL(dgemm)(transa, transb, &nrx, &ncy, &ncx, &one,
			    x, &nrx, y, &nry, &zero, z, &nrx);
    } else /* zero-extent operations should return zeroes */
	for(i = 0; i < nrx*ncy; i++) z[i] = 0;
}

/* computes matrix product between submatrix of M and vector v */
SEXP submatprod(SEXP M, SEXP v, SEXP am, SEXP nr, SEXP nc) {
  R_len_t  a = INTEGER(am)[0], nrm = INTEGER(nr)[0], ncm = INTEGER(nc)[0];
  SEXP ans;

  PROTECT(ans = allocVector(REALSXP, nrm));
  matprod(REAL(M) + a, nrm, ncm,
	  REAL(v), ncm, 1, REAL(ans));
  UNPROTECT(1);

  return(ans);
}

/* computes elementwise product between submatrix of M and vector v
   then puts result in submatrix of X */
SEXP subprod(SEXP X, SEXP M, SEXP v, SEXP a, SEXP z, SEXP nv) {
  R_len_t i = INTEGER(a)[0], j = 0, 
    last = INTEGER(z)[0], len_v = INTEGER(nv)[0];

  for ( ; i <= last; j = (++j == len_v) ? 0 : j) { 
    REAL(X)[i] = REAL(M)[i] * REAL(v)[j];
    i++;
  }
  return(X);
}

/* put results of nonlin localDesignFunction in X */

SEXP nonlin(SEXP X, SEXP a, SEXP z, SEXP expr, SEXP rho) {
  R_len_t i = INTEGER(a)[0], i1 = 0, last = INTEGER(z)[0];
  SEXP ans;

  PROTECT(ans = coerceVector(eval(expr, rho), REALSXP));
  for ( ; i <= last;) {
    REAL(X)[i++] = REAL(ans)[i1++];
  }
  UNPROTECT(1);
  return(X);
}

/* solves Ax = b using Fortran routine dgelsy */
void dgelsy(int *m, int *n, int *nrhs, double *a, double *b, double *rcond, 
	    int *rank, double *work, int *lwork, double *ans) {
  int i, j, ldb = max(*m, *n), jpvt[*n], info = 1;

  F77_CALL(dgelsy)(m, n, nrhs, a, m, b, &ldb, jpvt, rcond, rank, work, lwork, 
		   &info);

  for(i = 0; i < *n; i++) {
    for(j = 0; j < *nrhs; j++)
      ans[i + *n * j] = b[i + *m * j];
  }
}

