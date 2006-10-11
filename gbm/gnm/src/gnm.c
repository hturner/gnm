/* vector * matrix */

# include <Rinternals.h> /* for length */
# include <R_ext/Applic.h> /* for dgemm */

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

/* Computes elementwise products between submatrices of base matrix M and 
   columns of gradient matrix V, summing 'common' results an putting result 
   in submatrix of X. This version has start point in M, V and X for each
   "term" */
SEXP newsubprod(SEXP M, SEXP V, SEXP X, SEXP a, SEXP b, SEXP c, 
		  SEXP nt, SEXP lt, SEXP ls, SEXP nr, SEXP nc, SEXP max) {
  /* currently set up for single term so nt = 1 and all integers here */
  int i, j, k, l, start, end, common, nrow = INTEGER(nr)[0], 
    n = INTEGER(max)[0], final, jump;
  double *p[n], *q[n];

  final = INTEGER(nt)[0];
  for (i = 0; i < final; i++){
    common = INTEGER(nc)[i];
    jump = INTEGER(lt)[i];
    p[0] = &REAL(M)[INTEGER(a)[i]];
    q[0] = &REAL(V)[INTEGER(b)[i]];
    for (l = 1; l < common; l++) {
      p[l] = p[l - 1] + jump;
      q[l] = q[l - 1] + nrow;
    }
    start = INTEGER(c)[i];
    end = INTEGER(ls)[i];
    k = 0;
    for (j = start; j < end; j++, k = (++k == nrow) ? 0 : k) {
      REAL(X)[j] = *(p[0])++ * q[0][k];
      for (l = 1; l < common; l++){
	REAL(X)[j] += *(p[l])++ * q[l][k];
      }
    }
  }
  return(X);
}
/* computes single column of design matrix */
SEXP single(SEXP M, SEXP V, SEXP a, SEXP lt, SEXP nr, SEXP nc) {
  int j, k, l, nrow = INTEGER(nr)[0], common = INTEGER(nc)[0], jump;
  double *p[common], *q[common];
  SEXP col;

  jump = INTEGER(lt)[0];
  p[0] = &REAL(M)[INTEGER(a)[0]];
  q[0] = &REAL(V)[0];
  for (l = 1; l < common; l++) {
    p[l] = p[l - 1] + jump;
    q[l] = q[l - 1] + nrow;
  }
  k = 0;
  PROTECT(col = allocVector(REALSXP, nrow));
  for (j = 0; j < nrow; j++, k = (++k == nrow) ? 0 : k) {
    REAL(col)[j] = *(p[0])++ * q[0][k];
    for (l = 1; l < common; l++){
      REAL(col)[j] += *(p[l])++ * q[l][k];
    }
  }
  UNPROTECT(1);
  return(col);
}
