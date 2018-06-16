/*  Copyright (C) 2005, 2006, 2008-2010, 2017 Heather Turner */
/*  */
/*  This program is free software; you can redistribute it and/or modify */
/*  it under the terms of the GNU General Public License as published by */
/*  the Free Software Foundation; either version 2 or 3 of the License */
/*  (at your option). */
/*   */
/*  This program is distributed in the hope that it will be useful, */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/*  GNU General Public License for more details. */
/*   */
/*  A copy of the GNU General Public License is available at */
/*  http://www.r-project.org/Licenses/ */


/* vector * matrix */

# include <Rinternals.h> /* for length */
# include <R_ext/Applic.h> /* for dgemm */
# include <R_ext/Rdynload.h> /* for registering routines */

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
  double *dX, *dM, *dv;
  dX = REAL(X);
  dM = REAL(M);
  dv = REAL(v);
  for ( ; i <= last; j = (++j == len_v) ? 0 : j) { 
    dX[i] = dM[i] * dv[j];
    i++;
  }
  return(X);
}

/* Computes elementwise products between submatrices of base matrix M and 
   columns of gradient matrix V, summing 'common' results an putting result 
   in submatrix of X. This version has start point in M, V and X for each
   "term" */
SEXP newsubprod(SEXP M, SEXP V, SEXP X, SEXP a, SEXP b, SEXP c, 
		  SEXP nt, SEXP lt, SEXP ls, SEXP nr, SEXP nc, SEXP max) {
  /* currently set up for single term so nt = 1 and all integers here */
  int i, j, k, l, *start, *end, *common, nrow = INTEGER(nr)[0], 
	  n = INTEGER(max)[0], final = INTEGER(nt)[0], *jump, *ia, *ib;
  double *p[n], *q[n], *dM, *dV, *dX;

  dM = REAL(M);
  dV = REAL(V);
  dX = REAL(X);
  start = INTEGER(c);
  end = INTEGER(ls);
  common = INTEGER(nc);
  jump = INTEGER(lt);
  ia = INTEGER(a);
  ib = INTEGER(b);
  for (i = 0; i < final; i++){
    p[0] = &dM[ia[i]];
    q[0] = &dV[ib[i]];
    for (l = 1; l < common[i]; l++) {
      p[l] = p[l - 1] + jump[i];
      q[l] = q[l - 1] + nrow;
    }
    k = 0;
    for (j = start[i]; j < end[i]; j++, k = (++k == nrow) ? 0 : k) {
      dX[j] = *(p[0])++ * q[0][k];
      for (l = 1; l < common[i]; l++){
	dX[j] += *(p[l])++ * q[l][k];
      }
    }
  }
  return(X);
}
/* computes single column of design matrix */
SEXP onecol(SEXP M, SEXP V, SEXP a, SEXP lt, SEXP nr, SEXP nc) {
  int j, k, l, nrow = INTEGER(nr)[0], common = INTEGER(nc)[0], jump;
  double *p[common], *q[common], *dcol;
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
  dcol = REAL(col);
  for (j = 0; j < nrow; j++, k = (++k == nrow) ? 0 : k) {
    dcol[j] = *(p[0])++ * q[0][k];
    for (l = 1; l < common; l++){
      dcol[j] += *(p[l])++ * q[l][k];
    }
  }
  UNPROTECT(1);
  return(col);
}

/* register routines */

static const R_CallMethodDef callMethods[] = {
    {"submatprod", (DL_FUNC) &submatprod, 5},
    {"subprod", (DL_FUNC) &subprod, 6},
    {"newsubprod", (DL_FUNC) &newsubprod, 12},
    {"onecol", (DL_FUNC) &onecol, 6},
    {NULL}
};

void R_init_gnm(DllInfo *info)
{
    /* Register the C and .Call routines.
     No .C(), .Fortran() or .External() routines,
     so pass those arrays as NULL.
     */
    R_registerRoutines(info,
                       NULL, callMethods,
                       NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
    R_forceSymbols(info, TRUE);
}
