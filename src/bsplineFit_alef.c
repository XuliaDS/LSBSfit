/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Skin a series of BSpline Curves
 *
 *      Copyright 2011-2017, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#ifdef EGADS_DEBUG
#define _GNU_SOURCE
#include <fenv.h>
#endif

#include "egads.h"
#include "egadsSplineFit.h"
#include <time.h>
#include <string.h>
#include <stdio.h>
#include <math.h>


#define PI  3.1415926535897931159979635
#define MIN(A,B)        (((A) < (B)) ? (A) : (B))
#define MAX(A,B)        (((A) < (B)) ? (B) : (A))
#define DEPS 1.e-12

#define MAX_FALSE_ITERS 10
#define MIN_PTS_PER_SPAN 2
#define MIN_KNOT_DIST 1.e-5
#define MAX_KNOTS 1000
#define EPS2 1.e-10
#define EPS1 1.e-10
#define NEWTON_IT_MAX 50

#define SAVE_OPTI_AND_CP 1
#define SAVE_ITER_NORMS 1
#define PRINT_KNOTS
#define L2_INIT 100.0
#define DEBUG 1
static double
L2norm(double f[],                      /* (in)  vector */
    int    n)                        /* (in)  length of vector */
{
  double L2norm;                        /* (out) L2-norm */
  int    i;

  L2norm = 0.0;
  for (i = 0; i < n; i++) L2norm += f[i] * f[i];

  L2norm = sqrt(L2norm);

  return L2norm;
}

static int
matsol(double    A[],           /* (in)  matrix to be solved (stored rowwise) */
    /* (out) upper-triangular form of matrix */
    double    b[],           /* (in)  right hand side */
    /* (out) right-hand side after swapping */
    int       n,             /* (in)  size of matrix */
    double    x[])           /* (out) solution of A*x=b */
{

  int    ir, jc, kc, imax;
  double amax, swap, fact;
  double EPS12=1.E-12;

  /* --------------------------------------------------------------- */

  /* reduce each column of A */
  for (kc = 0; kc < n; kc++) {

    /* find pivot element */
    imax = kc;
    amax = fabs(A[kc*n+kc]);

    for (ir = kc+1; ir < n; ir++) {
      if (fabs(A[ir*n+kc]) > amax) {
        imax = ir;
        amax = fabs(A[ir*n+kc]);
      }
    }

    /* check for possibly-singular matrix (ie, near-zero pivot) */
    if (amax < EPS12) return EGADS_DEGEN;

    /* if diagonal is not pivot, swap rows in A and b */
    if (imax != kc) {
      for (jc = 0; jc < n; jc++) {
        swap         = A[kc  *n+jc];
        A[kc  *n+jc] = A[imax*n+jc];
        A[imax*n+jc] = swap;
      }

      swap    = b[kc  ];
      b[kc  ] = b[imax];
      b[imax] = swap;
    }

    /* row-reduce part of matrix to the bottom of and right of [kc,kc] */
    for (ir = kc+1; ir < n; ir++) {
      fact = A[ir*n+kc] / A[kc*n+kc];

      for (jc = kc+1; jc < n; jc++) A[ir*n+jc] -= fact * A[kc*n+jc];

      b[ir] -= fact * b[kc];

      A[ir*n+kc] = 0;
    }
  }

  /* back-substitution pass */
  x[n-1] = b[n-1] / A[(n-1)*n+(n-1)];
  for (jc = n-2; jc >= 0; jc--) {
    x[jc] = b[jc];
    for (kc = jc+1; kc < n; kc++) x[jc] -= A[jc*n+kc] * x[kc];
    x[jc] /= A[jc*n+jc];
  }
  return EGADS_SUCCESS;
}


static int
EG_distance_projection(int *ivec, double *data, double *xyz, double *t, double *rms)
{
  int deg, nCP, stat, it = 0, j;
  double t0, eps1, eps2, p[3], pp[3], spline_vals[9], dist_vec[3], arc_der, f, df, aux, t1, error;
  deg  = ivec[1];
  nCP  = ivec[2];
  t1 = *t;
  t0 = t1;
  while(it < NEWTON_IT_MAX) {
    // spline_vals contains spline point (positions 0,1,2)
    // 1st derivative (pos 3,4,5) and 2nd derivative (pos 6,7,8)
    stat = EG_spline1dDeriv(ivec, data, 2, t1, spline_vals);
    if ( stat != EGADS_SUCCESS) {
      printf(" EG_spline1dEval = %d\n", stat);
      return stat;
    }
    eps1    = 0.0;
    arc_der = 0.0;
    f		= 0.0;
    df      = 0.0;
    for ( j = 0; j < 3; ++j) {
      dist_vec[j] = (spline_vals[j]-xyz[j]);///d_euclid;
      eps1       += pow(dist_vec[j],2);
      arc_der    += pow(spline_vals[3+j],2);
      f		   += spline_vals[3+j]*dist_vec[j];
      df         += spline_vals[6+j]*dist_vec[j];
    }
    eps1    = sqrt(eps1);
    df     += arc_der;
    arc_der = sqrt(arc_der);

    if (eps1 <= EPS1 || arc_der < DEPS || fabs(df) < DEPS) {
      error = eps1;
      break;
    }

    eps2    = fabs(f)/(arc_der*eps1);
    if (eps2 <= EPS2) {
      error = eps1;
      break;
    }

    // Newton iteration step
    t0 = t1;
    t1 = t0 - f/df;
    if ( t1 < data[0] || t1 > data[ivec[3] -1])	t1 = t0;
    if (fabs((t1-t0)*arc_der) <= EPS1  || it == (NEWTON_IT_MAX -1) ) {
      stat  = EG_spline1dEval(ivec, data, t1, p);
      stat += EG_spline1dEval(ivec, data, *t, pp);
      if ( stat != EGADS_SUCCESS) {
        printf(" EG_spline1dEval = %d\n", stat);
        return stat;
      }
      error = 0.0;
      aux   = 0.0;
      for ( j = 0 ; j < 3; ++j)	{
        error += pow(xyz[j]-p[j],2);
        aux   += pow(xyz[j]-pp[j],2);
      }
      error = sqrt(error);
      aux   = sqrt(aux);
      if ( aux < error)	error = aux;
      break;
    }
    ++it;
  }
  *t = t1;
  *rms = error;
  return EGADS_SUCCESS;
}
/* Computes Bsplines values for a given knot sequence */
static double
OneBasisFun (int p, int m, double U[], int i, double u)
{
  int    j, k;
  double saved, Uleft, Uright, temp, *N;
  // Special cases
  if ( ((i == 0) && (u == U[0])) ||  ((i == m-p-1) && (u == U[m])) ) return 1.0;
  if ( (u < U[i]) || (u >= U[i+p+1]) ) {
    return 0.0;
  }
  N = (double *) EG_alloc((p+1) *sizeof(double));
  if (N == NULL) return EGADS_MALLOC;
  for (j = 0; j <= p; ++j) {
    if ( ( u >= U[i+j] ) && ( u < U[i+j+1] ) ) {
      N[j] = 1.0;
    } else {
      N[j] = 0.0;
    }
  }
  for ( k = 1; k <= p; ++k ) {
    if ( N[0] == 0.0 ) {
      saved = 0.0;
    } else {
      saved = ( (u-U[i]) * N[0] )/( U[i+k]-U[i] );
    }
    for ( j = 0; j < p-k+1; ++j ) {
      Uleft  = U[i+j+1];
      Uright = U[i+j+k+1];
      if (N[j+1] == 0.0) {
        N[j]  = saved;
        saved = 0.0;
      } else {
        temp  = N[j+1]/( Uright - Uleft );
        N[j]  = saved + ( Uright - u )*temp;
        saved = ( u - Uleft )*temp;
      }
    }
  }
  saved = N[0];
  EG_free(N);
  return saved;
}

static int
getParameterVector(double *data, int dimData,double *t, int paramType)
{
  int k = 0, j = 0;
  double arc = 0, diam = 0.0;
  if ( t == NULL ) return EGADS_EMPTY;
  switch (paramType)
  {
  default:
  {   //centripetal knots
    for (k = 1; k < dimData; k++) {
      arc = 0.0;
      for(j = 0; j < 3; ++j)
        arc += pow( (data[3*k+j]-data[3*(k-1) + j]),2);
      diam += pow(arc,0.25);
    }
    t[0] = 0.0; t[dimData-1] = 1.0;
    for (k = 1; k < dimData - 1; k++) {
      arc  = pow( (data[3*k]   - data[3*(k-1)]), 2);
      arc += pow( (data[3*k+1] - data[3*(k-1)+1]), 2);
      arc += pow( (data[3*k+2] - data[3*(k-1)+2]), 2);
      t[k] = t[k-1] + pow(arc,0.25)/diam;
    }
    break;
  }
  case 1:    //chord-length knots
  {
    for (k = 1; k < dimData; k++) {
      arc = 0.0;
      for(j = 0; j < 3; ++j)
        arc += pow( (data[3*k+j]-data[3*(k-1) + j]),2);
      diam += pow(arc,0.5);
    }
    t[0] = 0.0; t[dimData-1] = 1.0;
    for (k = 1; k < dimData - 1; k++) {
      arc  = pow( (data[3*k]   - data[3*(k-1)]),2);
      arc += pow( (data[3*k+1] - data[3*(k-1)+1]),2);
      arc += pow( (data[3*k+2] - data[3*(k-1)+2]),2);
      t[k] = t[k-1] + pow(arc,0.5)/diam;
    }
    break;
  }
  }
  return EGADS_SUCCESS;
}

static int
getKnotVector(double *t,int dimT, double *u, int nKnots, int nCP, int deg)
{
  double alpha, d;
  int  j, i ;
  d = (double)(dimT)/(double)(nCP - deg);
  if ( t == NULL || u == NULL ) return EGADS_EMPTY;
  for ( j = 0; j <= deg; ++j)
  {
    u[j]                   = t[0];
    u[nKnots -1 - deg + j] = t[dimT-1];
  }
  for ( j = 1; j < nCP - deg; ++j)
  {
    i        = floor((double)j*d);
    alpha    = (double)j*d - (double)i ;
    u[j+deg] = (1.0 - alpha)*t[i-1] + alpha*t[i];
  }
  return EGADS_SUCCESS;
}


static int
EG_Linear_LS_Bspline_Optimization(int m, double XYZ[], double t[], int header[],
    double *spline, double *error, double *rms)
{
  double **vecB = NULL, **matA = NULL, **matTransA = NULL, **vecR = NULL, *r = NULL, *cp = NULL, tCopy;
  double basis_0, basis_n;
  int i, j, k, stat, it, d, nKnots, n, deg;
  deg    = header[1];
  nKnots = header[3];
  n      = header[2];
  r		  = (double*) EG_alloc (3*(m-2)*sizeof(double));
  vecB      = (double**)EG_alloc (3      *sizeof(double*));
  vecR      = (double**)EG_alloc (3      *sizeof(double*));
  matA      = (double**)EG_alloc ((m-2)  *sizeof(double*));
  matTransA = (double**)EG_alloc ((n-2)  *sizeof(double*));
  if ( matA  == NULL || vecB == NULL || vecR == NULL || matTransA == NULL || r == NULL) {
    return  EGADS_MALLOC;
  }
  for ( i = 0; i < n-2; ++i) {
    matTransA[i] = (double*)EG_alloc ((m-2)*sizeof(double));
    if (matTransA[i] == NULL ) {
      EG_free(matTransA);
      EG_free(matA);
      EG_free(vecB);
      EG_free(vecR);
      EG_free(r);
      return EGADS_MALLOC;
    }
  }
  for ( i = 0; i < m-2; ++i) {
    matA[i] = (double*)EG_alloc ((n-2)*sizeof(double));
    if (matA[i] == NULL ) {
      for ( j = 0; j < n-2; ++j)	EG_free(matTransA[j]);
      EG_free(matTransA);
      EG_free(matA);
      EG_free(vecB);
      EG_free(vecR);
      EG_free(r);
      return EGADS_MALLOC;
    }
  }
  for ( i = 0; i < 3; ++i) {
    vecB[i] = (double*)EG_alloc ((n-2)*(n-2)*sizeof(double));
    vecR[i] = (double*)EG_alloc ((n-2)*sizeof(double));
    if ( vecB[i] == NULL || vecR[i] == NULL ) {
      for ( j = 0; j < m-2; ++j)	EG_free(matA[j]);
      for ( j = 0; j < n-2; ++j)	EG_free(matTransA[j]);
      EG_free(matTransA);
      EG_free(matA);
      EG_free(r);
      EG_free(vecB);
      EG_free(vecR);
      return EGADS_MALLOC;
    }
  }
  for ( k = 0; k < m-2; ++ k ) {
    basis_0 = OneBasisFun(deg, nKnots-1, spline, 0,   t[k+1]);
    basis_n = OneBasisFun(deg, nKnots-1, spline, n-1, t[k+1]);
    for (d = 0 ; d < 3; ++d) {
      r[3*k+d] = XYZ[3*(k+1)+d] - basis_0*XYZ[d] - basis_n*XYZ[3*(m-1)+d];
    }
  }
  for ( i = 0; i < m-2; ++i) {
    for ( j = 0 ; j < n-2; ++ j)	{
      matA[i][j]      = OneBasisFun(deg, nKnots-1, spline, j+1, t[i+1]) ;
      matTransA[j][i] = matA[i][j];
    }
  }
  // Multiply matrices
  for ( it = j = 0 ; j < n-2; ++ j) {
    for ( i = 0 ; i < n-2; ++ i, ++it) {
      vecB[0][it] = 0.0;
      for (k = 0; k < m-2; ++k) vecB[0][it] += matTransA[j][k]*matA[k][i]; // solves for x
      vecB[1][it] = vecB[0][it]; // solves for y
      vecB[2][it] = vecB[0][it]; // solves for z
    }
  }
  for( i = 0; i < (n-2); ++i ) {
    vecR[0][i] = 0.0;
    vecR[1][i] = 0.0;
    vecR[2][i] = 0.0;
    for (k = 0; k < m-2; ++k) {
      vecR[0][i] += matTransA[i][k]*r[3*k];   // x-variable
      vecR[1][i] += matTransA[i][k]*r[3*k+1]; // y-variable
      vecR[2][i] += matTransA[i][k]*r[3*k+2]; // z-variable
    }
  }
  EG_free(r);
  for ( i = 0 ; i < m-2; ++ i) EG_free(matA[i]);
  for ( i = 0 ; i < n-2; ++ i) EG_free(matTransA[i]);
  EG_free(matTransA);
  EG_free(matA);
  cp  = (double*) EG_alloc ((n-2)*sizeof(double));
  if (cp == NULL) {
    stat = EGADS_MALLOC;
    for ( j = 0 ; j < 3; ++j) {
      EG_free(vecB[j]);
      EG_free(vecR[j]);
    }
    EG_free(vecB);
    EG_free(vecR);
    return stat ;
  }
  for ( i = 0 ; i < 3; ++i) {
    stat = matsol(vecB[i], vecR[i], (n-2), cp);
    if (stat != EGADS_SUCCESS) {
      for ( j = 0 ; j < 3; ++j) {
        EG_free(vecB[j]);
        EG_free(vecR[j]);
      }
      EG_free(vecB);
      EG_free(vecR);
      EG_free(cp);
      return stat ;
    }
    EG_free(vecB[i]);
    EG_free(vecR[i]);
    spline[nKnots]             = XYZ[0];
    spline[nKnots+1]           = XYZ[1];
    spline[nKnots+2]           = XYZ[2];
    spline[3*(n-1) + nKnots]   = XYZ[3*(m-1)];
    spline[3*(n-1) + nKnots+1] = XYZ[3*(m-1)+1];
    spline[3*(n-1) + nKnots+2] = XYZ[3*(m-1)+2];
    for ( j = 1; j < n-1 ; j++) {
      spline[3*j+nKnots+i] = cp[j-1];
    }
  }
  EG_free(cp);
  EG_free(vecB);
  EG_free(vecR);
  // Get curve to data error
  for (j = 0; j < m ; ++j ) {
    tCopy = t[j];
    stat = EG_distance_projection(header, spline, &XYZ[j*3], &tCopy, &error[j]);
    if (stat != EGADS_SUCCESS) {
      printf("EG_distance_projection = %d\n!!",stat);
      return stat;
    }
  }
  *rms = L2norm(error, m) / (double)m;
  return EGADS_SUCCESS;
}

/*
static int
checkSpan(double u0,double u1,double u2, double *t, int m, int *n)
{
	int i, spanOK = 0, u1u2Empty = 0;
	double alpha = 0.9;
	if (fabs ( u1-u2) < DEPS )	u1u2Empty = 1;
	for ( i = *n; i < m; ++i)
	{
		if (spanOK == 0 ) {
			if( t[i] >= u0 ) {
				if( t[i] > u1 ) // knotSpan is empty
					t[i] = (1.0 - alpha)*u0 + alpha*u1;
				spanOK = u1u2Empty + 1;
			}
		} else {
			if ( t[i] >= u1 ) {
				if ( t[i] > u2 ) // knotSpan is empty
					t[i] = (1.0 - alpha)*u1 + alpha*u2;
				++spanOK;
			}
		}
		if ( spanOK == 2) {
 *n = i + 1;
			return EGADS_SUCCESS;
		}
	}
 *n = i+1;
	return EGADS_EMPTY;
}
 */


static int
place_knot(int minIdx, int maxIdx, int centre, double left, double right, double *tVals, double *res) {
  int k = 0,  *iter = NULL, plus, minus;
  int cands;
  if ( (maxIdx - minIdx) == 0) {
    cands = 1;
    iter  = (int*)EG_alloc(cands*sizeof(int));
    if ( iter == NULL ) return EGADS_MALLOC;
    iter[0] = minIdx;
  } else {
    cands = maxIdx - minIdx;
    iter  = (int*)EG_alloc(cands*sizeof(int));
    if ( iter == NULL ) return EGADS_MALLOC;
    plus  = 1;
    minus = 1;
    k = 0;
    if ( centre < maxIdx) iter[k++] = centre;
    while ( k < cands) {
      if ( (centre - minus) >= minIdx) {
        iter[k] = centre - minus;
        ++minus;
        ++k;
      }
      if (centre + plus < maxIdx) {
        iter[k] = centre + plus;
        ++plus;
      }
      ++k;
    }
  }
  for ( k = 0; k < cands; ++k) {
    *res = 0.5*( tVals[iter[k]] + tVals[iter[k]+1]);
    if ( (fabs(*res - left) >= MIN_KNOT_DIST ) && (fabs(*res - right) >= MIN_KNOT_DIST ) ){//&& ( (iter[k] - minIdx) >= MIN_PTS_PER_SPAN - 1) && ( (maxIdx - (iter[k]+1)) >= MIN_PTS_PER_SPAN -1) ){
      EG_free(iter);
      return EGADS_SUCCESS;
    }
  }
  EG_free(iter);
  return EGADS_EMPTY;
}


//#define KNOT_BISECTION


static int
get_max_varirance_knot(int offset, // = degree (start looking only at the interior knots)
    int nSpans,            // number of iterior knots
    int *badKnotsVec,     // excludes previously refined knots that gave bad results
    int *totBadKnots,  // number of bad knots
    double bSpline[], // knot vector
    double error[],   // error at tvals
    double *t,        // curve evaluation points
    int m,            // number of data points
    int *errId,       // stores the span with maximum error
    double *tVal)     // stores the t with maximum error at the maximum span
{
  double  *knotErr  = NULL, *weightedErr, errMax = 0.0, sum, tAux;
  int     *knotPool = NULL, *ids = NULL;
  int     i, j, checkKnot,  first_t, stat;
  knotPool     = (int*)   EG_alloc(nSpans*sizeof(int));
  ids          = (int*)   EG_alloc(nSpans*sizeof(int));
  knotErr      = (double*)EG_alloc(nSpans*sizeof(double));
  weightedErr  = (double*)EG_alloc(nSpans*sizeof(double));
  if (knotPool == NULL || knotErr == NULL || weightedErr == NULL || ids == NULL ) return EGADS_MALLOC;
  for ( i = 0; i < nSpans; ++ i) {
    knotPool[i]  =   0;
    ids[i]       =   0;
    checkKnot    =   1;
    first_t      =  -1;
    knotErr[i]   = 0.0;
    weightedErr[i]   = 0.0;
    errMax       = 0.0;
    for ( j = 0; j < *totBadKnots; ++j) { // check if we already have used that knot
      if ( (i+offset) == badKnotsVec[j])  {
        checkKnot = 0;
        break;
      }
    }
    if ( checkKnot == 1 ) {
      for (j = 1; j < m - 1 ; ++j )	{
        if ( ( t[j] >= bSpline[i+offset] ) && (t[j] < bSpline[i+1+offset] ) ) {
          knotPool[i]++;
          knotErr[i] += error[j]*error[j];
          if (first_t == -1) {
            ids[i]  = j;
            first_t = 0;
          }
        }
      }
      if (knotPool[i] >=  2*MIN_PTS_PER_SPAN )	{  // min points per span.
        weightedErr[i] = knotErr[i]/(double)knotPool[i];
      }
      else {
        weightedErr[i] = -1.0;
        badKnotsVec[*totBadKnots] = i + offset;
        ++(*totBadKnots);
      }
    }
  }
  *errId = -1;
  errMax = 0.0;
  for ( i = 0; i < nSpans; ++ i) {
    if ( (weightedErr[i] > errMax)  && (bSpline[i + offset + 1] - bSpline[i+offset] >= 2.*MIN_KNOT_DIST) ) {
      errMax = weightedErr[i];
      *errId = i;
    }
  }
#ifdef KNOT_BISECTION
  if (*errId != -1) {
    *tVal = 0.5*(bSpline[*errId +offset] + bSpline[*errId +offset + 1] );
    if      (*tVal < t[ids[*errId]] )                    *tVal = 0.5*(t[ids[*errId]] + t[ids[*errId]+1]);
    else if (*tVal > t[ids[*errId]+knotPool[*errId]-1] ) *tVal = 0.5*( t[ids[*errId]+knotPool[*errId]-2] + t[ids[*errId]+knotPool[*errId]-1]);
    stat  = EGADS_SUCCESS;
  }
#else
  if (*errId != -1) {
    if ( knotPool[*errId] == 2*MIN_PTS_PER_SPAN) {
      *tVal = 0.5*(t[ids[*errId]] +  t[ids[*errId] +1]);
      stat  = EGADS_SUCCESS;
    } else {
      sum = 0.0;
      for ( j = 0 ; j < MIN_PTS_PER_SPAN ; ++j) sum  += error[ids[*errId]+j]*error[ids[*errId]+j];
      for ( j = MIN_PTS_PER_SPAN -1 ; j <= knotPool[*errId] - MIN_PTS_PER_SPAN ; ++j) {
        if ( (sum >= knotErr[*errId]*0.5) || (j == knotPool[*errId] - MIN_PTS_PER_SPAN)) {
          stat = place_knot(ids[*errId] + MIN_PTS_PER_SPAN, ids[*errId] + knotPool[*errId] - MIN_PTS_PER_SPAN, ids[*errId] + j, bSpline[*errId+offset], bSpline[*errId+offset+1], t, &tAux);
          if ( stat == EGADS_SUCCESS) {
            j     = m;
            *tVal = tAux;
          }
        }
        else sum += error[ids[*errId]+j]*error[ids[*errId]+j];

      }
    }
  }
#endif
  EG_free(knotErr);
  EG_free(weightedErr);
  EG_free(knotPool);
  EG_free(ids);
  if (*errId == -1)    return EGADS_EMPTY;
  return EGADS_SUCCESS;
}



static intcreateKnotSpan( int header[], double *bSpline, int m, double t[], double XYZ[], double *error, int totBadSpans, int badSpanVec[],
    int *spanLocation,  double  *newKnot, double *inL2error, double updateTOL)
{
  double *errAtIt = NULL, *bestError = NULL, *localSpline = NULL, *bestSpline = NULL;
  double  prevL2error, L2error = 0.0,  bestL2error , bestNewKnot, tVal;
  int     nKnots, nCP, deg, stat, nSpans, locBadSpans;
  int     i, j, bestSpanLocation,  tPos;
  prevL2error = *inL2error;
  nKnots      = header[3];
  nCP         = header[2];
  deg         = header[1];
  nSpans      = nKnots - 2*(deg+1);
  errAtIt     = (double*)EG_alloc(m*sizeof(double));
  bestError   = (double*)EG_alloc(m*sizeof(double));
  localSpline = (double*)EG_alloc((nKnots+3*nCP)*sizeof(double));
  bestSpline  = (double*)EG_alloc((nKnots+3*nCP)*sizeof(double));
  if (errAtIt == NULL || bestError == NULL || localSpline == NULL || bestSpline == NULL) return EGADS_MALLOC;
  locBadSpans      = totBadSpans;
  bestL2error      = L2_INIT;
  bestSpanLocation = -1;
  do {
    stat = get_max_varirance_knot(deg, nSpans, badSpanVec, &locBadSpans, bSpline, error, t,  m, &tPos, &tVal);
    if (stat == EGADS_SUCCESS) {
      for ( j = 0 ; j <= deg; ++j)
      {  // clamped: copy first & last
        localSpline[j]             = bSpline[j];
        localSpline[nKnots -1 - j] = bSpline[ (nKnots -1) -1 -j];
      }
      for ( j = i = 0 ; i < nSpans; ++i, ++j) {
        if ( i != tPos ) localSpline[j + deg + 1] = bSpline[i + deg + 1];
        else {
          localSpline[j+ deg + 1]  = tVal;
          ++j;
          localSpline[j + deg + 1] = bSpline[i + deg + 1];
        }
      }
      stat = EG_Linear_LS_Bspline_Optimization(m, XYZ, t, header, localSpline,  errAtIt, &L2error);
      if (stat != EGADS_SUCCESS) {
        printf("Solution To Linear System = %d!!!\n",stat);
        goto cleanup;
      }
      if ( L2error < prevL2error - updateTOL) {
        *spanLocation = tPos +deg;
        *newKnot      = tVal;
        *inL2error    = L2error;
        for ( i = 0 ; i < m; ++i) error[i] = errAtIt[i];
        for ( i = 0; i < header[3] + 3*header[2]; ++i ) bSpline[i] = localSpline[i];
        goto cleanup;
      }
      if ( L2error < bestL2error ) {
        for ( i = 0; i < header[3] + 3*header[2]; ++i ) bestSpline[i] = localSpline[i];
        for ( i = 0 ; i < m; ++i) bestError[i] = errAtIt[i];
        bestL2error      = L2error;
        bestSpanLocation = tPos +deg;
        bestNewKnot      = tVal;
      }
    }
    badSpanVec[locBadSpans] = tPos +deg;
    ++locBadSpans;
  }
  while ( (locBadSpans < nSpans) && (stat == EGADS_SUCCESS));
  if (bestSpanLocation != -1 ){
    for ( i = 0; i < header[3] + 3*header[2]; ++i ) bSpline[i] = bestSpline[i];
    for ( i = 0 ; i < m; ++i) error[i] = bestError[i];
    stat                     = EGADS_SUCCESS;
    *spanLocation            = bestSpanLocation;
    *newKnot                 = bestNewKnot;
    *inL2error               = bestL2error;
  }
  cleanup:
  EG_free(localSpline);
  EG_free(bestSpline);
  EG_free(errAtIt);
  EG_free(bestError);

  return stat;
}


static int
get_minimal_knots (int header[], double spline[], int m, double *error, double t[], double XYZ[], int totNewKnots, double knotSeq[], int *outHeader, double **outSpline, double *oriL2, int *reductions)
{
  int        stat, i, j, nKnots, deg,  nSpans, pos, it ;
  double     newL2error,  *knotCopy = NULL, *backKnot = NULL, *newSpline = NULL;
  knotCopy   = (double*)EG_alloc((header[3] + totNewKnots )*sizeof(double));
  backKnot   = (double*)EG_alloc((header[3] + totNewKnots )*sizeof(double));
  newSpline  = (double*)EG_alloc( (3*header[2]+header[3])*sizeof(double));
  if (backKnot == NULL || knotCopy == NULL || newSpline == NULL) return EGADS_MALLOC;
  nKnots        = header[3] + 1;
  deg           = header[1];
  for ( i = 0; i < nKnots; ++i ) knotCopy[i] = spline[i];
  it  = 1;
  pos = totNewKnots - it;
  do
  {
    outHeader [3] = nKnots;
    outHeader [2] = nKnots - deg -1;
    newSpline     = (double*)EG_reall(newSpline, (3*outHeader[2] + outHeader[3])*sizeof(double));
    for ( j = 0 ; j <= deg; ++j)
    {  // clamped: copy first & last
      newSpline[j]              = knotCopy[j];
      newSpline[nKnots -1 - j]  = knotCopy[ (nKnots-1) - 1 - j];
    }
    nSpans = outHeader[3] - 2*(deg+1);
    j = 0;
    for (  i = 0 ; i < nSpans; ++i ) {
      if ( (knotSeq[pos] >= knotCopy[i+deg]) && (knotSeq[pos] < knotCopy[i+deg+1]) ) {
        newSpline[deg+j+1] = knotSeq[pos];
        ++j;
      }
      newSpline[deg+j+1] = knotCopy[deg+i+1];
      ++j;
    }
    if ( pos == totNewKnots - it ) {
      for ( i = 0 ; i < nKnots; ++i) {
        backKnot[i] = newSpline[i];
      }
    }
    stat = EG_Linear_LS_Bspline_Optimization(m, XYZ, t, outHeader, newSpline, error, &newL2error);
    if (stat != EGADS_SUCCESS ) {
      printf("EG_Linear_LS_Bspline_Optimization %d\n",stat);
      EG_free(knotCopy);
      EG_free(backKnot);
      EG_free(newSpline);
      return stat;
    }
    if ( newL2error <= *oriL2 ) {
      *outSpline = (double*)EG_reall(*outSpline, (3*outHeader[2] + outHeader[3])*sizeof(double));
      if (*outSpline == NULL ) {
        EG_free(knotCopy);
        EG_free(backKnot);
        EG_free(newSpline);
        return EGADS_MALLOC;
      }
      for ( i = 0; i < 3*outHeader[2] + outHeader[3]; ++i ) {
        (*outSpline)[i] = newSpline[i];
      }
      *oriL2 = newL2error;
      EG_free(newSpline);
      EG_free(backKnot);
      EG_free(knotCopy);
      if ( outHeader[3] != totNewKnots + header[3]) (*reductions)++;
      return EGADS_SUCCESS;
    }
    --pos;
    if ( pos == -1 ) {
      for ( i = 0 ; i < nKnots; ++ i ) {
        knotCopy[i] = backKnot[i];
      }
      ++nKnots;
      ++it;
      pos = totNewKnots - it;
    }
  }
  while( it <= totNewKnots);
  outHeader[3] = outHeader[3] + totNewKnots;
  outHeader[2] = outHeader[3] - deg - 1;
  *outSpline = (double*)EG_reall(*outSpline, (3*outHeader[2] + outHeader[3])*sizeof(double));
  if (*outSpline == NULL ) {
    EG_free(knotCopy);
    EG_free(newSpline);
    EG_free(backKnot);
    return EGADS_MALLOC;
  }
  for ( i =0; i < 3*outHeader[2] + outHeader[3]; ++i ) (*outSpline)[i] = newSpline[i];
  EG_free(knotCopy);
  EG_free(newSpline);
  EG_free(backKnot);
  return EGADS_SUCCESS;
}



int
EG_1dBsplineCurveFit(ego context, double *XYZcloud, int m, ego *curve, int n,
    int deg, double *rms, double LS_tol, int stepControl, int maxIter, double delta)
{
  double updateTOL, L2error, prevL2error, backL2error, bestL2error;
  double  *t = NULL,  *bSpline = NULL,  *bSplineCopy = NULL, *backSpline = NULL, bestKnots[MAX_KNOTS], *XYZ = NULL, *error= NULL, redundantKnotSeq[MAX_KNOTS];
  double point[3], scale, xmin, xmax, ymin, ymax, zmin, zmax, xcent, ycent, zcent,  newKnot;
  int    falseIter, headerCopy[4], header[4], backHeader[4], badKnotSeq[MAX_KNOTS];
  int    reductions, accepted, spanLocation,  badCount, redundantKnots, update_step, nBestKnots;
  int    iter, stat, i, j, k, nKnots;
  int    paramType  = 0; // 1 = chord length, otherwise it uses centripetal
  FILE   *fn, *fbSpline;
  ego    bSplineCurve;
  nKnots    = n + deg + 1;
  header[0] = 0 ;
  header[1] = deg;
  header[2] = n;
  header[3] = nKnots;
  if (XYZcloud == NULL ) {
    printf(" No Data...\n");
    return EGADS_NODATA;
  }
  if ( n > m ) {
    printf(" There are more Control Points than data!! %d > %d\n. I'm going to try and fit directly\n",header[2],m);
    header[2] = m ;
    header[3] = header[2] + header[1] + 1;
    bSpline   = (double*)EG_alloc((header[2]*3 + header[3]) * sizeof(double));
    if ( bSpline == NULL ) return EGADS_MALLOC;
    for ( j = 0 ; j <= header[1]; ++j)
    {  // clamped: copy first & last
      bSpline[j]                = 0.0;
      bSpline[header[3] -1 - j] = 1.0;
    }
    nKnots = header[3] - 2*(header[1]+1) + 1;
    for ( i = 0; i < nKnots; ++ i ) {
      bSpline[i + header[1]+1 ] = (double)(i+1.0)/(double)nKnots;
    }
    for ( i = 0 ; i < header[2]; ++ i ) {
      bSpline[header[3] + 3*i ]   = XYZcloud[3*i];
      bSpline[header[3] + 3*i +1] = XYZcloud[3*i+1];
      bSpline[header[3] + 3*i +2] = XYZcloud[3*i+2];
    }
    stat = EG_makeGeometry(context, CURVE, BSPLINE, NULL, header, bSpline, &bSplineCurve);
    if (stat != EGADS_SUCCESS) {
      printf(" EG_makeGeometry %d\n",stat);
      EG_free(bSpline);
      return stat;
    }
    t     = (double*)EG_alloc (m*sizeof(double)); // Curve parameter
    error = (double*)EG_alloc (m*sizeof(double));
    if (error == NULL || t == NULL) {
      EG_free(bSpline);
      return EGADS_MALLOC;
    }
    stat  = getParameterVector(XYZcloud, m, t, paramType);
    for ( i = 0 ; i < m; ++ i ) {
      EG_distance_projection(header, bSpline, &XYZcloud[3*i], &t[i], &error[i]);
    }
    *rms = L2norm(error, m) / (double)m;
    *curve = bSplineCurve;
    printf("\n ---- LEAST SQUARES FIT STATUS %d ----\n"
        "- BEST FIT: %d CONTROL POINTS\n"
        "- RMS %11.4e\n"
        "---------------------------------------\n",stat,header[2],*rms);
    EG_free(t);
    EG_free(error);
    EG_free(bSpline);
    return  stat;
  }

  XYZ = (double *)EG_alloc(3*m*sizeof(double));
  if (XYZ == NULL ) return EGADS_MALLOC;
  // Scale Data: get max length and scale from centre
  xmin = XYZcloud[0];
  xmax = XYZcloud[0];
  ymin = XYZcloud[1];
  ymax = XYZcloud[1];
  zmin = XYZcloud[2];
  zmax = XYZcloud[2];
  for (k = 1; k < m; k++) {
    if (XYZcloud[3*k  ] < xmin) xmin = XYZcloud[3*k  ];
    if (XYZcloud[3*k  ] > xmax) xmax = XYZcloud[3*k  ];
    if (XYZcloud[3*k+1] < ymin) ymin = XYZcloud[3*k+1];
    if (XYZcloud[3*k+1] > ymax) ymax = XYZcloud[3*k+1];
    if (XYZcloud[3*k+2] < zmin) zmin = XYZcloud[3*k+2];
    if (XYZcloud[3*k+2] > zmax) zmax = XYZcloud[3*k+2];
  }
  fn =fopen("scaled_data","w");
  for (k = 0; k < m; k++) {
    if ( (xmax - xmin) > 0 )
      XYZ[3*k]    = (XYZcloud[3*k  ]  - xmin)/(xmax - xmin);
    if ( (ymax - ymin) > 0 )
      XYZ[3*k +1] = (XYZcloud[3*k +1] - ymin)/(ymax - ymin);
    if ( (zmax - zmin) > 0 )
      XYZ[3*k +2] = (XYZcloud[3*k +2] - zmin)/(zmax - zmin);
    fprintf(fn,"%.15lf\t%.15lf\t%.15lf\n",XYZ[3*k],XYZ[3*k+1],XYZ[3*k+2]);
  }
  fclose(fn);
  /*scale = 1.0 / MAX(MAX(xmax-xmin, ymax-ymin), zmax-zmin);
  if ( scale < DEPS ) {
    printf("This data is not suitable for fitting.\n");
    EG_free(XYZ);
    return EGADS_EMPTY;
  }
  xcent = 0.5* (xmin + xmax);
  ycent = 0.5* (ymin + ymax);
  zcent = 0.5* (zmin + zmax);
  fn =fopen("scaled_data","w"); 
  for (k = 0; k < m; k++) {
    XYZ[3*k  ] = scale * (XYZcloud[3*k  ] - xcent);
    XYZ[3*k+1] = scale * (XYZcloud[3*k+1] - ycent);
    XYZ[3*k+2] = scale * (XYZcloud[3*k+2] - zcent);
    fprintf(fn,"%.15lf\t%.15lf\t%.15lf\n",XYZ[3*k],XYZ[3*k+1],XYZ[3*k+2]);
  }
  fclose(fn);*/
  // ALLOCATE MEMO
  t             = (double*) EG_alloc (m*sizeof(double));	// Curve parameter
  bSpline       = (double*) EG_alloc ( (nKnots+ n*3)*sizeof(double));
  bSplineCopy   = (double*) EG_alloc ( (nKnots+ n*3)*sizeof(double));
  backSpline    = (double*) EG_alloc ( (nKnots+ n*3)*sizeof(double));
  error         = (double*) EG_alloc (m *sizeof(double));
  if ( t    == NULL || bSpline == NULL || error == NULL || bSplineCopy == NULL || backSpline == NULL) {
    stat = EGADS_MALLOC;
    goto cleanup;
  }
  // Get parametrization
  stat = getParameterVector(XYZ, m, t, paramType);
  // Get associated knot vector
  stat    = getKnotVector(t, m, bSpline, nKnots, n, deg );
  if ( stat != EGADS_SUCCESS) goto cleanup;
  // LEAST SQUARES + KNOT INSERTION ROUTINE
#if SAVE_ITER_NORMS
  fn   = fopen("iter_norms","w");
  if (fn == NULL ) {
    stat = EGADS_MALLOC;
    goto cleanup;
  }
#endif
  iter          = 0;
  reductions    = 0;
  accepted      = 0;
  headerCopy[0] = 0;
  headerCopy[1] = deg;
  stat = EG_Linear_LS_Bspline_Optimization(m, XYZ, t, header, bSpline, error, &L2error);
  if (stat != EGADS_SUCCESS) {
    printf("Solution To Linear System = %d!!!\n",stat);
    goto cleanup;
  }
  for ( i = 0 ; i < header[3] + 3*header[2]; ++i) backSpline[i] = bSpline[i];
  for ( i = 0; i < 4; ++i) backHeader[i] = header[i];
  prevL2error    = L2error;
  printf(" STARTING AT A ERROR %lf\n",L2error);
  updateTOL      = 0.0;
  falseIter      = 0;
  spanLocation   = -1;
  badCount       = 0;
  redundantKnots = 0;
  bestL2error    = L2_INIT;
  //LS_tol        *= scale;
  while ( (L2error > LS_tol) && (iter < maxIter) && (n < m) ){
    ++iter;
    update_step   = 1;
    headerCopy[3] = header[3] + 1;
    headerCopy[2] = headerCopy[3]-deg-1;
    bSplineCopy   = (double*)EG_reall(bSplineCopy, (headerCopy[3]+ headerCopy[2]*3)*sizeof(double));
    for ( i = 0; i < header[3]; ++i ) {
      bSplineCopy[i] = bSpline[i];
    }
    stat = createKnotSpan(headerCopy, bSplineCopy, m, t, XYZ, error, badCount, badKnotSeq, &spanLocation, &newKnot, &L2error, updateTOL);
    if ( stat != EGADS_SUCCESS ) {
      printf(" Cannot refine more the Bspline. Leave Fitting with 'best' approximation %d\n",stat);
      for ( i = 0 ; i < 4; ++i) header[i] = backHeader[i];
      bSpline = (double*)EG_reall(bSpline,(header[3]+header[2]*3)*sizeof(double));
      for ( i = 0; i < header[3] + header[2]*3; ++i ) bSpline[i] = backSpline[i];
      stat    = EGADS_SUCCESS;
      break;

    } else {
      for ( i = 0 ; i < 4; ++i) header[i] = headerCopy[i];
      bSpline = (double*)EG_reall(bSpline,(header[3]+header[2]*3)*sizeof(double));
      for ( i = 0; i < header[3] + header[2]*3; ++i ) bSpline[i] = bSplineCopy[i];
      redundantKnotSeq[redundantKnots] = newKnot;
      ++redundantKnots;
    }/*
    if ( (L2error > prevL2error)  || (fabs(prevL2error - L2error) < updateTOL) ) {
      printf("\n BAD KNOT INSERTION L2 %11.4e  > %11.4e  FALSE IT %d   RED KNOTS %d  BAD COUNT %d\n",L2error,prevL2error, falseIter, redundantKnots,badCount);
      update_step = 0 ;  // try to refine and improve the error
      if ( L2error < bestL2error ) {
        bestL2error = L2error;
        nBestKnots  = header[3];
        for ( i = 0 ; i < header[3]; ++i) bestKnots[i] = bSpline[i];
      }
      for ( i = deg ; i < header[3] - (deg+1); ++i) {
        if ( (i < spanLocation - header[1])  || (i > spanLocation + 2 ) ){
          for ( j = 0 ; j < badCount;  ++j) {
            if ( badKnotSeq[j] == i ) {
              j = MAX_KNOTS;
            }
          }
          if ( (j < MAX_KNOTS) ) {
            badKnotSeq[badCount] = i;
            ++badCount;
          }
        }
      }
      if (  (badCount >= header[3] - 2*(deg + 1))) {
        header[3] = nBestKnots;
        header[2] = nBestKnots - deg - 1;
        bSpline   = (double*)EG_reall(bSpline,(header[3]+header[2]*3)*sizeof(double));
        for ( i = 0; i < header[3]; ++i ) bSpline[i] = bestKnots[i];
        L2error        = bestL2error;
        update_step    = 1;
        redundantKnots = 0;
        falseIter      = 0;
      }
    }*/
    if ( (L2error > prevL2error)  || ( fabs(prevL2error - L2error) < updateTOL ) ){
      ++falseIter;
      update_step = 0;
      if ( (falseIter >= MAX_FALSE_ITERS) && (prevL2error > L2error) ) update_step = 1;
      else {
        for ( i = 0 ; i < header[3]; ++i) {
          if ( (i < spanLocation - header[1])  || (i > spanLocation + 2 ) ){
            badKnotSeq[badCount] = i;
            ++badCount;
          }
        }
        if ( badCount > header[3] -1) update_step = 1;
      }
    }
    if ( update_step == 1 ) {
      //printf("\n\n ========= CLEAN SEQUENCE ======  %d\t",header[3]);
      if ( falseIter > 0) {
        stat = get_minimal_knots(backHeader, backSpline, m, error, t, XYZ, redundantKnots, redundantKnotSeq, header, &bSpline, &L2error, &reductions);
        if (stat != EGADS_SUCCESS) {
          printf(" stat GET_MINIMAL_KNOTS %d!!\n",stat);
          goto cleanup;
        }
      }
      redundantKnots = 0 ;
      if ( fabs ( L2error - prevL2error) >= updateTOL) falseIter      = 0 ;
      ++accepted;
      updateTOL      = MAX(L2error*fabs(log10(prevL2error/L2error)),0.1*LS_tol) ;
#if SAVE_ITER_NORMS
      fprintf(fn, "%d\t%11.5e\t%11.5e\t%11.4e\n",iter,L2error, updateTOL,fabs(log10(L2error/prevL2error)));
#endif
      bestL2error    = L2_INIT;
      spanLocation   = -1;
      badCount       = 0 ;
      badKnotSeq[0] = -1;
      prevL2error    = L2error;
      backL2error    = L2error;
      for ( i = 0 ; i < 4; ++i)	backHeader[i] = header[i];
      backSpline  = (double*)EG_reall(backSpline,(header[3]+header[2]*3)*sizeof(double));
      for ( i = 0; i < header[3] + header[2]*3; ++i )	backSpline[i] = bSpline[i];
    }

  }
  printf(" Leaving Solver with norm %11.4e < %11.4e\n TOTALS: accepted %d rejected %d \n",L2error, LS_tol, accepted, reductions);
#if SAVE_ITER_NORMS
  fclose(fn);
#endif
  if (stat != EGADS_SUCCESS) goto cleanup;
  /*for (j = 0; j < header[2]; j++) {
    bSpline[3*j   + header[3]] = xcent + bSpline[3*j   + header[3]] / scale;
    bSpline[3*j+1 + header[3]] = ycent + bSpline[3*j+1 + header[3]] / scale;
    bSpline[3*j+2 + header[3]] = zcent + bSpline[3*j+2 + header[3]] / scale;
  }*/
  for (j = 0; j < header[2]; j++) {
    bSpline[3*j + header[3]]   = (1.0 - bSpline[3*j + header[3]] )  *xmin + bSpline[3*j + header[3]]  *xmax;
    bSpline[3*j + header[3]+1] = (1.0 - bSpline[3*j + header[3]+1] )*ymin + bSpline[3*j + header[3]+1]*ymax;
    bSpline[3*j + header[3]+2] = (1.0 - bSpline[3*j + header[3]+2] )*zmin + bSpline[3*j + header[3]+2]*zmax;
  }
  stat   = EG_makeGeometry(context, CURVE, BSPLINE, NULL, header, bSpline, &bSplineCurve);
  if (stat != EGADS_SUCCESS) {
    printf(" EG_makeGeometry %d\n",stat);
    goto cleanup;
  }
  *curve   = bSplineCurve;
  //*rms     = L2error/scale;
  *rms     = L2error;
#if SAVE_OPTI_AND_CP
  fbSpline = fopen("OPTIMAL_BSPLINE","w");
  if( fbSpline == NULL) {
    stat = EGADS_MALLOC;
    goto cleanup;
  }
  for ( i = 0; i < m; ++i) {
    stat = EG_evaluate(bSplineCurve, &t[i], point);
    if (stat != EGADS_SUCCESS) {
      printf(" EG_evaluate %d\n",stat);
      goto cleanup;
    }
    else
      fprintf(fbSpline,"%.11f\t%.11f\t%.11f\t%.11f\n",point[0],point[1],point[2],t[i]);
  }

  fbSpline = fopen("CONTROL_POLYGON","w");
  if( fbSpline == NULL) {
    stat = EGADS_MALLOC;
    goto cleanup;
  }
  for ( j = 0; j < header[2]; ++j) {
    fprintf(fbSpline,"%.11f\t%.11f\t%.11f\n", bSpline[3*j   + header[3]], bSpline[3*j   + header[3]+1], bSpline[3*j   + header[3]+2]);
  }
  fclose(fbSpline);
#ifdef PRINT_KNOTS
  fbSpline = fopen("KNOT_SEQUENCE","w");
  if( fbSpline == NULL) {
    stat = EGADS_MALLOC;
    goto cleanup;
  }
  double t0[3], t1[3];
  stat = EG_evaluate(bSplineCurve, &t[0], t0);
  stat = EG_evaluate(bSplineCurve, &t[m-1], t1);

  for ( j = 0; j < header[3]; ++j) {
    ycent  = (t1[0] - t0[0]);//*(t1[0] - t0[0]);
    //ycent += (t1[1] - t0[1])*(t1[1] - t0[1]);
    //ycent += (t1[2] - t0[2])*(t1[2] - t0[2]);
    // ycent =  sqrt(ycent);
    xcent =  bSpline[j];
    fprintf(fbSpline,"%.11f -1.0\n",xcent);
  }
  fclose(fbSpline);
#endif
#endif
  cleanup:
  EG_free(XYZ);
  EG_free(t);
  EG_free(error);
  EG_free(bSpline);
  EG_free(bSplineCopy );
  EG_free(backSpline);
  printf("\n ---- LEAST SQUARES FIT STATUS %d ----\n"
      "- BEST FIT: %d CONTROL POINTS\n"
      "- RMS %11.4e\n"
      "- EXIT AT ITERATION %d ? %d \n"
      "---------------------------------------\n",stat,header[2],L2error,iter,maxIter);
  return  stat;
}



//    ------------  B spline Splitting Functions -----------------   //

static int
FindSpan(int n, int p, double u, double *U)
{
  int low, mid, high;
  if ( fabs ( u - U[n+1] ) < DEPS) return n;
  low  = p; high = n+1;
  mid  = floor(( low + high )/2);
  while (  u < U[mid] -DEPS || ( u >= U[mid+1] - DEPS ) ) {
    if ( u < U[mid] ) {
      high = mid;
    } else {
      low  = mid;
    }
    mid =( low + high )/2;
    if ( fabs (u - U[mid] ) < DEPS) {
      while ( mid > 0  && ( fabs( U[mid] - U[mid -1] ) < DEPS ) ) --mid;
      return mid;
    }
  }
  if ( fabs (u - U[mid] ) < DEPS) {
    while ( mid > 0  && ( fabs( U[mid] - U[mid -1] ) < DEPS ) ) --mid;
  }
  return mid;
}

static int
getMultiplicity(double u, double *vector, int idx, int size)
{
  int j, m;
  m = 0;
  for ( j = idx; j < size; ++j ) {
    if ( fabs ( u - vector[j] ) < DEPS ) ++m;
  }
  return m;
}

static int
EG_CurveKnotIns(int *header,
    double **spline, //spline= knots + cps GETS OVERWRITTEN !
    double u,        // knot to be inserted
    int k,           // knot span where u lives
    int s,           // original multiplicity of u
    int r)           // # of insertions
{
  int    mp, nq, i, j, L, p, np, knotOffset;
  double *copySpline = NULL, *newSpline = NULL, *res = NULL, alpha;
  p  = header[1];
  np = header[2] - 1;
  mp = np + p +1;
  nq = np + r;      // # new control points
  copySpline = EG_alloc(((np+1)*3 + mp+1)  *sizeof(double));
  newSpline  = EG_alloc(((nq+1)*3 + mp+r+1)*sizeof(double));
  res        = EG_alloc(3*(p+1)            *sizeof(double));
  if (copySpline == NULL || res == NULL || newSpline == NULL ) return EGADS_MALLOC;
  for ( i = 0; i < ((np+1)*3 + mp+1); ++ i) copySpline[i] = (*spline)[i];
  // new knot sequence
  for (i = 0; i <= k; ++i)   newSpline[i]    = copySpline[i];
  for (i = 1; i <= r; ++i)   newSpline[k+i]  = u;
  for (i = k+1; i <= mp; ++i) newSpline[r+i] = copySpline[i];
  // Unaltered control points
  knotOffset = mp + 1 + r ;
  for ( i = 0 ; i <= k - p ; ++i) {  // unaltered CPs
    newSpline[knotOffset + i*3 ]   = copySpline[mp +1 + i*3 ];
    newSpline[knotOffset + i*3 +1] = copySpline[mp +1 + i*3 +1];
    newSpline[knotOffset + i*3 +2] = copySpline[mp +1 + i*3 +2];
  }

  for ( i = k-s ; i <= np ; ++i) {  // unaltered CPs
    newSpline[knotOffset + (i+r)*3 ]   = copySpline[mp+1 + i*3 ];
    newSpline[knotOffset + (i+r)*3 +1] = copySpline[mp+1 + i*3 +1];
    newSpline[knotOffset + (i+r)*3 +2] = copySpline[mp+1 + i*3 +2];
  }
  for ( i = 0; i <= p-s ; ++i) {
    res[3*i]   = copySpline[mp+1 + 3*(k-p+i)];
    res[3*i+1] = copySpline[mp+1 + 3*(k-p+i)+1];
    res[3*i+2] = copySpline[mp+1 + 3*(k-p+i)+2];
  }
  for ( j = 1; j <= r; ++j){
    L = k - p + j;
    for ( i = 0; i <= p-j-s; ++i){
      alpha    = (u - copySpline[L+i])/(copySpline[i+k+1]-copySpline[i+L]);
      res[3*i]   = alpha*res[3*(i+1)]   + (1.0-alpha)*res[3*i];
      res[3*i+1] = alpha*res[3*(i+1)+1] + (1.0-alpha)*res[3*i+1];
      res[3*i+2] = alpha*res[3*(i+1)+2] + (1.0-alpha)*res[3*i+2];
    }
    newSpline[knotOffset + 3*L]            = res[0];
    newSpline[knotOffset + 3*L+1]          = res[1];
    newSpline[knotOffset + 3*L+2]          = res[2];
    newSpline[knotOffset + 3*(k+r-j-s)]    = res[3*(p-j-s)];
    newSpline[knotOffset + 3*(k+r-j-s) +1] = res[3*(p-j-s) +1];
    newSpline[knotOffset + 3*(k+r-j-s) +2] = res[3*(p-j-s) +2];
  }
  for ( i = L+1; i < k - s; ++i) {
    newSpline[knotOffset + 3*i]   = res[3*(i-L)];
    newSpline[knotOffset + 3*i+1] = res[3*(i-L)+1];
    newSpline[knotOffset + 3*i+2] = res[3*(i-L)+2];
  }
  EG_free(res);
  EG_free(copySpline);
  header[3] = knotOffset;
  header[2] = header[3] - header[1]-1;
  (*spline) = EG_reall(*spline, (header[3] + header[2]*3)*sizeof(double));
  if ( *spline == NULL ) {
    EG_free(newSpline);
    return EGADS_MALLOC;
  }
  for ( i = 0; i < (header[3] + header[2]*3); ++i) {
    (*spline)[i] = newSpline[i];
  }
  EG_free(newSpline);
  return EGADS_SUCCESS;
}

static int
duplicateEndKnot(int *header,
    double **data, //spline= knots + cps GETS OVERWRITTEN !
    int locKnot,  // knot location
    int locCP)    // cp location
{
  int i,j;
  double *dataAux;
  dataAux = EG_alloc( (header[3] + 3* header[2])*sizeof(double) );
  if (dataAux == NULL ) {
    EG_free(dataAux);
    return EGADS_MALLOC;
  }
  for ( i = 0 ; i < (header[3] + 3*header[2]); ++ i) dataAux[i] = (*data)[i];
  *data = EG_reall(*data, (header[3] + 1 + 3* (header[2] +1))*sizeof(double));
  if ( *data == NULL ) {
    EG_free(dataAux);
    return EGADS_MALLOC;
  }
  for ( j = i = 0; i < header[3]; ++ i,++j) {
    (*data)[j] = dataAux[i];
    if ( i == locKnot){
      j++;
      (*data)[j] = dataAux[i];
    }
  }
  for ( j = i = 0; i < header[2]; ++ i,++j) {
    (*data)[ (header[3] + 1) + 3*j]    = dataAux[header[3] + 3*i];
    (*data)[ (header[3] + 1) + 3*j +1] = dataAux[header[3] + 3*i +1];
    (*data)[ (header[3] + 1) + 3*j +2] = dataAux[header[3] + 3*i +2];
    if ( i == locCP){
      j++;
      (*data)[ (header[3] + 1) + 3*j]    = dataAux[header[3] + 3*i];
      (*data)[ (header[3] + 1) + 3*j +1] = dataAux[header[3] + 3*i +1];
      (*data)[ (header[3] + 1) + 3*j +2] = dataAux[header[3] + 3*i +2];
    }
  }
  ++header[3];
  ++header[2];
  EG_free(dataAux);
  return EGADS_SUCCESS;
}

static int
EG_bisectBspline(ego context, ego bspline, double tVal, ego *piece)
{
  ego    geom;
  int    newHeaders[2][4], *header = NULL, ctype, oclass,  i,  mult, stat, offsetA, offsetB, span, inserts, initKnots, p1;
  double *data = NULL, *newData = NULL;
  stat = EG_getGeometry(bspline, &oclass, &ctype, &geom, &header, &data);
  if ( stat != EGADS_SUCCESS  || header == NULL || data == NULL) {
    printf(" Invalid curve for splitting... %d\n",stat);
    EG_free(data);
    EG_free(header);
    return stat;
  }
  p1        = 0;
  span      = FindSpan( header[2]-1, header[1], tVal, data);
  mult      = getMultiplicity(tVal, data, span, header[3]);
  inserts   = header[1] - mult;
  initKnots = header[3];
  stat      = EG_CurveKnotIns(header, &data, tVal, span, mult, inserts);
  // Now insert "dumb knot"
  if ( mult > 0 ) p1 = 1;
  stat = duplicateEndKnot(header, &data, span + header[1] - p1, span -header[1] + inserts);
  if ( stat != EGADS_SUCCESS) {
    EG_free(data);
    EG_free(header);
    return stat;
  }
  span     += header[1]+1 - p1;
  // First Spline
  newHeaders[0][0] = 0;
  newHeaders[0][1] = header[1];
  newHeaders[0][3] = span+1;
  newHeaders[0][2] = newHeaders[0][3] - newHeaders[0][1] - 1;
  newData          = (double*)EG_alloc( (newHeaders[0][3] + newHeaders[0][2]*3)*sizeof(double));
  if ( newData == NULL ) {
    EG_free(header);
    EG_free(data);
    return EGADS_MALLOC;
  }
  offsetB = header[3];
  for ( i = 0 ; i < newHeaders[0][3];  ++i) newData[i] = data[i];
  offsetA = newHeaders[0][3];
  offsetB = header[3];
  for ( i =  0; i < newHeaders[0][2]; ++i) {
    newData[offsetA + 3*i ]   = data[offsetB+3*i];
    newData[offsetA + 3*i +1] = data[offsetB+3*i+1];
    newData[offsetA + 3*i +2] = data[offsetB+3*i+2];
  }
  stat = EG_makeGeometry(context, CURVE, BSPLINE, NULL, newHeaders[0], newData, &piece[0]);
  if ( stat != EGADS_SUCCESS) {
    EG_free(newData);
    EG_free(header);
    EG_free(data);
    return stat;
  }
  // Second Spline
  newHeaders[1][0] = 0;
  newHeaders[1][1] = header[1];
  newHeaders[1][3] = header[3] - span + header[1];
  newHeaders[1][2] = newHeaders[1][3] - newHeaders[1][1] - 1;
  EG_free(newData);
  newData = (double*)EG_alloc( (newHeaders[1][3] + newHeaders[1][2]*3)*sizeof(double));
  if ( newData == NULL ) {
    EG_free(header);
    EG_free(data);
    return EGADS_MALLOC;
  }
  offsetB = newHeaders[0][3] - (header[1] + 1);
  for ( i = 0 ; i < newHeaders[1][3]  ; ++i)  newData[i] = data[offsetB + i];
  offsetA = newHeaders[1][3];
  offsetB = header[3] + 3*(header[2] - newHeaders[1][2]);
  for ( i =  0; i < newHeaders[1][2]; ++i) {
    newData[offsetA + 3*i ]   = data[offsetB+3*i];
    newData[offsetA + 3*i +1] = data[offsetB+3*i+1];
    newData[offsetA + 3*i +2] = data[offsetB+3*i+2];
  }
  stat = EG_makeGeometry(context, CURVE, BSPLINE, NULL, newHeaders[1], newData, &piece[1]);
  if ( stat != EGADS_SUCCESS) {
    EG_free(newData);
    EG_free(header);
    EG_free(data);
    return stat;
  }
  EG_free(newData);
  EG_free(data);
  EG_free(header);
  return stat;
}

int
EG_BsplineSplit(ego context, ego bspline, int nC, double t_cut[], ego *piece)
{
  int n, stat;
  stat = EG_bisectBspline(context, bspline, t_cut[0], &piece[0]);
  if ( stat !=EGADS_SUCCESS) {
    printf(" EG_bisectBspline %d !!!\n",stat);
    return stat;
  }
  for ( n = 1; n < nC -1; ++n) {
    stat = EG_bisectBspline(context, piece[n], t_cut[n], &piece[n]);
    if ( stat != EGADS_SUCCESS){
      return stat;
    }
  }
  return EGADS_SUCCESS;
}



#ifdef SPLINESPLIT

int main(/*@unused@*/ int argc, /*@unused@*/ char *argv[])
{

  clock_t start_t, end_t;
  int     nPt, nCP, nC, degree, stat,  i, j, nIters, stepControl;
  ego     context, curve, *pieces = NULL;
  double  *xyz = NULL, t, fitting_tol, rms, *tVals = NULL, *err1 = NULL, *err2 = NULL, delta;
  double  range[2], dPtInv, sum, point[3];
  int     zero = 0, offset;
  char    buffer[32];
  FILE    *fIn = NULL;

  start_t = clock();
#ifdef EGADS_DEBUG
  feenableexcept(FE_INVALID | FE_OVERFLOW);
#endif
  nC          = 3;
  nPt         = 50;
  degree      = 3;
  nCP         = degree + 1;
  fitting_tol = 1.E-08;
  delta       = 1.e-08;
  nIters      = 100;
  printf(" EG_open          = %d\n", EG_open(&context));
  // Reading of number of points from file (if input filename given)
  xyz  = (double *) EG_alloc(3* nPt * sizeof(double));
  err1 = (double *) EG_alloc(3* nPt * sizeof(double));
  err2 = (double *) EG_alloc(3* nPt * sizeof(double));
  if (xyz == NULL || err1 == NULL || err2 == NULL) {
    printf("MALLOC ERROR IN MAIN  !! ");
    return EGADS_MALLOC;
  }
  for (i = 0; i < nPt; i++) {
    t = (double)i/(double)(nPt-1);
    xyz[3*i]   = t;
    xyz[3*i+1] = sin(PI*t);
    xyz[3*i+2] = cos(2.0*PI*t);
  }
  stat = EG_1dBsplineCurveFit(context, xyz, nPt, &curve, nCP, degree, &rms, fitting_tol, stepControl, nIters, delta);
  if (stat == EGADS_SUCCESS) {
    fIn = fopen("splineSplit.dat", "w");
    if ( fIn == NULL ) {
      printf(" Failed to open file ...");
      stat =  EGADS_EMPTY;
      goto cleanup;
    }
    dPtInv = 1.0 / ((double) (nPt - 1));
    for (i = 0; i < nPt; i++) {
      t    = ((double) i) * dPtInv;
      stat = EG_evaluate(curve, &t, point);
      err1[3*i]   = point[0];
      err1[3*i+1] = point[1];
      err1[3*i+2] = point[2];
      if (stat == EGADS_SUCCESS)
        fprintf(fIn, "%22.14le %22.14le %22.14le %22.14le\n", err1[3*i], err1[3*i+1], err1[3*i+2],  t);
      else {
        printf(" stat at EG_evaluate %d\n",stat);
        fclose(fIn);
        goto cleanup;
      }
    }
    pieces = (ego*)EG_alloc(nC*sizeof(ego));
    tVals  = (double*)EG_alloc((nC-1)*sizeof(double));
    if (pieces == NULL || tVals == NULL) {
      printf(" in main: MALLOC !!!!!!\n");
      stat = EGADS_MALLOC;
      goto cleanup;
    }
    for ( i = 0; i < nC-1; ++ i)  tVals[i] = (i+1)*(1.0)/(double)(nC);
    stat = EG_BsplineSplit(context, curve, nC, tVals, pieces);
    if ( stat != EGADS_SUCCESS){
      printf(" EG_BsplineSplit = %d!!\n",stat);
      goto cleanup;
    }
    zero   = 0;
    offset = 0;
    for ( i = 0; i < nC; ++ i) {
      snprintf(buffer, sizeof(char) * 32, "spline%i.txt", i);
      fIn = fopen(buffer, "w");
      if ( fIn == NULL ) {
        printf(" Failed to open file ...");
        goto cleanup;
      }
      stat = EG_getRange(pieces[i], range, &zero);
      for (j = 0; j < nPt; j++) {
        t = ((double) j) * dPtInv;
        if ( (t >= range[0]) && (t <= range[1]) ) {
          stat = EG_evaluate(pieces[i], &t, point);
          err2[3*offset]   = point[0];
          err2[3*offset+1] = point[1];
          err2[3*offset+2] = point[2];
          if ( stat != EGADS_SUCCESS) {
            printf(" EG_evaluate  = %d\n",stat);
            goto cleanup;
          }
          fprintf(fIn, "%22.14le %22.14le %22.14le %22.14le\n", err2[3*offset], err2[3*offset+1], err2[3*offset+2],t);
          ++offset;
          if ( offset == nPt) break;
        }
      }
      fclose(fIn);
    }
    sum = 0.0;
    for ( i = 0; i < nPt; ++i) {
      sum += ( fabs( err1[i*3] - err2[i*3] ) + fabs( err1[i*3+1] - err2[i*3+1] ) + fabs( err1[i*3+2] - err2[i*3+2] ) );
    }
    printf(" ERROR IN B SPLINE SPLIT %22.14le\n",sum);
    printf(" JUNCTIONS:\n");
    for ( i = 0; i < nC - 1; ++ i) {
      stat = EG_evaluate(pieces[i], &tVals[i], point);
      printf(" t %lf  Curve %d  %.15lf\t%.15lf\t%.15lf\t",tVals[i], i, point[0],point[1],point[2]);
      stat = EG_evaluate(pieces[i+1], &tVals[i], point);
      printf("Curve %d  %.15lf\t%.15lf\t%.15lf\n",i+1, point[0],point[1],point[2]);
    }
  }
  cleanup :
  EG_free(xyz);
  EG_free(pieces);
  EG_free(tVals);
  EG_free(err1);
  EG_free(err2);
  printf("\n");
  EG_close(context);
  end_t = clock();
  printf(" \n TIME  = %lf  segs\n",(double)(end_t - start_t) / CLOCKS_PER_SEC);
  return stat;


}

#endif
















#ifdef BSPLINEAPPROXIMATION

int main(/*@unused@*/ int argc, /*@unused@*/ char *argv[])
{
  clock_t start_t, end_t;
  int     nPt, nCP, degree, stat, nChar = 120, i, *splineHeader = NULL, nIters, eD, stepControl;
  ego     context, curve;
  double  rms, point[3], *xyz = NULL, *splineData = NULL, dPtInv, t, fitting_tol, delta;
  char    fileIn[nChar], fileOut[nChar], line[nChar];
  FILE    *fIn = NULL;
  if (argc < 5) {
    printf(" Usage is  name input file + name output file  + field dimension + tolerance (NCP)\n");
    return 1;
  }
  start_t = clock();

#ifdef EGADS_DEBUG
  feenableexcept(FE_INVALID | FE_OVERFLOW);
#endif
  eD          = atoi(argv[3]);
  degree      = 3;
  sscanf(argv[4], "%lf", &fitting_tol);
  stepControl = 1;
  if ( argc == 6) {
    nCP         = atoi(argv[5]);//degree + 1;
    nIters      = 0;
  }
  else {
    nCP         = degree+ 1;
    nIters      = 100;
  }

  delta       = 1.e-08;//fitting_tol*0.01;   // Controls when the algorithm has "stalled" ie,  inserting a new knot improves the L2 error by delta
  strcpy(fileIn,  argv[1]);
  strcpy(fileOut, argv[2]);
  printf(" EG_open          = %d\n", EG_open(&context));
  // Reading of number of points from file (if input filename given)
  fIn = fopen(fileIn, "r");
  if ( fIn == NULL ) {
    printf(" Failed to open file %s...",fileIn);
    return EGADS_EMPTY;
  }
  fscanf(fIn, "%s %d\n", line, &nPt);
  printf("  Fitting for %d data points\n",nPt);
  // Allocation of memory for points
  xyz = (double *) EG_alloc(3 * nPt * sizeof(double));
  if (xyz == NULL) {
    printf("MALLOC ERROR IN MAIN  !! ");
    return EGADS_MALLOC;
  }
  if (eD == 2) {
    for (i = 0; i < nPt; i++) {
      int posit      = 3 * i;
      fscanf(fIn, "%lf %lf \n", &(xyz[posit]), &(xyz[posit + 1]));
      xyz[posit+2] = 0.0;
    }
  }else {
    for (i = 0; i < nPt; i++) {
      int posit      = 3 * i;
      fscanf(fIn, "%lf %lf %lf\n", &(xyz[posit]), &(xyz[posit + 1]), &(xyz[posit + 2]));

    }
  }
  fclose(fIn);
  stat = EG_1dBsplineCurveFit(context,xyz, nPt, &curve, nCP, degree, &rms, fitting_tol, stepControl, nIters, delta);
  printf(" STAT =  %d\n",stat);
  if (stat == EGADS_SUCCESS) {
    fIn = fopen(fileOut, "w");
    if ( fIn == NULL ) {
      printf(" Failed to open file %s...",fileIn);
      EG_free(xyz);
      return EGADS_EMPTY;
    }
    dPtInv = 1.0 / ((double) (nPt - 1));
    for (i = 0; i < nPt; i++) {
      t    = ((double) i) * dPtInv;
      stat = EG_evaluate(curve, &t, point);
      fprintf(fIn, "%22.14le %22.14le %22.14le %22.14le \n", point[0], point[1], point[2], t);
    }
  }
  EG_free(splineData);
  EG_free(splineHeader);
  EG_free(xyz);
  EG_close(context);
  end_t = clock();
  printf(" TIME  = %lf  segs\n",(double)(end_t - start_t) / CLOCKS_PER_SEC);
  return (0);
}

#endif
