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

#include "egads.h"
#include "egadsSplineFit.h"
#include <time.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#define PI  3.1415926535897931159979635
#define           MIN(A,B)        (((A) < (B)) ? (A) : (B))
#define           MAX(A,B)        (((A) < (B)) ? (B) : (A))
#define DEPS 1.e-12

#define MIN_KNOT_DIST 1.e-4
#define MAX_KNOTS 1000
#define EPS2 1.e-07
#define EPS1 1.e-07
#define NEWTON_IT_MAX 50

#define SAVE_OPTI_AND_CP 1
#define SAVE_ITER_NORMS 1

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
		eps2    = fabs(f)/(arc_der*eps1);
		if (eps1 <= EPS1 || arc_der < DEPS || eps2 <= EPS2 || fabs(df) < DEPS) {
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


static int
checkSpan(double u0,double u1,double u2, double *t, int m, int *n)
{
	int i, spanOK = 0, u1u2Empty = 0;
	double alpha = 0.9;
	if (fabs ( u1-u2) < DEPS )
	{
		u1u2Empty = 1;
	}
	for ( i = *n; i < m; ++i)
	{
		if (spanOK == 0 ) {
			if( t[i] >= u0 ) {
				if( t[i] > u1 ) // knotSpan is empty
				{
					//	printf(" MODIFYING T VAL \n %lf\t",t[i]);
					t[i] = (1.0 - alpha)*u0 + alpha*u1;
					//	printf("%lf\n",t[i]);
				}
				spanOK = u1u2Empty + 1;
			}
		} else {
			if ( t[i] >= u1 ) {
				if ( t[i] > u2 ) // knotSpan is empty
				{
					//printf(" MODIFYING T VAL \n %lf\t",t[i]);
					t[i] = (1.0 - alpha)*u1 + alpha*u2;
					//	printf("%lf\n",t[i]);
				}
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

#define MIN_PTS_PER_SPAN 2

static int
get_max_error_knot(int offset, // = degree (start looking only at the interior knots)
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
	int     i, j, checkKnot,  first_t;
	knotPool     = (int*)   EG_alloc(nSpans*sizeof(int));
	ids          = (int*)   EG_alloc(nSpans*sizeof(int));
	knotErr      = (double*)EG_alloc(nSpans*sizeof(double));
	weightedErr  = (double*)EG_alloc(nSpans*sizeof(double));
	if (knotPool == NULL || knotErr == NULL || weightedErr == NULL || ids == NULL ) return EGADS_MALLOC;
	for ( i = 0; i < nSpans; ++ i) {
		knotPool[i] =   0;
		ids[i]      =   0;
		checkKnot   =   1;
		first_t     =  -1;
		knotErr[i]  = 0.0;
		errMax      = 0.0;
		for ( j = 0; j < *totBadKnots; ++j) { // check if we already have used that knot
			if ( (i+offset) == badKnotsVec[j])  {
				checkKnot = 0;
				break;
			}
		}
		if ( (checkKnot == 1) && (fabs( bSpline[i+offset] - bSpline[i+1+offset] ) >= 2.0*MIN_KNOT_DIST) ) {
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
		}
		if (knotPool[i] >=  MIN_PTS_PER_SPAN )	{  // min points per span.
			weightedErr[i] = sqrt ((knotErr[i])*(double)knotPool[i]/(double)(m-2));
		}
		else {
			weightedErr[i] = -1.0;
			badKnotsVec[*totBadKnots] = i + offset;
			++(*totBadKnots);
		}
	}
	*errId = -1;
	errMax = 0.0;
	for ( i = 0; i < nSpans; ++ i) {
		if ( weightedErr[i] > errMax ) {
			sum  = 0.0;
			tAux = bSpline[i+offset];
			for ( j = 0 ; j < knotPool[i] - 1; ++j) {
				sum += error[ids[i]+j]*error[ids[i]+j];
				if ( sum > knotErr[i]*0.5) {
					tAux = 0.5*( t[ids[i]+j] + t[ids[i]+j+1]);
					j    = m;
				}
			}
			if (  (fabs(tAux      - bSpline[i + offset])   >= MIN_KNOT_DIST )
					&& (fabs(tAux - bSpline[i + offset+1]) >= MIN_KNOT_DIST )) {
				*errId = i;
				errMax = weightedErr[i];
				*tVal  = tAux;
			}
		}
	}
	EG_free(knotErr);
	EG_free(weightedErr);
	EG_free(knotPool);
	EG_free(ids);
	if (*errId == -1) {
		printf(" Total Bad Spans %d. Possible spans %d. I can't refine anywhere\n",*totBadKnots, nSpans);
		return EGADS_EMPTY;
	}
	return EGADS_SUCCESS;
}



static intcreateKnotSpan( int header[], double *bSpline, int m, double t[], double XYZ[], double *error,
		int *spanLocation,  double  *newKnot, double *inL2error, double updateTOL)
{
	double *errAtIt = NULL, *bestError = NULL, *localSpline = NULL, *bestSpline = NULL;
	double  prevL2error, L2error = 0.0, LInferror =0.0, bestL2error , bestNewKnot, tVal;
	int     nKnots, nCP, deg, stat, totBadSpans, badSpanVec[MAX_KNOTS], nSpans;
	int     i,j, bestSpanLocation,  tPos, v, knotCount ;
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
	totBadSpans = 0;
	bestL2error = L2_INIT;
	do {
		stat = get_max_error_knot (deg, nSpans, badSpanVec, &totBadSpans, bSpline, error, t,  m, &tPos, &tVal);
		if (stat == EGADS_SUCCESS) {
			for ( j = 0 ; j <= deg; ++j)
			{  // clamped: copy first & last
				localSpline[j]             = bSpline[j];
				localSpline[nKnots -1 - j] = bSpline[ (nKnots -1) -1 -j];
			}
			v         = 1;
			knotCount = 0;	// checks if the inserted knot does not produce any problems for the LS-Solver
			for ( j = i = 0 ; i < nSpans; ++i, ++j) {
				if ( i != tPos ) {
					stat = checkSpan(bSpline[i + deg], bSpline[i + deg + 1], bSpline[i + deg +1], t, m-1, &v);
					if (stat != EGADS_SUCCESS ) {
						i = nKnots;
						break ;
					}
					localSpline[j + deg + 1] = bSpline[i + deg + 1];
					++knotCount;
				} else {
					localSpline[j+ deg + 1]  = tVal;
					stat = checkSpan(bSpline[i+deg], localSpline[j+deg+1], bSpline[i + deg + 1], t, m-1, &v);
					if (stat != EGADS_SUCCESS ) {
						i = nKnots;
						break;
					}
					++knotCount;
					++j;
					localSpline[j + deg + 1] = bSpline[i + deg + 1];
				}
			}
			if ( knotCount == nSpans ) {
				stat = EG_Linear_LS_Bspline_Optimization(m, XYZ, t, header, localSpline,  errAtIt, &L2error);
				if (stat != EGADS_SUCCESS) {
					printf("Solution To Linear System = %d!!!\n",stat);
					goto cleanup;
				}
				if ( L2error < prevL2error  && fabs ( L2error - prevL2error) >= updateTOL ) {
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
		}
		badSpanVec[totBadSpans] = tPos +deg;
		++totBadSpans;
	}
	while ( totBadSpans < nSpans );
	for ( i = 0; i < header[3] + 3*header[2]; ++i ) bSpline[i] = bestSpline[i];
	for ( i = 0 ; i < m; ++i) error[i] = bestError[i];
	*spanLocation = bestSpanLocation;
	*newKnot      = bestNewKnot;
	*inL2error    = bestL2error;
	stat          = EGADS_SUCCESS;
	cleanup:
	EG_free(localSpline);
	EG_free(bestSpline);
	EG_free(errAtIt);
	EG_free(bestError);
	return stat;
}


static int
get_minimal_knots (int header[], double spline[], int m, double *error, double t[], double XYZ[], int totNewKnots, double knotSeq[], int *outHeader, double **outSpline, double *oriL2)
{
	int    stat, i, j, nKnots, deg,  nSpans, pos;
	double newL2error,  *knotCopy = NULL, *newSpline = NULL;
	knotCopy   = (double*)EG_alloc((header[3] + totNewKnots )*sizeof(double));
	newSpline  = (double*)EG_alloc( (3*header[2]+header[3])*sizeof(double));
	if (knotCopy == NULL || newSpline == NULL) return EGADS_MALLOC;
	nKnots        = header[3];
	deg           = header[1];
	for ( i = 0; i < nKnots; ++i ) knotCopy[i] = spline[i];
	pos = totNewKnots - 1;
	do
	{
		++nKnots;
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
		stat = EG_Linear_LS_Bspline_Optimization(m, XYZ, t, outHeader, newSpline, error, &newL2error);
		if (stat != EGADS_SUCCESS ) {
			printf("EG_Linear_LS_Bspline_Optimization %d\n",stat);
			EG_free(knotCopy);
			EG_free(newSpline);
			return stat;
		}
		if ( newL2error <= *oriL2 ) {
			EG_free(knotCopy);
			*outSpline = (double*)EG_reall(*outSpline, (3*outHeader[2] + outHeader[3])*sizeof(double));
			for ( i =0; i < 3*outHeader[2] + outHeader[3]; ++i ) (*outSpline)[i] = newSpline[i];
			*oriL2     = newL2error;
			EG_free(newSpline);
			return EGADS_SUCCESS;
		}
		--pos;
		for ( i = 0 ; i < nKnots; ++ i ) knotCopy[i] = newSpline[i];
	}
	while(pos >= 0);
	outHeader[3] = outHeader[3] + totNewKnots;
	outHeader[2] = outHeader[3] - deg - 1;
	*outSpline = (double*)EG_reall(*outSpline, (3*outHeader[2] + outHeader[3])*sizeof(double));
	for ( i =0; i < 3*outHeader[2] + outHeader[3]; ++i ) (*outSpline)[i] = newSpline[i];
	EG_free(knotCopy);
	EG_free(newSpline);
	return EGADS_SUCCESS;
}


#define MAX_FALSE_ITERS 10
static int
EG_1dBsplineCurveFit(ego context, double *XYZcloud, int m, ego *curve, int n,
		int deg, double *rms, double LS_tol, int maxIter)
{
	double updateTOL, L2error, prevL2error, backL2error;
	double  *t = NULL,  *bSpline = NULL,  *bSplineCopy = NULL, *backSpline = NULL, *XYZ = NULL, *error= NULL, badKnotSeq[MAX_KNOTS];
	double point[3], scale, xmin, xmax, ymin, ymax, zmin, zmax, xcent, ycent, zcent,  newKnot;
	int    falseIter, normIter, headerCopy[4], header[4], backHeader[4];
	int    rejections, accepted, reductions, spanLocation,  badCount;
	int    iter, stat, i, j, k, nKnots;
	int    paramType  = 0; // 1 = chord length, otherwise it uses centripetal
	FILE   *fn, *fbSpline;
	ego    bSplineCurve;
	nKnots = n + deg + 1;
	header[0] = 0 ;
	header[1] = deg;
	header[2] = n;
	header[3] = nKnots;
	if (XYZcloud == NULL ) {
		printf(" No Data...\n");
		return EGADS_NODATA;
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
	scale = 1.0 / MAX(MAX(xmax-xmin, ymax-ymin), zmax-zmin);
	if ( scale < DEPS ) {
		printf("This data is not suitable for fitting.\n");
		EG_free(XYZ);
		return EGADS_EMPTY;
	}
	xcent = 0.5* (xmin + xmax);
	ycent = 0.5* (ymin + ymax);
	zcent = 0.5* (zmin + zmax);
	for (k = 0; k < m; k++) {
		XYZ[3*k  ] = scale * (XYZcloud[3*k  ] - xcent);
		XYZ[3*k+1] = scale * (XYZcloud[3*k+1] - ycent);
		XYZ[3*k+2] = scale * (XYZcloud[3*k+2] - zcent);
	}
	// ALLOCATE MEMO
	t             = (double*)  EG_alloc (m*sizeof(double));	// Curve parameter
	bSpline       = (double*)  EG_alloc ( (nKnots+ n*3)*sizeof(double));
	bSplineCopy   = (double*)  EG_alloc ( (nKnots+ n*3)*sizeof(double));
	backSpline    = (double*)  EG_alloc ( (nKnots+ n*3)*sizeof(double));
	error         = (double*)  EG_alloc (m *sizeof(double));
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
	rejections    = 0;
	accepted      = 0;
	reductions    = 0;
	headerCopy[0] = 0;
	headerCopy[1] = deg;
	stat = EG_Linear_LS_Bspline_Optimization(m, XYZ, t, header, bSpline, error, &L2error);
	if (stat != EGADS_SUCCESS) {
		printf("Solution To Linear System = %d!!!\n",stat);
		goto cleanup;
	}
	for ( i = 0 ; i < header[3] + 3*header[2]; ++i) backSpline[i] = bSpline[i];
	for ( i = 0; i < 4; ++i) backHeader[i] = header[i];
	prevL2error = L2error;
	updateTOL   = L2error*0.01;
	falseIter   = 0;
	normIter    = 0;
	spanLocation = -1;
	badCount     = 0;
	while ( (L2error > LS_tol) && (iter < maxIter) && (n < m) ){
		++iter;
		headerCopy[3] = header[3] + 1;
		headerCopy[2] = headerCopy[3]-deg-1;
		bSplineCopy   = (double*)EG_reall(bSplineCopy, (headerCopy[3]+ headerCopy[2]*3)*sizeof(double));

		for ( i = 0; i < header[3]; ++i ) {
			bSplineCopy[i] = bSpline[i];
		}
		stat  = createKnotSpan(headerCopy, bSplineCopy, m, t, XYZ, error, &spanLocation, &newKnot, &L2error, updateTOL);
		if ( stat != EGADS_SUCCESS ) {
			printf(" Cannot refine more the Bspline. Leave Fitting with 'best' approximation %d\n",stat);
			for ( i = 0 ; i < 4; ++i) header[i] = backHeader[i];
			bSpline     = (double*)EG_reall(bSpline,(header[3]+header[2]*3)*sizeof(double));
			for ( i = 0; i < header[3] + header[2]*3; ++i ) bSpline[i] = backSpline[i];
			stat = EGADS_SUCCESS;
			break;
		}
		else {
			for ( i = 0 ; i < 4; ++i) header[i] = headerCopy[i];
			bSpline = (double*)EG_reall(bSpline,(header[3]+header[2]*3)*sizeof(double));
			for ( i = 0; i < header[3] + header[2]*3; ++i ) bSpline[i] = bSplineCopy[i];
		}
		if ( L2error > prevL2error  && (normIter == 0 )) {
			++rejections;
			++falseIter;
			updateTOL = 0.0;
			badKnotSeq[badCount] = newKnot;
			++badCount;
		}
		else if ( fabs(L2error - prevL2error) < updateTOL) {
			++reductions;
			++normIter;
			if (normIter > 10) updateTOL = 0.0;
			if (falseIter == 0) {
				badKnotSeq[badCount] = newKnot;
				++badCount;
			}
		} else {
			if (falseIter >= 1 || normIter >=1 ) {
				badKnotSeq[badCount] = newKnot;
				++badCount;
				stat = get_minimal_knots(backHeader, backSpline, m, error, t, XYZ, badCount, badKnotSeq, header, &bSpline, &prevL2error);
				if (stat != EGADS_SUCCESS) {
					printf(" stat GET_MINIMAL_KNOTS %d!!\n",stat);
					goto cleanup;
				}
			}
			++accepted;
			updateTOL    = L2error*0.01;
			falseIter    = 0 ;
			normIter     = 0 ;
			spanLocation = -1;
			badCount     = 0 ;
			prevL2error  = L2error;
			backL2error  = L2error;
#if SAVE_ITER_NORMS
			fprintf(fn, "%d\t%.15lf\n",iter,L2error);
#endif
			for ( i = 0 ; i < 4; ++i)	backHeader[i] = header[i];
			backSpline  = (double*)EG_reall(backSpline,(header[3]+header[2]*3)*sizeof(double));
			for ( i = 0; i < header[3] + header[2]*3; ++i )	backSpline[i] = bSpline[i];
			while( fabs(log10(prevL2error/updateTOL)) < 1.0) updateTOL*=0.1;
		}
		if ( falseIter == MAX_FALSE_ITERS ) {
			for ( i = 0 ; i < 4; ++i) header[i] = backHeader[i];
			bSpline = (double*)EG_reall(bSpline,(3*backHeader[2]+ backHeader[3])*sizeof(double));
			for ( i = 0 ; i < 3*backHeader[2]+ backHeader[3]; ++ i)	bSpline[i] = backSpline[i];
			L2error = backL2error;
			printf(" I have done 10 iterations with no improvement. Leave Fitting with 'best' approximation %d\n",stat);
			stat = EGADS_SUCCESS;
			break;
		}
	}
	printf(" Leaving Solver with norm %11.4e < %11.4e\n TOTALS: accepted %d rejected %d reduced %d \n",L2error, LS_tol, accepted, rejections, reductions);
#if SAVE_ITER_NORMS
	fclose(fn);
#endif
	if (stat != EGADS_SUCCESS) goto cleanup;
	for (j = 0; j < header[2]; j++) {
		bSpline[3*j   + header[3]] = xcent + bSpline[3*j   + header[3]] / scale;
		bSpline[3*j+1 + header[3]] = ycent + bSpline[3*j+1 + header[3]] / scale;
		bSpline[3*j+2 + header[3]] = zcent + bSpline[3*j+2 + header[3]] / scale;
	}
	stat   = EG_makeGeometry(context, CURVE, BSPLINE, NULL, header, bSpline, &bSplineCurve);
	if (stat != EGADS_SUCCESS) {
		printf(" EG_makeGeometry %d\n",stat);
		goto cleanup;
	}
	*curve   = bSplineCurve;
	*rms     = L2error/scale;
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
	fclose(fbSpline);
	fbSpline = fopen("CONTROL_POLYGON","w");
	if( fbSpline == NULL) {
		stat = EGADS_MALLOC;
		goto cleanup;
	}
	for ( j = 0; j < header[2]; ++j) {
		fprintf(fbSpline,"%.11f\t%.11f\t%.11f\n", bSpline[3*j   + header[3]], bSpline[3*j   + header[3]+1], bSpline[3*j   + header[3]+2]);
	}
	fclose(fbSpline);
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
			"---------------------------------------\n",stat,header[2],L2error/scale,iter,maxIter);
	return  stat;
}


int main(/*@unused@*/ int argc, /*@unused@*/ char *argv[])
{
	clock_t start_t, end_t;
	int     nPt, nCP, degree, stat, nChar = 120, i, *splineHeader = NULL, nIters;
	ego     context, curve;
	double  rms, point[3], *xyz = NULL, *splineData = NULL, dPtInv, t, fitting_tol;
	char    fileIn[nChar], fileOut[nChar], line[nChar];
	FILE    *fIn = NULL;
	if (argc != 3) {
		printf(" Usage is  name input file + name output file \n");
		return 1;
	}
	start_t = clock();

	degree      = 3;
	nCP         = degree + 4;
	fitting_tol = 1.E-07;
	nIters      = 100;
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
	for (i = 0; i < nPt; i++) {
		int posit      = 3 * i;
		//fscanf(fIn, "%lf %lf%lf\n", &(xyz[posit]), &(xyz[posit + 1]), &(xyz[posit + 2]));
		fscanf(fIn, "%lf %lf\n", &(xyz[posit]), &(xyz[posit + 1]));
		xyz[posit+2] = 0.0;
	}
	fclose(fIn);
	stat = EG_1dBsplineCurveFit(context,xyz, nPt, &curve, nCP, degree, &rms, fitting_tol, nIters);
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
			fprintf(fIn, "%22.14le %22.14le %22.14le %22.14le\n", point[0], point[1], point[2], t);
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
