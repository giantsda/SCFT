# include <stdio.h>
#include <math.h>
#define NRANSI
#include "nrutil.h"

#define MAXITS 400
/* MAXITS is the maximum number of iterations */
#define EPS 1e-14
/* EPS is a number close to machine precision */
#define TOLX EPS
#define STPMX 100.0
#define FREERETURN {free_vector(fvec,1,n);free_vector(xold,1,n);\
	free_vector(w,1,n);free_vector(t,1,n);free_vector(s,1,n);\
	free_vector(p,1,n);free_vector(g,1,n);free_vector(fvcold,1,n);\
	free_vector(c,1,n);return;}

int nn;
float *fvec;
void
(*funcvbrd) (int n, float v[], float f[]);

/* The following must be defined in the calling program */
extern double **qt, **r, *d, /* qt[1:n,1:n], r[1:n,1:n] and d[1:n] must be allocated in the calling program */
err; /* Passed in as the convergence criterion, and returned with the actual residual error */
extern int funcerr, /* Flag for error in evaluating the function in Broyden method */
jc; /* For re-use of Jacobian. jc denotes the status of Jacobian calculation: 0 for not calculated,
 1 for previously calculated, 2 for currently calculated */
extern int PRINT;

float
fminbrd (float x[])
{
  int i;
  float sum;

  (*funcvbrd) (nn, x, fvec);
  if (funcerr)
    return 0;
  for (sum = 0.0, i = 1; i <= nn; i++)
    sum += SQR(fvec[i]);
  return 0.5 * sum;
}

void
broydn (float x[], int n, int *check, void
(*vecfunc) (int, float[], float[]))
{
  void
  fdjac (int n, float x[], float fvec[], float **df, void
  (*vecfunc) (int, float[], float[]));
  void
  lnsrch (int n, float xold[], float fold, float g[], float p[], float x[],
  float *f,
	  float stpmax, int *check, float
	  (*func) (float[]));
  void
  qrdcmp (float **a, int n, float *c, float *d, int *sing);
  void
  qrupdt (float **r, float **qt, int n, float u[], float v[]);
  void
  rsolv (float **a, int n, float d[], float b[]);
  int i, its, j, k, restrt = 0, sing, skip, ii;
  float den, f, fold, stpmax, sum, temp, *c, *fvcold, test;
  float *g, *p, *s, *t, *w, *xold;
  double TOLF = err, /* The maximum absolute value on function values for convergence */
  TOLMIN = TOLF; /* Used to determine if spurious convergence to a minimum of fminbrd has occurred */

  c = vector (1, n);
  fvcold = vector (1, n);
  g = vector (1, n);
  p = vector (1, n);
  s = vector (1, n);
  t = vector (1, n);
  w = vector (1, n);
  xold = vector (1, n);
  fvec = vector (1, n);
  nn = n;
  funcvbrd = vecfunc;
  /* PRINT = 1; */
  f = fminbrd (x); /* Function values are calculated */
  PRINT = 0;
  if (funcerr)
    goto error;
  err = 0.0;
  for (i = 1; i <= n; i++)
    if (fabs (fvec[i]) > err)
      {
	err = fabs (fvec[i]);
	ii = i;
      }
  printf ("   Broydn 0: |f|max = %e, imax=%d\n", err, ii);
  if (err < TOLF)
    {
      *check = 0;
      FREERETURN
    }
  for (sum = 0.0, i = 1; i <= n; i++)
    sum += SQR(x[i]);
  stpmax = STPMX * FMAX(sqrt (sum), (float)n);
  if (!jc)
    restrt = 1;

  for (its = 1; its <= MAXITS; its++)
    {
      /* Form Q^T and R in current step */
      if (restrt)
	{
	  fdjac (n, x, fvec, r, vecfunc);
	  printf ("r=:\n");
	  for (int i = 1; i <= n; i++)
	    {
	      for (int j = 1; j <= n; j++)
		printf ("%f,", r[i][j]);
	      printf ("\n");
	    }
	  int de;
//	  printf("continue?\n");
//	  scanf("%d",&de);
	  if (funcerr)
	    goto error;
	  qrdcmp (r, n, c, d, &sing);
	  if (sing)
	    {
	      printf ("It is singlular \n ");
//	      scanf ("%d", &de);
	    }
	  if (sing)
	    nrerror ("singular Jacobian in broydn");
	  for (i = 1; i <= n; i++)
	    {
	      for (j = 1; j <= n; j++)
		qt[i][j] = 0.0;
	      qt[i][i] = 1.0;
	    }
	  for (k = 1; k < n; k++)
	    {
	      if (c[k])
		{
		  for (j = 1; j <= n; j++)
		    {
		      sum = 0.0;
		      for (i = k; i <= n; i++)
			sum += r[i][k] * qt[i][j];
		      sum /= c[k];
		      for (i = k; i <= n; i++)
			qt[i][j] -= sum * r[i][k];
		    }
		}
	    }
	  for (i = 1; i <= n; i++)
	    {
	      r[i][i] = d[i];
	      for (j = 1; j < i; j++)
		r[i][j] = 0.0;
	    }
	  jc = 2;
	}
      else if (its > 1)
	{
	  for (i = 1; i <= n; i++)
	    s[i] = x[i] - xold[i];
	  for (i = 1; i <= n; i++)
	    {
	      for (sum = 0.0, j = i; j <= n; j++)
		sum += r[i][j] * s[j];
	      t[i] = sum;
	    }
	  skip = 1;
	  for (i = 1; i <= n; i++)
	    {
	      for (sum = 0.0, j = 1; j <= n; j++)
		sum += qt[j][i] * t[j];
	      w[i] = fvec[i] - fvcold[i] - sum;
	      if (fabs (w[i]) >= EPS * (fabs (fvec[i]) + fabs (fvcold[i])))
		skip = 0;
	      else
		w[i] = 0.0;
	    }
	  if (!skip)
	    {
	      for (i = 1; i <= n; i++)
		{
		  for (sum = 0.0, j = 1; j <= n; j++)
		    sum += qt[i][j] * w[j];
		  t[i] = sum;
		}
	      for (den = 0.0, i = 1; i <= n; i++)
		den += SQR(s[i]);
	      for (i = 1; i <= n; i++)
		s[i] /= den;
	      qrupdt (r, qt, n, t, s);
	      for (i = 1; i <= n; i++)
		{
		  if (r[i][i] == 0.0)
		    nrerror ("r singular in broydn");
		  d[i] = r[i][i];
		}
	    }
	}

      /* Find suitable step size */
      for (i = 1; i <= n; i++)
	{
	  for (sum = 0.0, j = 1; j <= n; j++)
	    sum += qt[i][j] * fvec[j];
	  p[i] = -sum;
	}
      for (i = n; i >= 1; i--)
	{
	  for (sum = 0.0, j = 1; j <= i; j++)
	    sum -= r[j][i] * p[j];
	  g[i] = sum;
	}
      for (i = 1; i <= n; i++)
	{
	  xold[i] = x[i];
	  fvcold[i] = fvec[i];
	}
      fold = f;
      rsolv (r, n, d, p);
      /*PRINT = 1; */
      lnsrch (n, xold, fold, g, p, x, &f, stpmax, check, fminbrd);
      PRINT = 0;

      err = 0.0;
      for (i = 1; i <= n; i++)
	if (fabs (fvec[i]) > err)
	  err = fabs (fvec[i]);
      if (its % 1 == 0)
	printf ("   Broydn %d: |f|max = %e\n", its, err);
      if (err < TOLF)
	{
	  printf ("   Broydn %d: |f|max = %e\n", its, err);
	  *check = 0;
	  jc = 1;
	  FREERETURN
	}

      if (*check)
	{ /* Line search failed to find a new x[] */
	  if (restrt)
	    FREERETURN
	      /* Broyden failed, return with jc=2 */
	  else
	    {
	      test = 0.0;
	      den = FMAX(f, 0.5 * n);
	      for (i = 1; i <= n; i++)
		{
		  temp = fabs (g[i]) * FMAX(fabs (x[i]), 1.0) / den;
		  if (temp > test)
		    test = temp;
		}
	      if (test < TOLMIN)
		{
		  printf ("   Broydn %d: |f|max = %e   |df/dx|max = %e\n", its,
			  err, test);
		  (*check) = 0;
		  jc = 1;
		  FREERETURN
		}
	      else
		restrt = 1;
	    }
	}
      else
	{
	  restrt = 0;
	  test = 0.0;
	  for (i = 1; i <= n; i++)
	    {
	      temp = (fabs (x[i] - xold[i])) / FMAX(fabs (x[i]), 1.0);
	      if (temp > test)
		test = temp;
	    }
	  if (test < TOLX)
	    {
	      printf ("   Broydn %d: |f|max = %e   |dx|max = %e\n", its, err,
		      test);
	      jc = 1;
	      FREERETURN
	    }
	}
    }

  error:
  /*	nrerror("MAXITS exceeded in broydn"); */
  funcerr = 0;
  *check = 1;
  FREERETURN
    /* with jc=0, 1, or 2 */
}

#undef MAXITS
#undef EPS
#undef TOLX
#undef STPMX
#undef FREERETURN
#undef NRANSI
