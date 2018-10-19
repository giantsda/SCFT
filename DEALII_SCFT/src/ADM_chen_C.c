/*
 ============================================================================
 Name        : ADM_chen_C.c
 Author      : Chen
 Version     :
 Copyright   : Your copyright notice
 Description : Andmixing C, Ansi-style
 ============================================================================
 */

#include "NR_chen.h"
#include <stdio.h>
#include <math.h>
#include "nr.h"
#include "nrutil.h"
#include <stdlib.h>

int
adm_chen (void
(*f) (int, double* in, double* out),
	  double* x_old, double tol, int maxIteration, int n, double lmd,
	  int nn)
/* It solves a function of type: void f (double* in, double* out, int n, struct parameterDumper* p) p is for some parameters you would like to pass
 * so that global variables can be avoided!
 * x_old is the initial guess; tol is the tolerance; maxIteration is the max iteration number you allowed.
 * lmd is the relaxzation factor
 * nn is the max size of matrix U
 *
 * How to use it:
 * to solve myfun:
 *
 * double x[] = { 1, 2, 3 };
 *   struct parameterDumper p;
 *   p.double_a = 2.5;
 *   p.double_b = 3.6;
 *   int fail = adm_chen (&myfun, x, 1e-15, 3000, 3,0.99,30);
 */

{
  // n is the size of the vector problem
  int k = 0; // kth iteration
  double lk = lmd;
  int m;
  int nm = DMIN(nn, n);
  double** U = dmatrix (1, nm, 1, nm);
  double** V = dmatrix (1, nm, 1, 1);
  int k_restart = 0; // k_restart is used when U is ill.
  double err = 9.9e99; // err
  double** X = Matcreate (maxIteration + 2, n);
  double** Y = Matcreate (maxIteration + 2, n); /* X is used to store the guessed solution and
   Y is the resulted rhs*/
  for (int i = 0; i < n; i++)
    X[0][i] = x_old[i];

  while (err > tol && k <= maxIteration)
    {
      f (n, X[k], Y[k]);
      double err = 0.; // evaluate err
      for (int i = 0; i < n; i++)
	{
	  if (fabs (Y[k][i]) >= err)
	    err = fabs (Y[k][i]);
	}
      printf ("And_chen:iter=%d;N=%d,lk=%e,err=%2.4E\n", k, n, lk, err);
      if (err < tol)
	{
	  for (int i = 0; i < n; i++)
	    x_old[i] = X[k][i];
	  printf (
	      "*****And_chen: Solved equation successfully!*****\nThe solution is:\n");
	  for (int i = 0; i < n; i++)
	    printf ("X[%d]=%2.15f\n", i, x_old[i]);
	  free_dmatrix (U, 1, nm, 1, nm);
	  free_dmatrix (V, 1, nm, 1, 1);
	  Matfree (X);
	  Matfree (Y);
	  return 0;
	}

      restart: m = DMIN(nm, k - k_restart);
      // Constuct matrix U and vector V

      for (int i = 0; i < m; i++)
	{
	  for (int j = 0; j < m; j++)
	    {
	      U[i + 1][j + 1] = 0.;
	      for (int t = 0; t < n; t++)
		U[i + 1][j + 1] += (Y[k][t] - Y[k - i - 1][t])
		    * (Y[k][t] - Y[k - j - 1][t]);
	    }
	  V[i + 1][1] = 0.;
	  for (int t = 0; t < n; t++)
	    V[i + 1][1] += (Y[k][t] - Y[k - i - 1][t]) * Y[k][t];
	}

      int fail = 0;
      fail = gaussj (U, m, V, 1);

      if (fail == 1)
	{
	  printf ("And_chen: Singular Matrix detected And_chen restarted!\n");
	  k_restart = k;
	  lk = lmd;
	  goto restart;
	}
      // update new X;
      for (int i = 0; i < n; i++)
	{
	  double cx = 0., cd = 0.;
	  for (int j = 0; j < m; j++)
	    {
	      cx += V[j + 1][1] * (X[k - j - 1][i] - X[k][i]);
	      cd += V[j + 1][1] * (Y[k - j - 1][i] - Y[k][i]);
	    }
	  X[k + 1][i] = X[k][i] + cx + (1 - lk) * (Y[k][i] + cd);
	}

      if (err < 0.03 && k > 100) // only modifiy lk if it is close to solution
	lk *= lmd;

      if (lk < 1e-5)  // reset lk if it is too small
	lk = lmd;

      k++;
    }

  printf (
      "And_chen failed after %d iterations :(  Try to increase max iteration allowed\n",
      maxIteration);
  for (int i = 0; i < n; i++)
    x_old[i] = X[k][i];
  free_dmatrix (U, 1, nm, 1, nm);
  free_dmatrix (V, 1, nm, 1, 1);
  Matfree (X);
  Matfree (Y);
  return 1;
}

