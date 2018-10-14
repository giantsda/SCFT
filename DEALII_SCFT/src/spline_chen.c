/*
 * spline_chen.cpp
 *
 *  Created on: Oct 11, 2018
 *      Author: chen
 */
#include "NR_chen.h"
#include "nr.h"
#include "nrutil.h"
#include <stdlib.h>

void
spline_chen (double* x, double* y, double* xp, double* yp, int Nx, int Nxp,
	     double* m)
/* x,y is the table waiting to be interpolated.
 x from [0 Nx];
 xp is the x vector you want to know yp at.
 is m==NULL use not-a-knot end condition
 m=0 for natural spline;
 otherwise m is the ddy at boundarys
 */
{
  double** A = dmatrix (1, Nx, 1, Nx); // NR requires it 1 based.

  double** B = dmatrix (1, Nx, 1, 1);

  for (int i = 1; i < Nx + 1; i++)  // initialize A and b
    {
      for (int j = 1; j < Nx + 1; j++)
	A[i][j] = 0.;
      B[i][1] = 0.;
    }

  for (int i = 2; i < Nx; i++)
    {
      A[i][i - 1] = (x[i - 1] - x[i - 2]) / 6.;
      A[i][i] = (x[i] - x[i - 2]) / 3.;
      A[i][i + 1] = (x[i] - x[i - 1]) / 6.;
      B[i][1] = (y[i] - y[i - 1]) / (x[i] - x[i - 1])
	  - (y[i - 1] - y[i - 2]) / (x[i - 1] - x[i - 2]);
    }

  if (m == NULL) // use not-a-knot
    {
      A[1][1] = 1 / (x[1] - x[0]);
      A[1][2] = -1 / (x[1] - x[0]) - 1 / (x[2] - x[1]);
      A[1][3] = 1 / (x[2] - x[1]);
      A[Nx][Nx - 2] = 1 / (x[Nx - 2] - x[Nx - 3]);
      A[Nx][Nx - 1] = -1 / (x[Nx - 2] - x[Nx - 3])
	  - 1 / (x[Nx - 1] - x[Nx - 2]);
      A[Nx][Nx] = 1 / (x[Nx - 1] - x[Nx - 2]);
    }
  else
    {
      A[1][1] = 1.;
      A[Nx][Nx] = 1.;
      B[1][1] = *m;
      B[Nx][1] = *m;
    }

  int fail = 0;
  fail = gaussj (A, Nx, B, 1);
  if (fail)
    {
      printf ("spline_chen:ERROR! gaussj cannot solve Ax=b\n");
      exit (-1);
    }
  int k;
  double h, a, b;

  for (int i = 0; i < Nxp; i++)
    {
      int klo = 0;
      int khi = Nx - 1;
      while (khi - klo > 1)
	{
	  k = (khi + klo) >> 1;
	  if (x[k] > xp[i])
	    khi = k;
	  else
	    klo = k;
	}
      h = x[khi] - x[klo];
      if (h == 0.0)
	{
	  printf ("spline_chen:ERROR! Bad x input x should be increasing\n");
	  exit (-1);
	}
      a = (x[khi] - xp[i]) / h;
      b = (xp[i] - x[klo]) / h;
      yp[i] = a * y[klo] + b * y[khi]
	  + ((a * a * a - a) * B[klo + 1][1] + (b * b * b - b) * B[khi + 1][1])
	      * (h * h) / 6.0;

    }
}

