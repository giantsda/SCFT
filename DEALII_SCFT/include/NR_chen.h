/*
 * NR_chen.h
 *
 *  Created on: Oct 6, 2018
 *      Author: chen
 */

#ifndef NR_CHEN_H_
#define NR_CHEN_H_

#include <stdio.h>

#include <stdlib.h>

struct parameterDumper
{
  int i1, i2, i3, i4, i5;
  double d1, d2, d3, d4, d5;
  int *pi1, *pi2, *pi3, *pi4, *pi5;
  double* pd1, *pd2, *pd3, *pd4, *pd5;
  char* s1, *s2, *s3, *s4, *s5;
};

double **
Matcreate (int r, int c); // The elements in the rows are next to each other.

void
Matfree (double** A);

int
adm_chen (void
(*f) (int, double* in, double* out),
	  double* x_old, double tol, int maxIteration, int n, double lmd,
	  int nn);

void
spline_chen (double* x, double* y, double* xp, double* yp, int Nx, int Nxp,
	     double* m);

#endif /* NR_CHEN_H_ */
