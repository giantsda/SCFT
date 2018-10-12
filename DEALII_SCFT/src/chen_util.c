/*
 * chen_util.c
 *
 *  Created on: Oct 12, 2018
 *      Author: chen
 */
#include "NR_chen.h"
#include <stdio.h>
#include <stdlib.h>

double **
Matcreate (int r, int c) // The elements in the rows are next to each other.
{
  double** A = (double **) malloc (sizeof(double *) * r);
  A[0] = (double *) malloc (sizeof(double) * c * r);
  if (A == NULL || A[0] == NULL)
    {
      printf ("Matcreate:ERROR: out of Memory!\n");
      exit (-1);
    }
  for (int i = 0; i < r; i++)
    A[i] = (*A + c * i);
  return A;
}

void
Matfree (double** A)
{
  free (A[0]);
  free (A);
}
