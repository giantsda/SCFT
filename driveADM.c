/*
 * driveADM.c
 *
 *  Created on: Oct 6, 2018
 *      Author: chen
 */
#include "NR_chen.h"

void
myfun (double* in, double* out, int n, struct parameterDumper* p)
{
  out[0] = in[0] * in[1] * in[2] - 12.;
  out[1] = in[1] * in[1] + in[0] * in[0] - 8.;
  out[2] = in[1] + in[0] + in[2] - 511.;

}

int
main ()
{
  double x[] =
    { 1, 2, 3 };

  struct parameterDumper p;
  p.d1 = 2.5;
  p.d2 = 3.6;
  int fail = adm_chen (&myfun, x, 1e-15, 3000, 3, &p);

  return 0;
}
