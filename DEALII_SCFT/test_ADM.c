#include "NR_chen.h"
#include <stdio.h>
#include <math.h>
#include "nr.h"
#include "nrutil.h"
#include <stdlib.h>

void
myfun (int n, double* in, double* out)
{
  double x = in[0];
  double y = in[1];
  double z = in[2];
  out[0] = x * y * z - 12.;
  out[1] = x * x + y * y - 8.;
  out[2] = x + y + z - 511.;

}

int
main ()
{
  double x[] =
    { 1, 2, 3 };
  double y[3];

  adm_chen (myfun, x, 1e-12, 50000, 3);

  return 0;

}
