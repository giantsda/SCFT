#include <stdio.h>

#include <fstream>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include "nrutil.h"
#include "nr.h"
/* Use Romberg integration to integrate (not volume-average !) array f[]
 having m+1 (m=2^(M-1)) data points with a uniform grid spacing of hh */
//  g++ testFiBar.cc src/romint.c src/nrutil.c src/polint.c -I ./include/ -lm
int
main ()
{
  double tau = 0.5302, L = 3.72374;
  int N = pow (2, 22) + 1;
  std::vector<double> f0_given;

  f0_given.resize (N);
  for (int i = 0; i < N; i++)
    f0_given[i] = 1.;
  std::vector<double> x;
  x.resize (N);
  for (int i = 0; i < N; i++)
    {
      x[i] = 1. / (N - 1) * i;
    }

  for (int i = 0; i < N; i++)
    {
      if (x[i] <= tau)
	{
	  f0_given[i] = pow (
	      (exp (4 * tau * x[i] / (tau * tau - x[i] * x[i])) - 1), 2)
	      / pow (((exp (4 * tau * x[i] / (tau * tau - x[i] * x[i])) + 1)),
		     2);
	  if (std::isnan (f0_given[i]))
	    f0_given[i] = 1.;

	  f0_given[N - i - 1] = f0_given[i];

	}
      else
	break;
    }

//  for (int i = 0; i < N; i++)
//    printf ("f0_given[%d]=%2.15f\n", i, f0_given[i]);

  double f0Bar = romint (&f0_given[0], N - 1, 1. / (N - 1)) / L;
  printf ("f0Bar=%2.15f;N=%d\n", f0Bar, N);
  return 0;
}
