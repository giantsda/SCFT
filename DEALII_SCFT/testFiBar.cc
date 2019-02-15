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
  double tau = 5.30252230020752e-01, L = 3.72374357332160;
  int N = pow (2, 16) + 1;
  std::vector<double> f0_given;

  f0_given.resize (N);
  for (int i = 0; i < N; i++)
    f0_given[i] = 1.;
  std::vector<double> x;
  x.resize (N);
  for (int i = 0; i < N; i++)
    {
      x[i] = L / (N - 1) * i;
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

  double f0Bar = romint (&f0_given[0], N - 1, L / (N - 1)) / L;
  printf ("f0Bar=%2.15f;N=%d\n", f0Bar, N);

  FILE *file;
  file = fopen ("inputFiles/Exp_m32_n2048_IE.res", "r");
  if (file == NULL)
    {
      fprintf (stderr, "Can't open input file!\n");
      exit (1);
    }
  char buff[255];
  int i;
  double value;
  N = 33;
  std::vector<double> yita_middle;
  yita_middle.resize (N, 0.);
  std::vector<double> xp (N, 0.);

  for (int i = 0; i < 9; i++)
    fgets (buff, 255, file);

  i = 0;
  while (fgets (buff, 255, file))
    {
      double xtmp;
      sscanf (buff, "%lf  %*lf  %lf", &xtmp, &value);
      xp[i] = xtmp * L;
      yita_middle[i] = value;
      i++;
    }

//  f0_given.resize (N);
  f0_given.assign (N, 1.0);
  for (int i = 0; i < N; i++)
    {
      x[i] = xp[i];
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

  for (int i = 0; i < N; i++)
    {
      yita_middle[i] = yita_middle[i] * f0_given[i];
    }

  double mean_field_free_energy = romint (&yita_middle[0], N - 1, L / (N - 1));

  mean_field_free_energy = (mean_field_free_energy / f0Bar / L + log (f0Bar))
      / (-1000.);
  printf ("integration=%2.15f\n", mean_field_free_energy);
  fclose (file);
  return 0;
}
