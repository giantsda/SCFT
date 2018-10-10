# include <math.h>
# include "nrutil.h"
# include "nr.h"
# include <stdio.h>
/* # include <conio.h> */
# include <stdlib.h>
# include <math.h>
/* # include <alloc.h> */
# include <stdarg.h>
# include <string.h>
# define DIM0(NAME,SIZE,TYPE) { (NAME)=(TYPE  *) calloc((SIZE), sizeof(TYPE));  if ((NAME)==NULL)  { printf("\7Out of Memory !\n"); exit(1); } }
# define DIM(NAME,SIZE,TYPE) { (NAME)=(TYPE  *) malloc((SIZE)*sizeof(TYPE)); if ((NAME)==NULL)  { printf("\7Out of Memory !\n");  exit(1); } }
# define DO(i,s,e) for ((i)=(s); (i)<=(e); ++(i))
# define UNTIL(f) while (! (f));
# define JUMP_REM(fp) fgets(buffer,800,(fp));
 

# define K 5
/* K is the number of points used in the extrapolation */

double romint (double *f, int m, double hh)
/* Use Romberg integration to integrate (not volume-average !) array f[]
   having m+1 (m=2^(M-1)) data points with a uniform grid spacing of hh */
{
	double *s, *h, sum, ss, dss;
	int M, np, ii, i, j;

	M = (int) (log(m*1.0)/0.6931471805599453 + 1.5);
	if (M < K)
	{
		printf ("m must be >= 2^%d !\n", K-1);
		exit(1);
	}
	s = vector (1, M);
	h = vector (1, M);

	h[1] = 1.0;
	s[1] = m*hh*(f[0]+f[m])/2;

	np = 1;
	for (j=2; j<=M; j++)
	{
		ii = m/np;
		for (sum=0,i=ii/2; i<m; i+=ii)  sum += f[i];
		s[j] = (s[j-1]+ii*hh*sum)/2;
		np += np;
		h[j] = h[j-1]/4;
	}

	polint (&h[M-K], &s[M-K], K, 0.0, &ss, &dss);
/* for (i=1; i<=M; ++i)  printf ("h[%d]=%22.14e  s[%d]=%22.14e\n", i, h[i], i, s[i]); */
/* printf ("Romberg: ss=%22.14e  |dss|/|ss|=%e\n", ss, fabs(dss)/fabs(ss)); */
	free_vector (s,1,M);
	free_vector (h,1,M);

	return ss;
}
# undef K


double romavg (double *f, int m, int C)
/* Use Romberg integration to calculate the volume-average of array f[] having m+1 (m=2^(M-1)) data points
   in Cartesian (C=0), cylindrical (C=1), or spherical (C=2) coordinates */
{
	double *tmp, s;
	register int i;

	if (C==0)  return (romint(f,m,1.0/m));

	DIM (tmp, m+1, double)
	tmp[0] = 0;
	if (C==1)
		DO (i,1,m)  tmp[i] = f[i]*i/m;
	else
		DO (i,1,m)  tmp[i] = f[i]*i*i/(m*m);
	s = romint (tmp,m,1.0/m);
	free (tmp);

	return ((C+1)*s);
}

 
