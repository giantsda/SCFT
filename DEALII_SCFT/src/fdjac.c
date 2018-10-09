# include <stdio.h>
#include <math.h>
#define NRANSI
#include "nrutil.h"
#define EPS 1.0e-7

extern int funcerr;


void fdjac(int n, float x[], float fvec[], float **df,
	void (*vecfunc)(int, float [], float []))
{
	int i,j;
	float h,temp,*f;

/*printf ("Calculating Jacobi...  "); */
	f=vector(1,n);
	for (j=1;j<=n;j++)
	{
		temp=x[j];
		h=EPS*temp;
		if (fabs(h) < EPS)  h = SIGN(EPS,temp);
		x[j]=temp+h;
		h=x[j]-temp;
		(*vecfunc)(n,x,f);
		if (funcerr)  goto retn;
		x[j]=temp;
		for (i=1;i<=n;i++) df[i][j]=(f[i]-fvec[i])/h;
	}
/*printf ("Done.\n");*/

/*printf ("Jacobian Matrix:\n");
output2d (df, n, n);
*/
retn:
	free_vector(f,1,n);
}
#undef EPS
#undef NRANSI
