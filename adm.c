#include <math.h>
#include "nr.h"
#include "nrutil.h"
#define PRTERR 0
#define MAXITS 1000000
#define NRMAX 10
#define ERRMAX 1000000
#define FREE_RETURN {free_vector(xhist,flag,n*(NRMAX+1)+shift);\
			free_vector(dhist,flag,n*(NRMAX+1)+shift);free_vector(errhist,0, ERRMAX-1);\
			free_matrix(u,1,NRMAX,1,NRMAX); free_matrix(b,1,NRMAX,1,1); free_vector(xorig,flag,n+shift);return xnew;}
//#include "CG.h"
void
gaussj (double **a, int n, double **b, int m);

/**********************************************************************************************
 Anderson mixing method to solve an equation x = f(x)
 flag==0: x[0,n-1];
 flag==1: x[1,n];
 n: total number of variables of x[];
 The dimension of the variables of funcvmix must be the same as x[];
 must assign a function like funcvmix(n,x,xnew) to get the new values of xnew from an old guess x
 **********************************************************************************************/

double*
adm (double *x, int n, int *check, void
(*funcvmix) (int, double *x, double *xnew),
     int flag)
{
  double err = 1e-10;
  int i, ib, ie, j, k, j1, k1, its = 1, nr = 1, shift, nc = 0, nnc, SMP = 0;
  double *xnew, *xorig, **u, **b, *xhist, *dhist, *errhist, lambda, tmp, tmp1,
      TOLF = err;
  int errmeth = 1;

  shift = flag - 1;
  xnew = vector (flag, n + shift);
  xorig = vector (flag, n + shift);
  xhist = vector (flag, n * (NRMAX + 1) + shift);
  dhist = vector (flag, n * (NRMAX + 1) + shift);
  errhist = vector (0, ERRMAX - 1);
  u = matrix (1, NRMAX, 1, NRMAX);
  b = matrix (1, NRMAX, 1, 1);
  lambda = 0.05;
  for (i = flag; i <= n + shift; ++i)
    xorig[i] = x[i];

  goto ADM;
  //goto SPM;

  /* Optional to use simple mixing */
  SPM: for (i = flag; i <= n + shift; ++i)
    x[i] = xorig[i];
  funcvmix (n, x, xnew);
  printf ("Start simple mixing...\n");
  lambda = 0.1;
  its = 0;
  for (;;)
    {
      ++its;
      for (i = flag; i <= n + shift; ++i)
	x[i] = (1 - lambda) * x[i] + lambda * xnew[i];
      funcvmix (n, x, xnew);
      if (errmeth == 1)
	{
	  for (tmp = 0.0, i = flag; i <= n + shift; ++i)
	    {
	      tmp1 = fabs (xnew[i] - x[i]);
	      tmp = tmp > tmp1 ? tmp : tmp1;
	    }
	  err = tmp;
	}
      else
	{
	  for (tmp = tmp1 = 0.0, i = flag; i <= n + shift; ++i)
	    {
	      tmp += SQR(xnew[i] - x[i]);
	      tmp1 += SQR(x[i]);
	    }
	  if (tmp1 == 0)
	    tmp1 = 1E-10;
	  err = sqrt (tmp / tmp1);
	}
      if (its % 1 == 0)
	printf ("simpmix: its= %d, err= %16.12e\n", its, err);
      if (err < 1E-3)
	{
	  printf ("simpmix: its= %d, err= %16.12e\nSwitch to anderson mixing\n",
		  its, err);
	  SMP = 1;
	  break;
	}
      if (its > 200)
	{
	  printf (
	      "simpmix: err = %16.12e at its= %d\nSwitch to anderson mixing\n",
	      err, its);
	  SMP = 1;
	  break;
	}
    }

  /* Anderson mixing */
  ADM: its = 1;
  lambda = 0.05;
  funcvmix (n, x, xnew);

  for (i = flag; i <= n + shift; ++i)
    {
      xhist[i] = x[i];
      //printf("xnew[%d]=%f\n",i,xnew[i]);
      dhist[i] = xnew[i] - x[i];
    }
  nc = 1;
  if (errmeth == 1)
    {
      for (tmp = 0.0, i = flag; i <= n + shift; ++i)
	{
	  tmp1 = fabs (dhist[i]);
	  tmp = tmp > tmp1 ? tmp : tmp1;
	}
      err = tmp;
    }
  else
    {
      for (tmp = tmp1 = 0.0, i = flag; i <= n + shift; ++i)
	{
	  tmp += SQR(dhist[i]);
	  tmp1 += SQR(x[i]);
	}
      if (tmp1 == 0)
	tmp1 = 1E-10;
      err = sqrt (tmp / tmp1);
    }

  printf ("andersmix: its= %d, err= %16.12e\n", its, err);

  if (err < TOLF)
    {
      //printf("andersmix: its= %d, err= %16.12e\n", its, err);
      *check = 0;

      FREE_RETURN
    }

  for (i = flag; i <= n + shift; ++i)
    x[i] = xhist[i] + lambda * dhist[i];

  for (its = 2; its <= MAXITS; ++its)
    {
      nr = IMIN(its - 1, NRMAX);
      lambda = 1.0 - pow (0.95, its);
      funcvmix (n, x, xnew);
      if (nc == NRMAX + 1)
	nc = 0;
      nnc = n * nc;
      ib = nnc + flag;
      ie = ib + n - 1;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //if (its % 10 == 0)
      //printf("andersmix: its= %d, err= %16.12e\n", its, err);
      //for (i = 0; i < n; i++)
      //printf("%d,xnew= %22.14e, x= %22.14e error=%22.14e  \n",i, xnew[i], x[i],xnew[i]-x[i]);

      /////////////////////////////////////////////////////////////////////////////////////////////////////

      int debug;

      for (i = ib; i <= ie; ++i)
	{
	  dhist[i] = xnew[i - nnc] - x[i - nnc];
	  xhist[i] = x[i - nnc];
	}

      debug = 11;
      ++nc;
      if (errmeth == 1)
	{
	  for (tmp = 0.0, i = ib; i <= ie; ++i)
	    {
	      tmp1 = fabs (dhist[i]);
	      tmp = tmp > tmp1 ? tmp : tmp1;
	    }
	  err = tmp;
	}
      else
	{
	  for (tmp = tmp1 = 0.0, i = ib; i <= ie; ++i)
	    {
	      tmp += SQR(dhist[i]);
	      tmp1 += SQR(xhist[i]);
	    }
	  if (tmp1 == 0)
	    tmp1 = 1E-10;
	  err = sqrt (tmp / tmp1);
	}

//		if (its % 1 == 0)
//			printf("andersmix: its= %d, err= %16.12e\n", its, err);

      if (err < TOLF)
	{

	  printf ("andersmix: its= %d, err= %16.12e\n", its, err);

	  *check = 0;

//			scanf("%d",&debug);
//			debug=100;
	  //printf("------------------------------------------\n");
	  //for (i = 0; i < n; i++)
	  //printf("%d,xnew= %22.14e, x= %22.14e error=%22.14e  \n",i, xnew[i], x[i],xnew[i]-x[i]);
	  printf (
	      "Solve sucessfully !!!!!!!!!!!! error =%21.14e,TOLF=%21.14e!!!!!!!!!!!! \n",
	      err, TOLF);

	  FREE_RETURN
	}
      if (its >= ERRMAX * 2)
	errhist[its % ERRMAX] = err;
      if (its >= ERRMAX * 3)
	{
	  i = its % ERRMAX;
	  if (errhist[i] > errhist[i == ERRMAX - 1 ? 0 : i + 1] || err >= 1E-4)
	    {
	      *check = 1;
	      printf (
		  "solve failed  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
	      FREE_RETURN
	    }
	}

      for (k = nc - 1; k >= 1; --k)
	{
	  k1 = n * (nc - k);
	  for (j = k; j >= 1; --j)
	    {
	      j1 = n * (nc - j);
	      for (tmp = 0.0, i = ib; i <= ie; ++i)
		tmp += (dhist[i] - dhist[i - k1]) * (dhist[i] - dhist[i - j1]);
	      u[nc - k][nc - j] = u[nc - j][nc - k] = tmp;
	    }
	  for (j = nr + 1; j >= nc + 1; --j)
	    {
	      j1 = n * (j - nc);
	      for (tmp = 0.0, i = ib; i <= ie; ++i)
		tmp += (dhist[i] - dhist[i - k1]) * (dhist[i] - dhist[i + j1]);
	      u[nc - k][nr + 1 + nc - j] = u[nr + 1 + nc - j][nc - k] = tmp;
	    }
	}
      for (k = nr + 1; k >= nc + 1; --k)
	{
	  k1 = n * (k - nc);
	  for (j = k; j >= nc + 1; --j)
	    {
	      j1 = n * (j - nc);
	      for (tmp = 0.0, i = ib; i <= ie; ++i)
		tmp += (dhist[i] - dhist[i + k1]) * (dhist[i] - dhist[i + j1]);
	      u[nr + 1 + nc - k][nr + 1 + nc - j] = u[nr + 1 + nc - j][nr + 1
		  + nc - k] = tmp;
	    }
	}

      for (j = nc - 1; j >= 1; --j)
	{
	  j1 = n * (nc - j);
	  for (tmp = 0.0, i = ib; i <= ie; ++i)
	    tmp += (dhist[i] - dhist[i - j1]) * dhist[i];
	  b[nc - j][1] = tmp;
	}
      for (j = nr + 1; j >= nc + 1; --j)
	{
	  j1 = n * (j - nc);
	  for (tmp = 0.0, i = ib; i <= ie; ++i)
	    tmp += (dhist[i] - dhist[i + j1]) * dhist[i];
	  b[nr + 1 + nc - j][1] = tmp;
	}

      gaussj (u, nr, b, 1);
      if (0)
	{
	  printf ("Adm: gaussj: Singular Matrix\n");
	  *check = 1;
	  printf (
	      "solve failed  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
	  FREE_RETURN
	    //exit(1);
	}

      for (i = ib; i <= ie; ++i)
	{
	  tmp = 0.0;
	  tmp1 = 0.0;
	  for (j = nc - 1; j >= 1; --j)
	    {
	      j1 = n * (nc - j);
	      tmp += b[nc - j][1] * (xhist[i - j1] - xhist[i]);
	      tmp1 += b[nc - j][1] * (dhist[i - j1] - dhist[i]);
	    }
	  for (j = nr + 1; j >= nc + 1; --j)
	    {
	      j1 = n * (j - nc);
	      tmp += b[nr + 1 + nc - j][1] * (xhist[i + j1] - xhist[i]);
	      tmp1 += b[nr + 1 + nc - j][1] * (dhist[i + j1] - dhist[i]);
	    }
	  x[i - nc * n + n] = xhist[i] + tmp + lambda * (dhist[i] + tmp1);
	}
    }
  printf ("MAXITS exceeded in Anderson mixing\n");
  *check = 1;
  printf (
      "solve failed  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
  FREE_RETURN
}

#undef FREE_RETURN
#undef MAXITS
#undef NRMAX
#undef ERRMAX
