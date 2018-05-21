#include <stdio.h>
#include <petscksp.h>
#include <petscdm.h>
#include <petscdmda.h>
#include <petscsnes.h>
#include <stdlib.h>
#include <math.h>
#include "nr.h"
#include "nrutil.h"

void
myfun (int n, double *x, double *xnew);
double*
adm (double *x, int n, int *check, void
(*myfun) (int, double *x, double *xnew),
     int flag);

double **
Matcreate (int r, int c)
{
  double** A = (double **) malloc (sizeof(double *) * r);
  A[0] = (double *) malloc (sizeof(double) * c * r);
  for (int i = 0; i < r; i++)
    A[i] = (*A + c * i);
  return A;
}

void
Matfree (double** A, int r, int c)
{
  free (A[0]);
  free (A);
}

int
simple_FEM_1D_transient (int N, double* yita_middle, double * f0,
			 double* f0_given)
{
  // N is the total number of nodes in 1D problem.
  KSP ksp;
  PC pc;
  Mat A, B, C, D;
  Vec X, b, xold;
  MPI_Comm comm;
  PetscErrorCode ierr;

  ierr = PetscInitialize (NULL, NULL, 0, NULL);
  if (ierr)
    return ierr;
  comm = MPI_COMM_SELF;

  double L = 1., dt = 0.01, t = 0, h = L / (N - 1); // L is the length of the space domain
  double solution[N], yita[N];
  int de;
  // add two fixed point of yita
  yita[0] = 0.;
  yita[N] = 0.;
  for (int i = 1; i < N - 1; i++)
    {
      yita[i] = yita_middle[i - 1];
    }
//  for (int i = 0; i < N; i++)
//    printf ("%f\n", yita[i]);
//  scanf ("%d", &de);

  int time_step = 10;

  PetscInt ix[N];
  double** solution_store = Matcreate (N + 1, time_step + 1);
  ierr = MatCreateSeqAIJ (comm, N, N, 3, NULL, &A);
  ierr = MatCreateSeqAIJ (comm, N, N, 3, NULL, &B);
  ierr = MatCreateSeqAIJ (comm, N, N, 3, NULL, &C);
  CHKERRQ(ierr);
  ierr = VecCreateSeq (comm, N, &b);
  CHKERRQ(ierr);
  ierr = VecDuplicate (b, &X);
  CHKERRQ(ierr);
  ierr = VecDuplicate (b, &xold);
  CHKERRQ(ierr);

  PetscScalar vforA[3] =
    { h / 6, 2. / 3 * h, h / 6 };
  PetscScalar vforB[3] =
    { -1 / h, 2. / h, -1 / h };

  for (PetscInt i = 0; i < N; i++)
    {
      PetscInt column_number[3] =
	{ i - 1, i, i + 1 };
      PetscScalar vforC[3] =
	{ h / 6 * yita[i], 2. / 3 * h * yita[i], h / 6 * yita[i] };
      if (i == 0)
	{
	  MatSetValues (A, 1, &i, 2, &column_number[1], &vforA[1],
			INSERT_VALUES);
	  MatSetValues (B, 1, &i, 2, &column_number[1], &vforB[1],
			INSERT_VALUES);
	  MatSetValues (C, 1, &i, 2, &column_number[1], &vforC[1],
			INSERT_VALUES);
	  VecSetValue (xold, i, 0., INSERT_VALUES);
	}
      else if (i == N - 1)
	{
	  MatSetValues (A, 1, &i, 2, column_number, vforA, INSERT_VALUES);
	  MatSetValues (B, 1, &i, 2, column_number, vforB, INSERT_VALUES);
	  MatSetValues (C, 1, &i, 2, column_number, vforC, INSERT_VALUES);
	  VecSetValue (xold, i, 0., INSERT_VALUES);
	}
      else
	{
	  MatSetValues (A, 1, &i, 3, column_number, vforA, INSERT_VALUES);
	  MatSetValues (B, 1, &i, 3, column_number, vforB, INSERT_VALUES);
	  MatSetValues (C, 1, &i, 3, column_number, vforC, INSERT_VALUES);
	  VecSetValue (xold, i, 1., INSERT_VALUES);
	}
    }

  for (PetscInt i = 0; i < N; i++)
    ix[i] = i;
  ierr = VecGetValues (xold, N, ix, solution);
  CHKERRQ(ierr);

  solution_store[0][0] = 0;
  for (int i = 1; i < N; i++)
    solution_store[i + 1][0] = solution[i];

  ierr = MatAssemblyBegin (A, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  ierr = MatAssemblyEnd (A, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  ierr = MatAssemblyBegin (B, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  ierr = MatAssemblyEnd (B, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  ierr = MatAssemblyBegin (C, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  ierr = MatAssemblyEnd (C, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  ierr = VecAssemblyBegin (xold);
  CHKERRQ(ierr);
  ierr = VecAssemblyEnd (xold);
  CHKERRQ(ierr);
  MatDuplicate (A, MAT_COPY_VALUES, &D);

//  printf ("A:\n");
//  ierr = MatView (A, PETSC_VIEWER_STDOUT_WORLD);
//  CHKERRQ (ierr);
//  printf ("B:\n");
//  ierr = MatView (B, PETSC_VIEWER_STDOUT_WORLD);
//  CHKERRQ (ierr);
//  printf ("C:\n");
//  ierr = MatView (C, PETSC_VIEWER_STDOUT_WORLD);
//  CHKERRQ (ierr);
//  printf ("D:\n");
//  ierr = MatView (D, PETSC_VIEWER_STDOUT_WORLD);
//  CHKERRQ (ierr);
//  printf ("xold:\n");
//  ierr = VecView (xold, 0);
//  CHKERRQ (ierr);

  // TODO: change this to a better time marching method.

  ierr = MatAXPY (D, dt, B, SAME_NONZERO_PATTERN);  // A=A+dt*B
  ierr = MatAXPY (D, dt, C, SAME_NONZERO_PATTERN);  // A=A+dt*C
  ierr = MatSetValue (D, 0, 0, 1.0, INSERT_VALUES);
  ierr = MatSetValue (D, 0, 1, 0.0, INSERT_VALUES);
  ierr = MatSetValue (D, N - 1, N - 1 - 1, 0.0, INSERT_VALUES);
  ierr = MatSetValue (D, N - 1, N - 1, 1.0, INSERT_VALUES);
  ierr = MatAssemblyBegin (D, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);
  ierr = MatAssemblyEnd (D, MAT_FINAL_ASSEMBLY);
  CHKERRQ(ierr);

//  printf ("D:\n");
//  ierr = MatView (D, PETSC_VIEWER_STDOUT_WORLD);

  ierr = KSPCreate (comm, &ksp);
  CHKERRQ(ierr);
  ierr = KSPSetOperators (ksp, D, D);
  CHKERRQ(ierr);
  ierr = KSPSetType (ksp, KSPCG);
  CHKERRQ(ierr);
  ierr = KSPSetInitialGuessNonzero (ksp, PETSC_TRUE);
  CHKERRQ(ierr);
  ierr = KSPGetPC (ksp, &pc);
  CHKERRQ(ierr);
  ierr = PCSetType (pc, PCICC);
  CHKERRQ(ierr);
  ierr = KSPSetFromOptions (ksp);
  CHKERRQ(ierr);
  ierr = KSPSetUp (ksp);
  CHKERRQ(ierr);

  for (int i = 1; i < time_step + 1; i++)
    {
      t = t + dt;
      printf ("time step:%d ; t= %f \n", i, t);
      ierr = MatMult (A, xold, b);  // b=A*xold;
      ierr = VecSetValue (b, 0, 0., INSERT_VALUES);
      ierr = VecSetValue (b, N - 1, 0., INSERT_VALUES);
      ierr = KSPSolve (ksp, b, X);
      VecAssemblyBegin (X);
      VecAssemblyEnd (X);
      CHKERRQ(ierr);
      ierr = VecGetValues (X, N, ix, solution);
      CHKERRQ(ierr);
      // store X to solution_store
      for (int j = 0; j < N; j++)
	{
	  solution_store[0][i] = t;
	  solution_store[j + 1][i] = solution[j];
	}
      VecCopy (X, xold);
    }

  ierr = KSPDestroy (&ksp);
  CHKERRQ(ierr);
  ierr = MatDestroy (&A);
  CHKERRQ(ierr);
  ierr = MatDestroy (&B);
  CHKERRQ(ierr);
  ierr = MatDestroy (&C);
  CHKERRQ(ierr);
  ierr = MatDestroy (&D);
  CHKERRQ(ierr);
  ierr = VecDestroy (&b);
  CHKERRQ(ierr);
  ierr = VecDestroy (&X);
  CHKERRQ(ierr);
  ierr = VecDestroy (&xold);
  CHKERRQ(ierr);

//   Plot results
  printf ("solution_store: \n");
  for (int i = 0; i < N + 1; i++)
    {
      for (int j = 0; j < time_step + 1; j++)
	printf ("%f, ", solution_store[i][j]);
      printf ("\n");
    }

//   integrate for f0
//   TODO:change this to better integration method

  printf ("f0: \n");
  for (int i = 0; i < N; i++)
    {
      f0[i] = 0.0;
      for (int j = 0; j < time_step; j++)
	{
	  double value_left = solution_store[i + 1][j]
	      * solution_store[i + 1][time_step - j];
	  double value_right = solution_store[i + 1][j + 1]
	      * solution_store[i + 1][time_step - j - 1];
	  f0[i] = f0[i] + 0.5 * dt * (value_left + value_right);
	}
      printf ("%f \n", f0[i]);
    }
  ierr = PetscFinalize ();
  printf ("have a nice day! \n");
  return ierr;
}

int
main ()
{
  printf ("Hello world!\n");
  int N = 5;
  double yita_middle[N - 2]; // initial guess, the ends are bounded
  for (int i = 0; i < N - 2; i++)
    {
      yita_middle[i] = i + 1;
    }
  double f0[N]; // results
  // calculate f0_given so that I can use andmix
  double tau = 0.2;
  double f0_given[N];
  double x[N];
  for (int i = 0; i < N; i++)
    {
      f0_given[i] = 1.;
      x[i] = i * 1. / (N - 1);
    }
  int x_left = ceil (N * tau / 1.);
  for (int i = 0; i < x_left; i++)
    {
      f0_given[i] = pow ((exp (4 * tau * x[i] / (tau * tau - x[i] * x[i])) - 1),
			 2)
	  / pow (((exp (4 * tau * x[i] / (tau * tau - x[i] * x[i])) + 1)), 2);
      f0_given[N - i - 1] = f0_given[i];
    }
//  for (int i = 0; i < N; i++)
//    {
//      printf (">>%f \n ", f0_given[i]);
//    }
  simple_FEM_1D_transient (N, yita_middle, f0, f0_given);







//  double in[2] =
//    { 0, 0 };
//  int check = 1;
//  double* haha = adm (in, 2, &check, myfun, 0);
//  for (int i = 0; i < 2; i++)
//    printf ("%f \n", haha[i]);
//  return 100;
}

void
myfun (int n, double * x, double * xnew)
{
  xnew[0] = x[1] * 0.5 - 2;
  xnew[1] = x[0] + 3;
  printf ("x=%.14f, xnew=%.14f \n", *x, *xnew);
}

