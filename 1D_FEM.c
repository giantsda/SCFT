#include <stdio.h>
#include <petscksp.h>

int
main ()
{
  KSP ksp;
  PC pc;
  Mat A, B, C, D;
  Vec X, b, xold, x_initial;
  MPI_Comm comm;
  //  PetscScalar        v;
  //  KSPConvergedReason reason;
  //  PetscInt           i,j,its;
  PetscErrorCode ierr;

  //  PetscFunctionBegin;
  ierr = PetscInitialize (NULL, NULL, 0, NULL);
  if (ierr)
    return ierr;
  comm = MPI_COMM_SELF;
  printf ("Hello");
  int N = 8; // N is the total number of nodes in 1D problem.
  double L = 1.9, tau = 0.4, dt = 0.01, t = 0, h = L / (N - 1); // L is the length of the space domain
  int i;
  int time_step = 7;
  double yita[N];
  double solution_store[N + 1][time_step];
  for (int i = 0; i < N; i++)
    {
      yita[i] = 1.;
    }

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
  ierr = VecDuplicate (b, &x_initial);
  CHKERRQ(ierr);

  double vforA[3] =
    { h / 6, 2. / 3 * h, h / 6 };
  double vforB[3] =
    { -1 / h, 2. / h, -1 / h };

  for (int i = 0; i < N; i++)
    {
      int column_number[3] =
	{ i - 1, i, i + 1 };
      double vforC[3] =
	{ h / 6 * yita[i], 2. / 3 * h * yita[i], h / 6 * yita[i] };
      printf ("i=%d\n", i);
      if (i == 0)
	{
	  MatSetValues (A, 1, &i, 2, &column_number[1], &vforA[1],
			INSERT_VALUES);
	  MatSetValues (B, 1, &i, 2, &column_number[1], &vforB[1],
			INSERT_VALUES);
	  MatSetValues (C, 1, &i, 2, &column_number[1], &vforC[1],
			INSERT_VALUES);
	  VecSetValue (xold, i, 1., INSERT_VALUES);
	}
      else if (i == N - 1)
	{
	  MatSetValues (A, 1, &i, 2, &column_number, &vforA, INSERT_VALUES);
	  MatSetValues (B, 1, &i, 2, &column_number, &vforB, INSERT_VALUES);
	  MatSetValues (C, 1, &i, 2, &column_number, &vforC, INSERT_VALUES);
	  VecSetValue (xold, i, 1., INSERT_VALUES);
	}
      else
	{
	  MatSetValues (A, 1, &i, 3, &column_number, &vforA, INSERT_VALUES);
	  MatSetValues (B, 1, &i, 3, &column_number, &vforB, INSERT_VALUES);
	  MatSetValues (C, 1, &i, 3, &column_number, &vforC, INSERT_VALUES);
	  VecSetValue (xold, i, 0., INSERT_VALUES);
	}
    }
  VecCopy (xold,x_initial);


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

  printf ("A:\n");
  ierr = MatView (A, PETSC_VIEWER_STDOUT_WORLD);
  CHKERRQ(ierr);
  printf ("B:\n");
  ierr = MatView (B, PETSC_VIEWER_STDOUT_WORLD);
  CHKERRQ(ierr);
  printf ("C:\n");
  ierr = MatView (C, PETSC_VIEWER_STDOUT_WORLD);
  CHKERRQ(ierr);
  printf ("D:\n");
  ierr = MatView (D, PETSC_VIEWER_STDOUT_WORLD);
  CHKERRQ(ierr);
  printf ("xold:\n");
  ierr = VecView (xold, 0);
  CHKERRQ(ierr);

for (i=0;i<=time_step;i++)
  {
    t=t+dt;






  }

























//  ierr = VecView (B, 0);
//  CHKERRQ (ierr);

//  ierr = KSPCreate (comm, &ksp);
//  CHKERRQ (ierr);
//  ierr = KSPSetOperators (ksp, A, A);
//  CHKERRQ (ierr);
//  ierr = KSPSetType (ksp, KSPCG);
//  CHKERRQ (ierr);
//  ierr = KSPSetInitialGuessNonzero (ksp, PETSC_TRUE);
//  CHKERRQ (ierr);
//  ierr = KSPGetPC (ksp, &pc);
//  CHKERRQ (ierr);
//  ierr = PCSetType (pc, PCICC);
//  CHKERRQ (ierr);
//  ierr = KSPSetFromOptions (ksp);
//  CHKERRQ (ierr);
//  ierr = KSPSetUp (ksp);
//  CHKERRQ (ierr);
//  ierr = KSPSolve (ksp, B, X);
//  CHKERRQ (ierr);
//  ierr = VecView (X, 0);
//  CHKERRQ (ierr);

  //  for (i=0; i<4; i++) {
  //    v    = 3;
  //    ierr = MatSetValues(A,1,&i,1,&i,&v,INSERT_VALUES);CHKERRQ(ierr);
  //    v    = 1;
  //    ierr = VecSetValues(B,1,&i,&v,INSERT_VALUES);CHKERRQ(ierr);
  //    ierr = VecSetValues(X,1,&i,&v,INSERT_VALUES);CHKERRQ(ierr);
  //  }
  //
  //  i    =0; v=0;
  //  ierr = VecSetValues(X,1,&i,&v,INSERT_VALUES);CHKERRQ(ierr);
  //
  //  for (i=0; i<3; i++) {
  //    v    = -2; j=i+1;
  //    ierr = MatSetValues(A,1,&i,1,&j,&v,INSERT_VALUES);CHKERRQ(ierr);
  //    ierr = MatSetValues(A,1,&j,1,&i,&v,INSERT_VALUES);CHKERRQ(ierr);
  //  }
  //  i=0; j=3; v=2;
  //
  //  ierr = MatSetValues(A,1,&i,1,&j,&v,INSERT_VALUES);CHKERRQ(ierr);
  //  ierr = MatSetValues(A,1,&j,1,&i,&v,INSERT_VALUES);CHKERRQ(ierr);
  //  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  //  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  //  ierr = VecAssemblyBegin(B);CHKERRQ(ierr);
  //  ierr = VecAssemblyEnd(B);CHKERRQ(ierr);
  //  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nThe Kershaw matrix:\n\n");CHKERRQ(ierr);
  //  ierr = MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  //
  //
  //  ierr = KSPCreate(comm,&ksp);CHKERRQ(ierr);
  //  ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
  //
  //  ierr = KSPSetType(ksp,KSPCG);CHKERRQ(ierr);
  //  ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);CHKERRQ(ierr);
  //
  //  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  //  ierr = PCSetType(pc,PCICC);CHKERRQ(ierr);
  //  /* ierr = PCFactorSetShiftType(prec,MAT_SHIFT_POSITIVE_DEFINITE);CHKERRQ(ierr); */
  //
  //  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
  //  ierr = KSPSetUp(ksp);CHKERRQ(ierr);
  //
  //  ierr = PCFactorGetMatrix(pc,&M);CHKERRQ(ierr);
  //  ierr = VecDuplicate(B,&D);CHKERRQ(ierr);
  //  ierr = MatGetDiagonal(M,D);CHKERRQ(ierr);
  //  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nPivots:\n\n");CHKERRQ(ierr);
  //  ierr = VecView(D,0);CHKERRQ(ierr);
  //
  //
  //  ierr = KSPSolve(ksp,B,X);CHKERRQ(ierr);
  //  ierr = KSPGetConvergedReason(ksp,&reason);CHKERRQ(ierr);
  //  if (reason==KSP_DIVERGED_INDEFINITE_PC) {
  //    ierr = PetscPrintf(PETSC_COMM_WORLD,"\nDivergence because of indefinite preconditioner;\n");CHKERRQ(ierr);
  //    ierr = PetscPrintf(PETSC_COMM_WORLD,"Run the executable again but with -pc_factor_shift_positive_definite option.\n");CHKERRQ(ierr);
  //  } else if (reason<0) {
  //    ierr = PetscPrintf(PETSC_COMM_WORLD,"\nOther kind of divergence: this should not happen.\n");CHKERRQ(ierr);
  //  } else {
  //    ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
  //    ierr = PetscPrintf(PETSC_COMM_WORLD,"\nConvergence in %d iterations.\n",(int)its);CHKERRQ(ierr);
  //  }
  //  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRQ(ierr);
  //
  //  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
  ierr = MatDestroy (&A);
  CHKERRQ(ierr);
  ierr = MatDestroy (&B);
  CHKERRQ(ierr);
  ierr = MatDestroy (&C);
  CHKERRQ(ierr);
  ierr = MatDestroy (&D);
  CHKERRQ(ierr);
//  ierr = MatDestroy (&B);
//  CHKERRQ (ierr);
  //  ierr = VecDestroy(&B);CHKERRQ(ierr);
  //  ierr = VecDestroy(&X);CHKERRQ(ierr);
  //  ierr = VecDestroy(&D);CHKERRQ(ierr);

  ierr = PetscFinalize ();
  return ierr;
}
