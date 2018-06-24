/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2013 - 2017 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Wolfgang Bangerth, Texas A&M University, 2013
 */

#include <deal.II/base/utilities.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/numerics/matrix_tools.h>

#include <fstream>
#include <iostream>

#include "nr.h"
#include "nrutil.h"

double* f0_given;
double **qt, **r, *d, /* qt[1:n,1:n], r[1:n,1:n] and d[1:n] must be allocated in the calling program */
err; /* Passed in as the convergence criterion, and returned with the actual residual error */
int funcerr, /* Flag for error in evaluating the function in Broyden method */
jc; /* For re-use of Jacobian. jc denotes the status of Jacobian calculation: 0 for not calculated,
 1 for previously calculated, 2 for currently calculated */
int PRINT;

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
Matfree (double** A)
{
  free (A[0]);
  free (A);
}

namespace Step26
{
  using namespace dealii;

  template<int dim>
    class HeatEquation
    {
    public:
      HeatEquation (int N, int total_time_step, double L, double* yita);
      double *
      run ();

    private:
      void
      setup_system ();
      void
      solve_time_step ();
      void
      output_results () const;
      void
      refine_mesh (const unsigned int min_grid_level,
		   const unsigned int max_grid_level);

      Triangulation<dim> triangulation;
      FE_Q<dim> fe;
      DoFHandler<dim> dof_handler;

      ConstraintMatrix constraints;

      SparsityPattern sparsity_pattern;
      SparseMatrix<double> mass_matrix;
      SparseMatrix<double> laplace_matrix;
      SparseMatrix<double> system_matrix;
      SparseMatrix<double> tmp;

      Vector<double> solution;
      Vector<double> old_solution;
      Vector<double> system_rhs;

      double time;
      double time_step;
      int timestep_number;
      int N, total_time_step;
      double** solution_store;
      double* f0;
      double L;
      double* yita_middle_1D;
      double* yita_full_1D;
      double* yita_full_2D;
      double* out;
    };

  template<int dim>
    class Initial_condition : public Function<dim>
    {
    public:
      Initial_condition () :
	  Function<dim> ()
      {
      }
      virtual double
      value (const Point<dim> &p, const unsigned int component = 0) const;
    };

  template<int dim>
    double
    Initial_condition<dim>::value (const Point<dim> &p,
				   const unsigned int component) const
    {
      (void) component;
      Assert(component == 0, ExcIndexRange(component, 0, 1));
      Assert(dim == 2, ExcNotImplemented());
      if (p[0] == 0. || p[0] == 3.72374)
	return 0.;
      else
	return 1.;
    }

  template<int dim>
    HeatEquation<dim>::HeatEquation (int N, int total_time_step, double L,
				     double* yita) :
	fe (1), dof_handler (triangulation), time (0.0), time_step (1. / 2047), timestep_number (
	    0), N (N), total_time_step (total_time_step), L (L), yita_middle_1D (
	    yita)
    {
      time_step = 1. / (total_time_step - 1);
      solution_store = Matcreate (N + 1, total_time_step + 1);
      f0 = (double*) malloc (sizeof(double) * N);
      yita_full_1D = (double*) malloc (sizeof(double) * N);
      yita_full_2D = (double*) malloc (sizeof(double) * 2 * N);
      out = (double*) malloc (sizeof(double) * (N - 1));
    }

  template<int dim>
    void
    HeatEquation<dim>::setup_system ()
    {
      dof_handler.distribute_dofs (fe);

      std::cout << std::endl << "==========================================="
	  << std::endl << "Number of active cells: "
	  << triangulation.n_active_cells () << std::endl
	  << "Number of degrees of freedom: " << dof_handler.n_dofs ()
	  << std::endl << std::endl;

      DynamicSparsityPattern dsp (dof_handler.n_dofs ());
      DoFTools::make_sparsity_pattern (dof_handler, dsp, constraints,
      /*keep_constrained_dofs = */true);
      sparsity_pattern.copy_from (dsp);

      mass_matrix.reinit (sparsity_pattern);
      laplace_matrix.reinit (sparsity_pattern);
      system_matrix.reinit (sparsity_pattern);
      tmp.reinit (sparsity_pattern);

      MatrixCreator::create_mass_matrix (dof_handler,
					 QGauss<dim> (fe.degree + 1),
					 mass_matrix);
      MatrixCreator::create_laplace_matrix (dof_handler,
					    QGauss<dim> (fe.degree + 1),
					    laplace_matrix);

      solution.reinit (dof_handler.n_dofs ());
      old_solution.reinit (dof_handler.n_dofs ());
      system_rhs.reinit (dof_handler.n_dofs ());

    }

  template<int dim>
    void
    HeatEquation<dim>::solve_time_step ()
    {
      SolverControl solver_control (8000, 1e-13);
      SolverCG<> solver (solver_control);
      solver.solve (system_matrix, solution, system_rhs,
		    PreconditionIdentity ());

//      std::cout << "     " << solver_control.last_step () << " CG iterations."
//	  << std::endl;
    }

  template<int dim>
    void
    HeatEquation<dim>::output_results () const
    {
      DataOut<dim> data_out;

      data_out.attach_dof_handler (dof_handler);
      data_out.add_data_vector (solution, "U");

      data_out.build_patches ();

      const std::string filename = "solution-"
	  + Utilities::int_to_string (timestep_number, 3) + ".vtk";
      std::ofstream output (filename.c_str ());
      data_out.write_vtk (output);
    }

  template<int dim>
    double *
    HeatEquation<dim>::run ()
    {
      int de;
      std::vector<unsigned int> repetitions;
      repetitions.push_back (N - 1);
      repetitions.push_back (1);
      GridGenerator::subdivided_hyper_rectangle (triangulation, repetitions,
						 Point<2> (0.0, 0.0),
						 Point<2> (L, L / N), true);
      std::cout << "Number of active cells: " << triangulation.n_active_cells ()
	  << std::endl;

      // convert yita_middle_1D to yita_full_2D;
      for (int i = 1; i < N - 1; i++)
	yita_full_1D[i] = yita_middle_1D[i - 1];
      yita_full_1D[0] = 0.;
      yita_full_1D[N] = 0.;
      for (int i = 2; i < N; i++)
	{
	  yita_full_2D[2 * i] = yita_full_1D[i];
	  yita_full_2D[2 * i + 1] = yita_full_1D[i];
	}
      yita_full_2D[0] = yita_full_1D[0];
      yita_full_2D[2] = yita_full_1D[0];
      yita_full_2D[1] = yita_full_1D[1];
      yita_full_2D[3] = yita_full_1D[1];

      printf (">>>>>>>>>>>>>>>>>>>>>>>>>>>..1\n");
//      std::ofstream out ("grid-1.vtk");
//      GridOut grid_out;
//      grid_out.write_vtk (triangulation, out);
//      std::cout << "Grid written to grid-1.vtk" << std::endl;

      setup_system ();
      printf (">>>>>>>>>>>>>>>>>>>>>>>>>>>..5\n");
      VectorTools::interpolate (dof_handler, Initial_condition<dim> (),
				old_solution);
      solution = old_solution;
      output_results ();

      printf (">>>>>>>>>>>>>>>>>>>>>>>>>>>..2\n");
      std::vector<int> solution_table (N);
      for (int i = 0; i < N; i++)
	{
	  if (i == 0)
	    solution_table[i] = 0;
	  else if (i == 1)
	    solution_table[i] = 1;
	  else
	    solution_table[i] = i * 2;
	}

      for (int i = 2; i < N; i++)
	solution_store[i][0] = 1.;

      int period = 1;
      for (timestep_number = 1; timestep_number < total_time_step;
	  timestep_number++)
	{
	  time += time_step;

	  std::cout << "Time step " << timestep_number << " at t=" << time
	      << std::endl;

	  system_matrix.copy_from (mass_matrix);
	  system_matrix.add (time_step, laplace_matrix);
	  tmp.copy_from (mass_matrix);

//	  printf ("tmp is a %d by %d mat \n", tmp.m (), tmp.m ());

//	  std::ofstream out("haha.txt");
//	  tmp.print (out);

	  for (unsigned int i = 0; i < tmp.m (); i++)
	    {
	      SparseMatrix<double>::iterator begin = tmp.begin (i), end =
		  tmp.end (i);
	      for (; begin != end; ++begin)
		{
		  begin->value () *= yita_full_2D[i];
		}
	    }

	  system_matrix.add (time_step, tmp);
	  mass_matrix.vmult (system_rhs, old_solution);
	  std::map<types::global_dof_index, double> boundary_values;
	  VectorTools::interpolate_boundary_values (dof_handler, 0,
						    ZeroFunction<2> (),
						    boundary_values);
	  MatrixTools::apply_boundary_values (boundary_values, system_matrix,
					      solution, system_rhs);
	  VectorTools::interpolate_boundary_values (dof_handler, 1,
						    ConstantFunction<2> (0.),
						    boundary_values);
	  MatrixTools::apply_boundary_values (boundary_values, system_matrix,
					      solution, system_rhs);
	  solve_time_step ();
	  output_results ();
	  old_solution = solution;
	  solution_store[0][timestep_number] = time;
	  for (int i = 0; i < N; i++)
	    {
	      solution_store[i + 1][timestep_number] =
		  solution[solution_table[i]];
	    }
	}
      printf (">>>>>>>>>>>>>>>>>>>>>>>>>>>..3\n");

      // write solution;
      FILE * fp;
      fp = fopen ("solution_store.txt", "w+");
      for (int i = 0; i < N + 1; i++)
	{
	  for (int j = 0; j < total_time_step; j++)
	    fprintf (fp, "%2.15f,", solution_store[i][j]);
	  fprintf (fp, "\n");
	}

      fclose (fp);

      //   integrate for f0
      //   TODO:change this to better integration method

      printf ("f0: \n");
      for (int i = 0; i < N; i++)
	{
	  f0[i] = 0.0;
	  for (int j = 0; j < total_time_step - 1; j++)
	    {
	      double value_left = solution_store[i + 1][j]
		  * solution_store[i + 1][total_time_step - j - 1];
	      double value_right = solution_store[i + 1][j + 1]
		  * solution_store[i + 1][total_time_step - j - 1 - 1];
	      f0[i] = f0[i] + 0.5 * time_step * (value_left + value_right);
	    }
//	  printf ("%0.16f \n", f0[i]);
	}

      printf (">>>>>>>>>>>>>>>>>>>>>>>>>>>..6\n");

      Matfree (solution_store);
      for (int i = 1; i < N - 1; i++)
	{
	  out[i] = f0_given[i] - f0[i];
	}
      printf (">>>>>>>>>>>>>>>>>>>>>>>>>>>..7\n");
      return out;
    }

}

void
get_f0_given (double tau, double L, int N)
{
  double x[N];
  for (int i = 0; i < N; i++)
    {
      f0_given[i] = 1.;
      x[i] = i * L / (N - 1);
    }
  int x_left = ceil (N * tau / L);
  for (int i = 0; i < x_left; i++)
    {
      f0_given[i] = pow ((exp (4 * tau * x[i] / (tau * tau - x[i] * x[i])) - 1),
			 2)
	  / pow (((exp (4 * tau * x[i] / (tau * tau - x[i] * x[i])) + 1)), 2);
      f0_given[N - i - 1] = f0_given[i];
    }
}

void
myfun (int n, double * x, double * xnew)
{
  xnew[1] = x[1] * 0.5 - 2.;
  xnew[2] = x[2] * 0.5 + 3.;
//  printf ("x[1]=%.14f,x[2]=%.14f, xnew[1]=%.14f xnew[2]=%.14f \n", x[1], x[2],
//	  xnew[1], xnew[2]);
}

void
SCFT_wrapper (int N, double * in, double * out)
{
  double L = 3.72374;

  Step26::HeatEquation<2> heat_equation_solver (N, 2048, L, in);

  double* res = heat_equation_solver.run ();
  for (int i = 1; i < N - 1; i++)
    out[i] = res[i];

}

int
main ()
{
  try
    {
      using namespace dealii;
      using namespace Step26;
      int N = 33;
      int de;
      double* yita_1D = (double*) malloc (N * sizeof(double)); // this is the yita for 1D, length=N;
      double* yita_2D = (double*) malloc (N * sizeof(double) * 2); // need to convert it to 2D because mat is 2N by 2N;
      double yita_middle[N - 2]; // initial guess, the ends are bounded   // this is the middle of yita_1D, because the boundary are fixed.
      f0_given = (double*) malloc (N * sizeof(double)); // this is the ideal f0;
      double tau = 0.5302, L = 3.72374; // tau is for calculating f0_given, L is the length.

      for (int i = 0; i < N; i++)
	yita_1D[i] = i;
      // Convert yita_given to yita, because we are at a 2D problem, every node needs a associate yita value.

      get_f0_given (tau, L, N);
      printf ("f0_given\n");
      for (int i = 0; i < N; i++)
	printf ("%f \n", f0_given[i]);
      printf ("\n");

      // read data from file:
      FILE *file;
      file = fopen ("Exp_m32_n2048_IE.res", "r");
      if (file == NULL)
	{
	  fprintf (stderr, "Can't open input file in.list!\n");
	  return 1;
	}
      char buff[255];
      int line = 0;
      double c, e;
      for (int i = 0; i < 10; i++)
	fgets (buff, 255, (FILE*) file);
      while (fgets (buff, 255, (FILE*) file))
	{
	  sscanf (buff, "%lf %lf %lf %lf %lf", &c, &c, &e, &c, &c);
	  yita_middle[line] = e;
	  line++;
	  if (line == N - 2)
	    break;
	}
      fclose (file);

      int check = 1;
      qt = dmatrix (1, N - 2, 1, N - 2);
      r = dmatrix (1, N - 2, 1, N - 2);
      d = dvector (1, N - 2);
      jc = 0;
      err = 0.00000001;
      double* x_nr = dvector (1, N - 2);
      for (int i = 1; i < N - 1; i++)
	x_nr[i] = yita_middle[i - 1];

      double* out = (double*) malloc (sizeof(double) * (N - 1));
//      SCFT_wrapper (N, yita_middle, out);
//      SCFT_wrapper (N, yita_middle, out);

      for (int i = 0; i < 4; i++)
	{
	  Step26::HeatEquation<2> heat_equation_solver (N, 2048, L, yita_middle);
	   heat_equation_solver.run ();
	}

//      broydn (x_nr, N - 2, &check, SCFT_wrapper);

//      broydn (x_nr, 2, &check, myfun);

//      void broydn(float x[], int n, int *check,
//          void (*vecfunc)(int, float [], float []));

//      for (int i = 1; i < 5; i++)
//	printf (">>>>%f \n", x_nr[i]);
//      scanf ("%d", &de);

//      HeatEquation<2> heat_equation_solver (N, 2048, L, yita_middle);

//      double* out = heat_equation_solver.run ();

      for (int i = 1; i < N - 1; i++)
	printf ("out[%d]=%0.16f \n", i, out[i]);

      free (yita_2D);
      free (yita_1D);
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
	  << "----------------------------------------------------"
	  << std::endl;
      std::cerr << "Exception on processing: " << std::endl << exc.what ()
	  << std::endl << "Aborting!" << std::endl
	  << "----------------------------------------------------"
	  << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
	  << "----------------------------------------------------"
	  << std::endl;
      std::cerr << "Unknown exception!" << std::endl << "Aborting!" << std::endl
	  << "----------------------------------------------------"
	  << std::endl;
      return 1;
    }

  return 0;
}

