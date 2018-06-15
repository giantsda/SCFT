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
      HeatEquation (int N, int total_time_step);
      void
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
      if (p[0] == 0. || p[0] == 1.)
	return 0.;
      else
	return 1.;
    }

  template<int dim>
    HeatEquation<dim>::HeatEquation (int N, int total_time_step) :
	fe (1), dof_handler (triangulation), time (0.0), time_step (1. / 2047), timestep_number (
	    0), N (N), total_time_step (total_time_step)
    {
      solution_store = Matcreate (N + 1, total_time_step + 1);
      f0 = (double*) malloc (sizeof(double) * N);
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
      SolverControl solver_control (8000, 1e-6);
      SolverCG<> solver (solver_control);
      solver.solve (system_matrix, solution, system_rhs,
		    PreconditionIdentity ());

      std::cout << "     " << solver_control.last_step () << " CG iterations."
	  << std::endl;
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
    void
    HeatEquation<dim>::run ()
    {
      std::vector<unsigned int> repetitions;
      repetitions.push_back (N - 1);

      repetitions.push_back (1);
      GridGenerator::subdivided_hyper_rectangle (triangulation, repetitions,
						 Point<2> (0.0, 0.0),
						 Point<2> (1, 1./N), true);

      std::cout << "Number of active cells: " << triangulation.n_active_cells ()
	  << std::endl;

//      std::ofstream out ("grid-1.vtk");
//      GridOut grid_out;
//      grid_out.write_vtk (triangulation, out);
//      std::cout << "Grid written to grid-1.vtk" << std::endl;

      int de;

      setup_system ();

      Vector<double> yita;

      yita.reinit (solution.size ());
      yita = 1.;


      VectorTools::interpolate (dof_handler, Initial_condition<dim> (),
     				old_solution);


      solution = old_solution;

      output_results ();

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

//      for (int i = 0; i < N; i++)
//	printf ("%d>>%d  \n", i, solution_table[i]);

//      scanf ("%d", &de);

      for (int i = 2; i < N; i++)
	solution_store[i][0] = 1.;


      int period=1;


      for (timestep_number=1; timestep_number < total_time_step; timestep_number++)
	{
	  time += time_step;

	  std::cout << "Time step " << timestep_number << " at t=" << time
	      << std::endl;

	  system_matrix.copy_from (mass_matrix);
	  system_matrix.add (time_step, laplace_matrix);
	  tmp.copy_from (mass_matrix);

//	  std::ofstream out("haha.txt");
//	  tmp.print (out);
	  for (unsigned int i = 0; i < tmp.m (); i++)
	    {
	      SparseMatrix<double>::iterator begin = tmp.begin (i), end =
		  tmp.end (i);
	      for (; begin != end; ++begin)
		{
		  begin->value () *= yita[i];
		}
	    }

//	  tmp.print (std::cout);
//	  old_solution.print (std::cout);
//	  scanf("%d",&de);
	  system_matrix.add (time_step, tmp);

	  mass_matrix.vmult (system_rhs, old_solution);

//
//	  if (period == 1)
//	    {
//	      std::vector<
//		  GridTools::PeriodicFacePair<
//		      typename Triangulation<2>::cell_iterator> > periodicity_vector;
//	      GridTools::collect_periodic_faces (triangulation, 2, 3, 1,
//						 periodicity_vector);
//	      triangulation.add_periodicity (periodicity_vector);
//	      period = 0;
//	    }

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

//	  scanf("%d",&de);

	  old_solution = solution;

	  solution_store[0][timestep_number] = time;
	  for (int i = 0; i < N; i++)
	    {
	      solution_store[i + 1][timestep_number] =
		  solution[solution_table[i]];
	    }

	}

      FILE * fp;

      fp = fopen ("fhaha.txt", "w+");

      // plot solution;
      for (int i = 0; i < N + 1; i++)
	{
	  for (int j = 0; j < total_time_step; j++)
	    fprintf (fp, "%f,", solution_store[i][j]);
	  fprintf (fp, "\n");
	}

      fclose (fp);

      //   integrate for f0
      //   TODO:change this to better integration method

      printf ("f0: \n");
      for (int i = 0; i < N; i++)
	{
	  f0[i] = 0.0;
	  for (int j = 0; j < total_time_step; j++)
	    {
	      double value_left = solution_store[i + 1][j]
		  * solution_store[i + 1][total_time_step - j];
	      double value_right = solution_store[i + 1][j + 1]
		  * solution_store[i + 1][total_time_step - j - 1];
	      f0[i] = f0[i] + 0.5 * time_step * (value_left + value_right);
	    }
	  printf ("%f \n", f0[i]);
	}

      Matfree (solution_store);
    }

}

int
main ()
{
  try
    {
      using namespace dealii;
      using namespace Step26;

      HeatEquation<2> heat_equation_solver (33, 2048);
//      HeatEquation<2> heat_equation_solver (33, 10);
      heat_equation_solver.run ();

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
