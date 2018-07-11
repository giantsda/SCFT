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
#include <utility>
#include "nr.h"
#include "nrutil.h"

double* f0_given;
double **qt, **r, *d, /* qt[1:n,1:n], r[1:n,1:n] and d[1:n] must be allocated in the calling program */
err; /* Passed in as the convergence criterion, and returned with the actual residual error */
int funcerr, /* Flag for error in evaluating the function in Broyden method */
jc; /* For re-use of Jacobian. jc denotes the status of Jacobian calculation: 0 for not calculated,
 1 for previously calculated, 2 for currently calculated */
int PRINT, local_iteration = 0;

int de; // My debug varaibel

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
      HeatEquation ();

      HeatEquation (HeatEquation&& other);
      HeatEquation&
      operator = (HeatEquation &other);
      void
      print_point ();
      ~HeatEquation () = default;
      int
      get_refine_times ();
      void
      get_x (std::vector<double> & in);
      void
      set_refine_times (int a);
      int
      get_N ();
      void
      set_yita_middle_1D (double * in);
      double *
      run ();
      double *
      run_experiemnt ();
      void
      refine_mesh ();
      void
      update_internal_data ();
    private:
      void
      setup_system ();
      void
      solve_time_step ();
      void
      output_results () const;
      void
      build_solution_table ();
      void
      output_results_for_yita_full_2D_t () const;
      void
      print_and_save_yita_1D ();

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
      int refine_times;
      std::map<int, int> solution_table_1D_to_2D;
      std::map<int, int> solution_table_2D_to_1D;
      std::map<double, int> solution_table_x_to_2D;
      std::map<int, double> solution_table_2D_to_x;
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
	refine_times (0), fe (1), dof_handler (triangulation), time (0.0), time_step (
	    1. / 2047), timestep_number (0), N (N), total_time_step (
	    total_time_step), L (L), yita_middle_1D (yita)
    {
      time_step = 1. / (total_time_step - 1);
      solution_store = Matcreate (N + 1, total_time_step + 1);
      f0 = (double*) malloc (sizeof(double) * N);
      yita_full_1D = (double*) malloc (sizeof(double) * N);
      yita_full_2D = (double*) malloc (sizeof(double) * 2 * N);
      out = (double*) malloc (sizeof(double) * (N - 1));
    }

  template<int dim>
    HeatEquation<dim>::HeatEquation () :
	fe (1), dof_handler (triangulation)
    {
    }

  template<int dim>
    HeatEquation<dim>&
    HeatEquation<dim>::operator = (HeatEquation & other)

    {
      refine_times = other.refine_times;
      time = other.time;
      time_step = 1. / 2047;
      timestep_number = other.timestep_number;
      N = other.N;
      total_time_step = other.total_time_step;
      L = other.L;
      yita_middle_1D = other.yita_middle_1D;
      time_step = 1. / (total_time_step - 1);
      solution_store = Matcreate (N + 1, total_time_step + 1);
      f0 = other.f0;
      yita_full_1D = other.yita_full_1D;
      yita_full_2D = other.yita_full_2D;
      out = other.out;

      other.solution_store = nullptr;
      other.yita_full_1D = nullptr;
      other.yita_full_2D = nullptr;
      other.out = nullptr;

    }

  template<int dim>
    int
    HeatEquation<dim>::get_refine_times ()
    {
      return refine_times;
    }

  template<int dim>
    void
    HeatEquation<dim>::get_x (std::vector<double> & in)
    {
      typename Triangulation<dim>::active_vertex_iterator vertex =
	  triangulation.begin_active_vertex (), endv =
	  triangulation.end_vertex ();
      int i = 0;
      Point<2> p;
      for (; vertex != endv; vertex++)
	{
	  p = vertex->vertex (1);
	  if (p (1) == 0)
	    {
	      in[i] = p (0);
	      i++;
	    }
	}
      std::sort (in.begin (), in.end ());
    }

  template<int dim>
    void
    HeatEquation<dim>::set_refine_times (int a)
    {
      refine_times = a;
    }

  template<int dim>
    int
    HeatEquation<dim>::get_N ()
    {
      return N;
    }

  template<int dim>
    void
    HeatEquation<dim>::set_yita_middle_1D (double * in)
    {
      yita_middle_1D = in;
    }

  template<int dim>
    HeatEquation<dim>::HeatEquation (HeatEquation&& other) :
	refine_times (other.refine_times), fe (other.fe), dof_handler (
	    other.dof_handler), time (other.time), time_step (1. / 2047), timestep_number (
	    other.timestep_number), N (other.N), total_time_step (
	    other.total_time_step), L (other.L), yita_middle_1D (
	    other.yita_middle_1D)
    {
      time_step = 1. / (total_time_step - 1);
      solution_store = Matcreate (N + 1, total_time_step + 1);
      f0 = other.f0;
      yita_full_1D = other.yita_full_1D;
      yita_full_2D = other.yita_full_2D;
      out = other.out;

      other.solution_store = nullptr;
      other.yita_full_1D = nullptr;
      other.yita_full_2D = nullptr;
      other.out = nullptr;

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
      N = triangulation.n_active_cells () + 1;

    }

  template<int dim>
    void
    HeatEquation<dim>::solve_time_step ()
    {
      SolverControl solver_control (80000, 1e-13);
      SolverCG<> solver (solver_control);
      solver.solve (system_matrix, solution, system_rhs,
		    PreconditionIdentity ());
//      std::cout << "     " << solver_control.last_step () << " CG iterations."
//	  << std::endl;
    }
  template<int dim>
    void
    HeatEquation<dim>::refine_mesh ()
    {

//      output_results_for_yita_full_2D_t ();
      print_and_save_yita_1D ();

      scanf ("%d", &de);
#undef float
      Vector<float> estimated_error_per_cell (triangulation.n_active_cells ());
      Vector<double> yita_full_2D_t (N * 2);

      for (int i = 0; i < 2 * N; i++)
	yita_full_2D_t[i] = yita_full_2D[i];

      KellyErrorEstimator<dim>::estimate (dof_handler,
					  QGauss<dim - 1> (fe.degree + 1),
					  typename FunctionMap<dim>::type (),
					  yita_full_2D_t,
					  estimated_error_per_cell);

#define float double

//      GridRefinement::refine_and_coarsen_fixed_fraction (
//	  triangulation, estimated_error_per_cell, 1., 0.0);

      GridRefinement::refine (triangulation, estimated_error_per_cell, 0.50,
			      200);

      typename DoFHandler<dim>::active_cell_iterator cell =
	  dof_handler.begin_active (), endc = dof_handler.end ();

      for (cell = dof_handler.begin_active (); cell != endc; ++cell)
	{
	  if (cell->refine_flag_set ())
	    cell->set_refine_flag (RefinementCase<dim>::cut_axis (0));
	}

      SolutionTransfer<dim> solution_trans (dof_handler);
      Vector<double> previous_solution, new_solution;
      previous_solution = solution;
      triangulation.prepare_coarsening_and_refinement ();
      solution_trans.prepare_for_coarsening_and_refinement (previous_solution);
      triangulation.execute_coarsening_and_refinement ();
      setup_system ();
      new_solution.reinit (dof_handler.n_dofs ());
      solution_trans.interpolate (previous_solution, new_solution);

      new_solution.print (std::cout);
      //	scanf("%d", &de);
      refine_times++;

    }

  template<int dim>
    void
    HeatEquation<dim>::print_point ()
    {
      printf ("The solution_store is %x\n", solution_store);
      printf ("The solution_store[2][0] is %x\n", &solution_store[2][0]);

//      scanf ("%d", &de);
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
      printf ("%s is written \n", filename.c_str ());
    }

  template<int dim>
    void
    HeatEquation<dim>::build_solution_table ()
    {

      Point<2, double> P;
      std::vector<types::global_dof_index> loc_dof_indices (fe.dofs_per_cell);
      typename DoFHandler<2>::active_cell_iterator cell =
	  dof_handler.begin_active (), endc = dof_handler.end ();
      for (cell = dof_handler.begin_active (); cell != endc; cell++)
	{
	  cell->get_dof_indices (loc_dof_indices);
	  for (int i = 0; i < fe.dofs_per_cell; i++)
	    {
	      P = cell->vertex (i);
	      if (P (1) == 0)
		{
//		  std::cout << P (0) << "=" << loc_dof_indices[i] << std::endl;
		  solution_table_x_to_2D.insert (
		      std::pair<double, int> (P (0), loc_dof_indices[i]));
		}
	    }
	}
      int ii = 0;
      for (auto itr = solution_table_x_to_2D.begin ();
	  itr != solution_table_x_to_2D.end (); ++itr)
	{
	  solution_table_1D_to_2D.insert (
	      std::pair<int, int> (ii, itr->second));
	  ii++;
	}
// Now build table_x_to_1D

      std::map<double, int> solution_table_x_to_1D;
      std::vector<double> x;
      x.resize (get_N ());
      get_x (x);
      for (int i = 0; i < x.size (); i++)
	solution_table_x_to_1D.insert (std::pair<double, int> (x[i], i));

      for (cell = dof_handler.begin_active (); cell != endc; cell++)
	{
	  cell->get_dof_indices (loc_dof_indices);
	  for (int i = 0; i < fe.dofs_per_cell; i++)
	    {
	      P = cell->vertex (i);
	      solution_table_2D_to_x.insert (
		  std::pair<int, double> (loc_dof_indices[i], P (0)));
	    }
	}

      for (int i = 0; i < 2 * N; i++)
	{
	  double x = solution_table_2D_to_x.find (i)->second;
	  solution_table_2D_to_1D.insert (
	      std::pair<int, int> (i, solution_table_x_to_1D.find (x)->second));
	}

//      printf("print Table: \n");

//      printf ("N=%d \n", get_N ());
//      printf ("solution_table_x_to_2D\n");
//      for (auto itr = solution_table_x_to_2D.begin ();
//	  itr != solution_table_x_to_2D.end (); ++itr)
//	{
//	  std::cout << '\t' << itr->first << "---" << itr->second << '\n';
//	}
//      std::cout << std::endl;
//
//      printf ("solution_table_2D_to_x\n");
//      for (auto itr = solution_table_2D_to_x.begin ();
//	  itr != solution_table_2D_to_x.end (); ++itr)
//	{
//	  std::cout << '\t' << itr->first << "---" << itr->second << '\n';
//	}
//      std::cout << std::endl;
//
//      printf ("solution_table_1D_to_2D\n");
//      for (auto itr = solution_table_1D_to_2D.begin ();
//	  itr != solution_table_1D_to_2D.end (); ++itr)
//	{
//	  std::cout << '\t' << itr->first << "---" << itr->second << '\n';
//	}
//      std::cout << std::endl;
//
//      printf ("solution_table_2D_to_1D\n");
//      for (auto itr = solution_table_2D_to_1D.begin ();
//	  itr != solution_table_2D_to_1D.end (); ++itr)
//	{
//	  std::cout << '\t' << itr->first << "---" << itr->second << '\n';
//	}
//      std::cout << std::endl;
//
//      for (cell = dof_handler.begin_active (); cell != endc; cell++)
//	{
//	  cell->get_dof_indices (loc_dof_indices);
//	  for (int i = 0; i < fe.dofs_per_cell; i++)
//	    {
//	      P = cell->vertex (i);
//	      std::cout << P << "->" << loc_dof_indices[i] << std::endl;
//	    }
//	}
//
//      scanf ("%d", &de);

    }

  template<int dim>
    void
    HeatEquation<dim>::output_results_for_yita_full_2D_t () const
    {
      Vector<double> yita_full_2D_t (N * 2);

      for (int i = 0; i < 2 * N; i++)
	{
	  yita_full_2D_t[i] = yita_full_2D[i];
	}
      DataOut<dim> data_out;

      data_out.attach_dof_handler (dof_handler);
      data_out.add_data_vector (yita_full_2D_t, "U");

      data_out.build_patches ();

      const std::string filename = "yita_full_2D-"
	  + Utilities::int_to_string (refine_times, 3) + ".vtk";
      std::ofstream output (filename.c_str ());
      data_out.write_vtk (output);
    }

  template<int dim>
    void
    HeatEquation<dim>::print_and_save_yita_1D ()
    {

      printf ("Solved ! \n middle_solution: N=%d \n", get_N ());

      // print and write  solution;

      FILE * fp;
      char filename[64], num[64];
      sprintf (num, "%d", get_N ());

      strcpy (filename, "solution_yita_1D_N= ");
      strcat (filename, num);
      strcat (filename, ".txt");

      fp = fopen (filename, "w+");
      fprintf (fp, "N= %d \n", get_N ());
      for (int i = 0; i < N; i++)
	{
	  printf ("solution[%d]=%2.15f \n", i, yita_full_1D[i]);
	  fprintf (fp, "%d, %2.15f\n", i, yita_full_1D[i]);
	}

      fclose (fp);

      printf ("%s is written. \n", filename);
//      scanf ("%d", &de);

    }

  template<int dim>
    void
    HeatEquation<dim>::update_internal_data ()
    {
      time = 0.;
      timestep_number = 0;
//      N = triangulation.n_active_cells () + 1;
      printf ("N=%d  \n", N);
      Matfree (solution_store);
      solution_store = Matcreate (N + 1, total_time_step + 1);
      f0 = (double*) realloc (f0, sizeof(double) * N);
      yita_full_1D = (double*) realloc (yita_full_1D, sizeof(double) * N);
      yita_full_2D = (double*) realloc (yita_full_2D, sizeof(double) * 2 * N);
      out = (double*) realloc (out, sizeof(double) * (N - 1));

      free_dmatrix (qt, 1, N - 2, 1, N - 2);
      qt = dmatrix (1, N - 2, 1, N - 2);
      free_dmatrix (r, 1, N - 2, 1, N - 2);
      r = dmatrix (1, N - 2, 1, N - 2);
      free_dvector (d, 1, N - 2);
      d = dvector (1, N - 2);
      jc = 0;
      err = 0.00000001;
      solution_table_1D_to_2D.clear ();
      solution_table_2D_to_1D.clear ();
      solution_table_x_to_2D.clear ();
      solution_table_2D_to_x.clear ();
    }

  template<int dim>
    double *
    HeatEquation<dim>::run ()
    {

      if (refine_times == 0)
	{
	  std::vector<unsigned int> repetitions;
	  repetitions.push_back (N - 1);
	  repetitions.push_back (1);
	  GridGenerator::subdivided_hyper_rectangle (triangulation, repetitions,
						     Point<2> (0.0, 0.0),
						     Point<2> (L, L / N), true);
	  refine_times++;
	}
      std::cout << "Number of active cells: " << triangulation.n_active_cells ()
	  << std::endl;

      //      std::ofstream out ("grid-1.vtk");
      //      GridOut grid_out;
      //      grid_out.write_vtk (triangulation, out);
      //      std::cout << "Grid written to grid-1.vtk" << std::endl;

      setup_system ();

      VectorTools::interpolate (dof_handler, Initial_condition<dim> (),
				old_solution);
      solution = old_solution;
//      output_results ();

      if (local_iteration == 0)
	build_solution_table ();

      // convert yita_middle_1D to yita_full_2D;
      for (int i = 1; i < N - 1; i++)
	yita_full_1D[i] = yita_middle_1D[i];
      yita_full_1D[0] = 0.;
      yita_full_1D[N - 1] = 0.;

      for (int i = 0; i < 2 * N; i++)
	{
	  int j = solution_table_2D_to_1D.find (i)->second;
	  yita_full_2D[i] = yita_full_1D[j];
//	  printf ("%d -> %d\n", i, j);
	}

//      for (int i = 0; i < 2 * N; i++)
//	printf ("yita_full_2D[%d]=%2.15f\n", i, yita_full_2D[i]);

      for (timestep_number = 1; timestep_number < total_time_step;
	  timestep_number++)
	{
	  time += time_step;
	  system_matrix.copy_from (mass_matrix);
	  system_matrix.add (time_step, laplace_matrix);
	  tmp.copy_from (mass_matrix);

	  //	  if (refine_times == 2)
	  //	    {
	  //	      system_matrix.print (std::cout);
	  //	      system_rhs.print (std::cout);
	  //	      scanf ("%d", &de);
	  //	    }

	  //	  std::ofstream out("haha.txt");
	  //	  tmp.print (out);

	  for (unsigned int i = 0; i < tmp.m (); i++)
	    {
	      SparseMatrix<double>::iterator begin = tmp.begin (i), end =
		  tmp.end (i);
	      for (; begin != end; ++begin)
		{
		  //		  printf ("yita_full_2D[%d]=%f \n", i, yita_full_2D[i]);
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

//	  if (refine_times == 2 && local_iteration >= 41)
//	    {
//	      system_matrix.print (std::cout);
//	      system_rhs.print (std::cout);
//	      scanf ("%d", &de);
//	    }

	  solve_time_step ();

//	  if (refine_times == 2)
//	    {
//	  output_results ();
//	  scanf ("%d", &de);
//	    }

//	  for (auto itr = solution_table.begin (); itr != solution_table.end ();
//	      ++itr)
//	    {
//	      std::cout << '\t' << itr->first << '\t' << itr->second << '\n';
//	    }
//	  std::cout << std::endl;
//
//	  scanf ("%d", &de);

	  old_solution = solution;
	  solution_store[0][timestep_number] = time;
	  for (int i = 0; i < N; i++)
	    {
	      solution_store[i + 1][timestep_number] =
		  solution[solution_table_1D_to_2D.find (i)->second];
	    }
	}

      // write solution;
//      if (refine_times == 2)
//	{
//	  FILE * fp;
//	  fp = fopen ("solution_store.txt", "w+");
//	  for (int i = 0; i < N + 1; i++)
//	    {
//	      for (int j = 0; j < total_time_step; j++)
//		fprintf (fp, "%2.15f,", solution_store[i][j]);
//	      fprintf (fp, "\n");
//	    }
//
//	  fclose (fp);
//	  scanf ("%d", &de);
//	}

      //   integrate for f0
      //   TODO:change this to better integration method

      //      printf ("f0: \n");
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

      for (int i = 1; i < N - 1; i++)
	{
	  out[i] = f0_given[i] - f0[i];
	}
      return out;
    }

  template<int dim>
    double *
    HeatEquation<dim>::run_experiemnt ()
    {

      if (refine_times == 0)
	{
	  std::vector<unsigned int> repetitions;
	  repetitions.push_back (5);
	  repetitions.push_back (1);
	  GridGenerator::subdivided_hyper_rectangle (triangulation, repetitions,
						     Point<2> (0.0, 0.0),
						     Point<2> (1, 1), true);
	  refine_times++;
	}
      std::cout << "Number of active cells: " << triangulation.n_active_cells ()
	  << std::endl;

      std::ofstream out ("grid-1.vtk");
      GridOut grid_out;
      grid_out.write_vtk (triangulation, out);
      std::cout << "Grid written to grid-1.vtk" << std::endl;

      setup_system ();

      for (int i = 0; i < solution.size (); i++)
	solution[i] = i;

      output_results ();

      build_solution_table ();

//      scanf ("%d", &de);

      // convert yita_middle_1D to yita_full_2D;
      for (int i = 1; i < N - 1; i++)
	yita_full_1D[i] = yita_middle_1D[i];
      yita_full_1D[0] = 0.;
      yita_full_1D[N - 1] = 0.;

      for (int i = 0; i < 2 * N; i++)
	{
	  int j = solution_table_2D_to_1D.find (i)->second;
	  yita_full_2D[i] = yita_full_1D[j];
//	  printf ("%d -> %d\n", i, j);
	}

      typename DoFHandler<dim>::active_cell_iterator cell =
	  dof_handler.begin_active (), endc = dof_handler.end ();
	{
	  int i = 0;
	  for (cell = dof_handler.begin_active (); cell != endc; ++cell)
	    {
	      if (i == 0 || i == 4)
		cell->set_refine_flag (RefinementCase<dim>::cut_axis (0));
	      i++;
	    }
	}

//      triangulation.prepare_coarsening_and_refinement ();
      triangulation.execute_coarsening_and_refinement ();
      setup_system ();
      for (int i = 0; i < solution.size (); i++)
	solution[i] = i;
      update_internal_data ();
      timestep_number = 1;

      printf ("-------------------------------------------------\n");
      build_solution_table ();
      output_results ();

      double *a;
      return a;
    }

}

void
get_f0_given (double tau, double L, int N, std::vector<double> & x)
{

  for (int i = 0; i < N; i++)
    f0_given[i] = 1.;
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
}

void
myfun (int n, double * x, double * xnew)
{
  xnew[1] = x[1] * 0.5 - 2.;
  xnew[2] = x[2] * 0.5 + 3.;
}

Step26::HeatEquation<2> heat_equation_solver;

void
SCFT_wrapper (int N, double * in, double * out)
{

  N = N + 2;
//  for (int i = 1; i < N - 1; i++)
//    printf ("in[%d]=%2.15f \n", i, in[i]);
//  printf ("CONTINUE>>>>>>>>>>>>>>>\n");
  double* res = heat_equation_solver.run ();

  for (int i = 1; i < N - 1; i++)
    out[i] = res[i];

  for (int i = 1; i < N - 1; i++)
    printf (
	"in[%d]=%2.15f ; out[%d]=%2.15f ; local_iteration:%d refine_times:%d \n",
	i, in[i], i, out[i], local_iteration,
	heat_equation_solver.get_refine_times ());
  fflush (stdout);
  local_iteration++;

//  if (heat_equation_solver.get_refine_times () == 2)
//    scanf ("%d", &de);

}

int
main ()
{
  try
    {
      using namespace dealii;
      using namespace Step26;
      int N = 33;

      double* yita_middle = (double*) malloc ((N - 1) * sizeof(double) * 2); // initial guess, the ends are bounded   // this is the middle of yita_1D, because the boundary are fixed.
      f0_given = (double*) malloc (N * sizeof(double)); // this is the ideal f0;
      double tau = 0.5302, L = 3.72374; // tau is for calculating f0_given, L is the length.
      std::vector<double> x;

      // read data from file:
	{
	  FILE *file;
	  file = fopen ("N=33_for_read.txt", "r");
	  if (file == NULL)
	    {
	      fprintf (stderr, "Can't open input file in.list!\n");
	      return 1;
	    }
	  char buff[255];
	  int line = 1;
	  int i;
	  double value;

	  fgets (buff, 255, (FILE*) file);
	  sscanf (buff, "N= %d", &N);
	  while (fgets (buff, 255, (FILE*) file))
	    {
	      sscanf (buff, "%d , %lf", &i, &value);
	      yita_middle[i] = value;
	    }
	  fclose (file);
	}

//      for (int i = 0; i < N; i++)
//	printf ("yita_middle[%d]=%2.15f \n", i, yita_middle[i]);
//      scanf ("%d", &de);

      int check = 1;
      qt = dmatrix (1, N - 2, 1, N - 2);
      r = dmatrix (1, N - 2, 1, N - 2);
      d = dvector (1, N - 2);
      jc = 0;
      err = 0.00000001;
      double* x_nr = dvector (1, N - 2);
      for (int i = 1; i < N - 1; i++)
	x_nr[i] = yita_middle[i];
//	x_nr[i] = 0.;

      for (int i = 1; i < N - 1; i++)
	printf ("x_nr[%d]=%f \n", i, x_nr[i]);

      Step26::HeatEquation<2> other (N, 2048, L, x_nr); // This tells the yita_middle_1D and x_nr has the same address;

      heat_equation_solver = other;
      /*The assignment constructor transfer the address of x_nr
       to heat_equation_solver.yita_middle_1D*/

      // calculate f0_given
      x.resize (N);
      for (int i = 0; i < N; i++)
	{
	  x[i] = i * L / (N - 1);
	}
      get_f0_given (tau, L, N, x);

      //      printf ("f0_given\n");
      //      for (int i = 0; i < N; i++)
      //	printf ("%f \n", f0_given[i]);
      //      printf ("\n");

      /*--------------------------------------------------------------*/
//      for (int u = 0; u < 500; u++)
//	{
//	  broydn (x_nr, N - 2, &check, SCFT_wrapper);
//	  print_solution (x_nr, N);
//	  heat_equation_solver.refine_mesh ();
//	  heat_equation_solver.update_internal_data ();
//	  free_dvector (x_nr, 1, N - 2);
//	  N = heat_equation_solver.get_N ();
//	  x_nr = dvector (1, N - 2);
//	  for (int i = 1; i < N - 1; i++)
//	    x_nr[i] = 0.;
//	  heat_equation_solver.set_yita_middle_1D (x_nr);
//	  f0_given = (double*) realloc (f0_given, N * sizeof(double));
//	  x.resize (heat_equation_solver.get_N ());
//	  heat_equation_solver.get_x (x);
//	  for (auto u : x)
//	    printf ("%f \n", u);
//
//	  get_f0_given (tau, L, N, x);
//	  local_iteration = 0;
//	}

      broydn (x_nr, N - 2, &check, SCFT_wrapper);
      heat_equation_solver.refine_mesh ();
      heat_equation_solver.update_internal_data ();
      free_dvector (x_nr, 1, N - 2);
      N = heat_equation_solver.get_N ();
      x_nr = dvector (1, N - 2);
      for (int i = 1; i < N - 1; i++)
	x_nr[i] = 0.;
      heat_equation_solver.set_yita_middle_1D (x_nr);
      f0_given = (double*) realloc (f0_given, N * sizeof(double));
      x.resize (heat_equation_solver.get_N ());
      heat_equation_solver.get_x (x);
      get_f0_given (tau, L, N, x);
      local_iteration = 0;
      broydn (x_nr, N - 2, &check, SCFT_wrapper);
      heat_equation_solver.refine_mesh ();

//      heat_equation_solver.run_experiemnt ();

//      double haha[N];
//      SCFT_wrapper (N - 2, x_nr, haha);

      free (yita_middle);

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
