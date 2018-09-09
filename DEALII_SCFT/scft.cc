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
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>

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

double*
adm (double *x, int n, int *check, void
(*funcvmix) (int, double *x, double *xnew),
     int flag);

double **
Matcreate (int r, int c) // The elements in the rows are next to each other.
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

  void
  get_f0_given (double tau, int N, std::vector<double> & x);

  template<int dim>
    class HeatEquation
    {
    public:

      HeatEquation (double tau, int N, int total_time_step, double L,
		    double* yita);
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
      get_N () const;
      void
      set_yita_middle_1D (double * in);
      double *
      run ();
      double *
      run_experiemnt ();
      void
      refine_mesh (std::vector<double> & in);
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
      output_mesh ();

      void
      print_and_save_yita_1D ();

      Triangulation<dim> triangulation;
      FE_Q<dim> fe;
      DoFHandler<dim> dof_handler;

      ConstraintMatrix constraints;

      SparsityPattern sparsity_pattern;
      SparseMatrix<double> A; // A: mass matrix (fi[i],fi[j])
      SparseMatrix<double> B; // B: laplace_matrix (d(fi[i]),d(fi[j]))
      SparseMatrix<double> system_matrix;
      SparseMatrix<double> C; // C: (yita[i]*fi[i],fi[j])

      Vector<double> Xnp1;
      Vector<double> Xn;
      Vector<double> X_internal_1;
      Vector<double> X_internal_2;
      Vector<double> X_internal_3;
      Vector<double> system_rhs;

      double tau;
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
      Assert(component == 0, ExcIndexRange (component, 0, 1));
      Assert(dim == 2, ExcNotImplemented ());
      if (p[0] == 0. || p[0] == 3.72374)
	return 0.;
      else
	return 1.;
    }

  template<int dim>
    HeatEquation<dim>::HeatEquation (double tau, int N, int total_time_step,
				     double L, double* yita) :
	tau (tau), fe (1), dof_handler (triangulation), time (0.0), time_step (
	    1. / (total_time_step - 1)), timestep_number (0), N (N), total_time_step (
	    total_time_step), L (L), yita_middle_1D (yita), refine_times (0)
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
      tau = other.tau;
      refine_times = other.refine_times;
      time = other.time;
      time_step = other.time_step;
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
// old version, for better performance. make sure, all vertex is loop over once or will get a seg error.
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

//      typename Triangulation<dim>::active_vertex_iterator vertex =
//	  triangulation.begin_active_vertex (), endv =
//	  triangulation.end_vertex ();
//      Point<2> P;
//      std::map<double, int> map;
//      for (vertex = triangulation.begin_active_vertex (); vertex != endv;
//	  vertex++)
//	{
//	  P = vertex->vertex (1);
//
//	  if (P (1) == 0)
//	    {
//	      map.insert (std::pair<double, int> (P (0), 0));
//	    }
//	}
//
//      if (map.size () != in.size ())
//	{
//	  printf ("map.size () != in.size () in get_x()");
//	  exit (EXIT_FAILURE);
//	}
//      int i = 0;
//      for (auto itr = map.begin (); itr != map.end (); ++itr, i++)
//	{
//	  in[i] = itr->first;
//	}

    }

  template<int dim>
    void
    HeatEquation<dim>::set_refine_times (int a)
    {
      refine_times = a;
    }

  template<int dim>
    int
    HeatEquation<dim>::get_N () const
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
	tau (other.tau), refine_times (other.refine_times), fe (other.fe), dof_handler (
	    other.dof_handler), time (other.time), time_step (
	    1. / other.time_step), timestep_number (other.timestep_number), N (
	    other.N), total_time_step (other.total_time_step), L (other.L), yita_middle_1D (
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

      A.reinit (sparsity_pattern);

      B.reinit (sparsity_pattern);
      system_matrix.reinit (sparsity_pattern);
      C.reinit (sparsity_pattern);

      MatrixCreator::create_mass_matrix (dof_handler,
					 QGauss<dim> (fe.degree + 1), A);
      MatrixCreator::create_laplace_matrix (dof_handler,
					    QGauss<dim> (fe.degree + 1), B);

      Xnp1.reinit (dof_handler.n_dofs ());
      Xn.reinit (dof_handler.n_dofs ());
      system_rhs.reinit (dof_handler.n_dofs ());
      N = triangulation.n_active_cells () + 1;
      f0_given = (double*) realloc (f0_given, N * sizeof(double));
      std::vector<double> x;
      x.resize (N);
      get_x (x);
      get_f0_given (tau, N, x);
      update_internal_data ();
      build_solution_table ();
    }

  template<int dim>
    void
    HeatEquation<dim>::solve_time_step ()
    {
      SolverControl solver_control (80000, 1e-17);
      SolverCG<> solver (solver_control);
      solver.solve (system_matrix, Xnp1, system_rhs, PreconditionIdentity ());
//      std::cout << "     " << solver_control.last_step () << " CG iterations."
//	  << std::endl;
    }
  template<int dim>
    void
    HeatEquation<dim>::refine_mesh (std::vector<double> & in)
    {
      print_and_save_yita_1D ();
      output_results_for_yita_full_2D_t ();
      output_mesh ();
//      scanf ("%d", &de);

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
      previous_solution = yita_full_2D_t;
//      triangulation.prepare_coarsening_and_refinement ();
      solution_trans.prepare_for_coarsening_and_refinement (previous_solution);
      triangulation.execute_coarsening_and_refinement ();
      setup_system ();
      new_solution.reinit (dof_handler.n_dofs ());
      solution_trans.interpolate (previous_solution, new_solution);

      in.resize (get_N ());
      for (unsigned int i = 0; i < in.size (); i++)
	in[i] = new_solution[solution_table_1D_to_2D.find (i)->second];

      refine_times++;

//      for (int i = 0; i < in.size (); i++)
//	printf ("in[%d]=%0.5f \n", i, in[i]);
//      scanf ("%d", &de);

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
      data_out.add_data_vector (Xnp1, "U");

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
	  for (unsigned int i = 0; i < fe.dofs_per_cell; i++)
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
      for (unsigned int i = 0; i < x.size (); i++)
	{
	  solution_table_x_to_1D.insert (std::pair<double, int> (x[i], i));
	}

      for (cell = dof_handler.begin_active (); cell != endc; cell++)
	{
	  cell->get_dof_indices (loc_dof_indices);
	  for (unsigned int i = 0; i < fe.dofs_per_cell; i++)
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

//      printf ("print Table: \n");
//
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

      const std::string filename = "yita_full_2D_N="
	  + Utilities::int_to_string (get_N (), 3) + ".vtk";
      std::ofstream output (filename.c_str ());
      data_out.write_vtk (output);
      printf ("%s is written. \n", filename.c_str ());
    }

  template<int dim>
    void
    HeatEquation<dim>::output_mesh ()
    {

      const std::string filename = "yita_full_2D_N="
	  + Utilities::int_to_string (get_N (), 3) + ".msh";
      std::ofstream output (filename.c_str ());

      GridOutFlags::Msh msh_flags (true, true);

      GridOut gridOut;
      gridOut.set_flags (msh_flags);
      gridOut.write_msh (triangulation, output);
      printf ("%s is written. \n", filename.c_str ());
    }

  template<int dim>
    void
    HeatEquation<dim>::print_and_save_yita_1D ()
    {

      printf ("Solved ! \n middle_solution: N=%d \n", get_N ());

      // print and write  solution;
      std::vector<double> x;
      x.resize (N);
      get_x (x);
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
	  fprintf (fp, "%d,%2.15f,%2.15f\n", i, x[i], yita_full_1D[i]);
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

  void
  get_f0_given (double tau, int N, std::vector<double> & x)
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

  template<int dim>
    double *
    HeatEquation<dim>::run ()
    {

      if (refine_times == 0)
	{
	  if (1)
	    {
	      std::vector<unsigned int> repetitions;
	      repetitions.push_back (N - 1);
	      repetitions.push_back (1);
	      GridGenerator::subdivided_hyper_rectangle (triangulation,
							 repetitions,
							 Point<2> (0.0, 0.0),
							 Point<2> (L, L / N),
							 true);
	    }
	  else
	    {
	      GridIn<dim> grid_in;
	      grid_in.attach_triangulation (triangulation);
	      std::ifstream input_file ("yita_full_2D_N=043.msh");
	      Assert(dim == 2, ExcInternalError ());
	      grid_in.read_msh (input_file);
	    }
	  setup_system (); // The first time, the triangulation is generated and system is set up. The
	  // Next time, it is setup in the refine();
	  refine_times++;
	}
      std::cout << "Number of active cells: " << triangulation.n_active_cells ()
	  << std::endl;

      VectorTools::interpolate (dof_handler, Initial_condition<dim> (), Xn);
      Xnp1 = Xn;

      // convert yita_middle_1D to yita_full_2D;
      for (int i = 1; i < N - 1; i++)
	yita_full_1D[i] = yita_middle_1D[i];
      yita_full_1D[0] = 0.;
      yita_full_1D[N - 1] = 0.;

      for (int i = 0; i < 2 * N; i++)
	{
	  int j = solution_table_2D_to_1D.find (i)->second;
	  yita_full_2D[i] = yita_full_1D[j];
	}

//      for (int i = 0; i < 2 * N; i++)
//	printf ("yita_full_2D[%d]=%2.15f \n", i, yita_full_2D[i]);
//      scanf ("%d", &d);

      system_matrix.copy_from (A);
      system_matrix.add (time_step, B);
      C.copy_from (A);
      for (unsigned int i = 0; i < C.m (); i++)
	{
	  SparseMatrix<double>::iterator begin = C.begin (i), end = C.end (i);
	  for (; begin != end; ++begin)
	    {
	      begin->value () *= yita_full_2D[i];
	    }
	}
      system_matrix.add (time_step, C);

      for (int i = 2; i < N; i++)
	solution_store[i][0] = 1.;
      solution_store[0][0] = 0.;
      solution_store[N][0] = 0.;

      std::map<types::global_dof_index, double> boundary_values_l,
	  boundary_values_r;
      VectorTools::interpolate_boundary_values (dof_handler, 0,
						ZeroFunction<2> (),
						boundary_values_l);
      MatrixTools::apply_boundary_values (boundary_values_l, system_matrix,
					  Xnp1, system_rhs);
      VectorTools::interpolate_boundary_values (dof_handler, 1,
						ConstantFunction<2> (0.),
						boundary_values_r);
      MatrixTools::apply_boundary_values (boundary_values_r, system_matrix,
					  Xnp1, system_rhs);

      for (timestep_number = 1; timestep_number < total_time_step;
	  timestep_number++)
	{
	  time += time_step;
	  A.vmult (system_rhs, Xn);

	  // Manually apply BC;
	  for (auto itr = boundary_values_l.begin ();
	      itr != boundary_values_l.end (); ++itr)
	    {
	      system_rhs[itr->first] = itr->second;
	    }
	  for (auto itr = boundary_values_r.begin ();
	      itr != boundary_values_r.end (); ++itr)
	    {
	      system_rhs[itr->first] = itr->second;
	    }

	  solve_time_step ();
	  Xn = Xnp1;
	  solution_store[0][timestep_number] = time;
	  for (int i = 0; i < N; i++)
	    {
	      solution_store[i + 1][timestep_number] =
		  Xnp1[solution_table_1D_to_2D.find (i)->second];
	    }

	}

//       write solution;

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
//	}
//      scanf ("%d", &de);

      /*   integrate for f0 use romint   */
      double v_for_romint[total_time_step];
      for (int i = 0; i < N; i++)
	{
	  for (int j = 0; j < total_time_step; j++)
	    {
	      v_for_romint[j] = solution_store[i + 1][j]
		  * solution_store[i + 1][total_time_step - j];
	    }
	  f0[i] = romint (v_for_romint, total_time_step,
			  1. / (total_time_step - 1));
//	  printf ("f0[%d]=%2.15f\n", i, f0[i]);
	}
//      scanf ("%d", &de);

      for (int i = 1; i < N - 1; i++)
	{
	  out[i] = f0_given[i] - f0[i]+yita_full_2D[i];
	}
      return out;
    }

  template<int dim>
    double *
    HeatEquation<dim>::run_experiemnt ()
    {

      GridIn<dim> grid_in;
      grid_in.attach_triangulation (triangulation);
      std::ifstream input_file ("yita_full_2D_N=033.msh");
      Assert(dim == 2, ExcInternalError ());
      grid_in.read_msh (input_file);

      refine_times++;

      std::cout << "Number of active cells: " << triangulation.n_active_cells ()
	  << std::endl;

      f0_given = (double*) realloc (f0_given, N * sizeof(double));
      std::vector<double> x;
      x.resize (N);
      get_x (x);
      get_f0_given (tau, N, x);

      setup_system ();

      double *a;
      return a;
    }
}

void
read_yita_middle_1D (std::vector<double>& yita_middle, char* filename, int& N)
{

  FILE *file;
  file = fopen (filename, "r");
  if (file == NULL)
    {
      fprintf (stderr, "Can't open input file!\n");
      exit (EXIT_FAILURE);
    }
  char buff[255];
  int i;
  double value, x;
  fgets (buff, 255, (FILE*) file);
  sscanf (buff, "N= %d", &N);
  yita_middle.resize (N);
  while (fgets (buff, 255, (FILE*) file))
    {
      sscanf (buff, "%d ,%lf, %lf", &i, &x, &value);
      yita_middle[i] = value;
    }
  fclose (file);
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

      int N;
      std::vector<double> yita_middle; // initial guess, the ends are bounded   // this is the middle of yita_1D, because the boundary are fixed.
      double tau = 0.5302, L = 3.72374; // tau is for calculating f0_given, L is the length.

      // read data from file, also set N;  N=33_for_read.txt
      read_yita_middle_1D (yita_middle, "N=33_for_read.txt", N);
//      N=5;
      f0_given = (double*) malloc (N * sizeof(double)); // this is the ideal f0;

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

//      for (int i = 1; i < N - 1; i++)
//	printf ("x_nr[%d]=%f \n", i, x_nr[i]);
//      scanf ("%d", &de);

      Step26::HeatEquation<2> other (tau, N, 2049, L, x_nr); /* This tells the yita_middle_1D and x_nr has the same address;
       and 2049 are points, 2048 intervals */
      heat_equation_solver = other;
      /*The assignment constructor transfer the address of x_nr
       to heat_equation_solver.yita_middle_1D*/

      /*--------------------------------------------------------------*/
      std::vector<double> interpolated_solution_yita_1D;

      int ii;
      adm (x_nr, N - 2, &ii, SCFT_wrapper, 0);

      return 0;

//      for (int i = 0; i < 1; i++)
//	{
//	  broydn (x_nr, N - 2, &check, SCFT_wrapper);
//	  heat_equation_solver.refine_mesh (interpolated_solution_yita_1D);
//	  free_dvector (x_nr, 1, N - 2);
//	  N = heat_equation_solver.get_N ();
//	  x_nr = dvector (1, N - 2);
//	  for (int i = 1; i < N - 1; i++)
//	    x_nr[i] = interpolated_solution_yita_1D[i];
//	  heat_equation_solver.set_yita_middle_1D (x_nr);
//	  f0_given = (double*) realloc (f0_given, N * sizeof(double));
//	  local_iteration = 0;
//	}

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
