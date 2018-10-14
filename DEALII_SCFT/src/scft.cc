/*
 * scft.cc
 *
 *  Created on: Oct 8, 2018
 *      Author: chen
 */

#include "SCFT.h"
#include "NR_chen.h"
#include <vector>
#include <fstream>

namespace SCFT
{
  template<int dim>
    HeatEquation<dim>::HeatEquation (double tau, int N, int total_time_step,
				     double L) :
	fe (1), dof_handler (triangulation), tau (tau), time (0.0), time_step (
	    1. / (total_time_step - 1)), timestep_number (0), N (N), total_time_step (
	    total_time_step), L (L), refine_times (0)
    {
      time_step = 1. / (total_time_step - 1);
      solution_store = Matcreate (N + 1, total_time_step + 1);
      f0.reinit (N);
      yita_middle_1D = NULL;
      yita_full_1D.reinit (N);
      yita_full_2D.reinit (2 * N);
      out.reinit (N);
      f0_given.reinit (N);
      local_iteration = 0;
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
      time_step = 1. / (total_time_step - 1);
      solution_store = Matcreate (N + 1, total_time_step + 1);
      f0 = other.f0;
      yita_full_1D = other.yita_full_1D; //TODO: make sure it is notdeep copy
      yita_full_2D = other.yita_full_2D;
      local_iteration = other.local_iteration;
      out = other.out;

      return *this;
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
      in.resize (N);
      typename Triangulation<dim>::active_vertex_iterator vertex =
	  triangulation.begin_active_vertex (), endv =
	  triangulation.end_vertex ();
      int i = 0;
      Point<dim> p;
      for (; vertex != endv; vertex++)
	{
	  p = vertex->vertex ();
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
    HeatEquation<dim>::set_refine_times (int in)
    {
      refine_times = in;
    }

  template<int dim>
    void
    HeatEquation<dim>::set_N (int in)
    {
      N = in;
    }

  template<int dim>
    int
    HeatEquation<dim>::get_N () const
    {
      return N;
    }
  template<int dim>
    void
    HeatEquation<dim>::set_local_iteration (int in)
    {
      local_iteration = in;
    }

  template<int dim>
    int
    HeatEquation<dim>::get_local_iteration () const
    {
      return local_iteration;
    }

  template<int dim>
    void
    HeatEquation<dim>::refine_mesh (std::vector<double> oldSolution,
				    std::vector<double> & newSolution)
    {
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

      std::vector<double> x, xp;
      get_x (x);
      triangulation.execute_coarsening_and_refinement ();
      setup_system ();
      get_x (xp);
      newSolution.resize (N);
      spline_chen (&x[1], &oldSolution[1], &xp[1], &newSolution[1],
		   x.size () - 2, xp.size () - 2, NULL);

      refine_times++;
    }

  template<int dim>
    void
    HeatEquation<dim>::update_internal_data ()
    {
      time = 0.;
      timestep_number = 0;
      N = triangulation.n_active_cells () + 1;
      printf ("N=%d  \n", N);
      Matfree (solution_store);
      solution_store = Matcreate (N + 1, total_time_step + 1);
      f0.reinit (N);
      yita_full_1D.reinit (N);
      yita_full_2D.reinit (2 * N);
      out.reinit (N);
      solution_table_1D_to_2D.clear ();
      solution_table_2D_to_1D.clear ();
      solution_table_x_to_2D.clear ();
      solution_table_2D_to_x.clear ();
    }

  template<int dim>
    void
    HeatEquation<dim>::get_f0_given ()
    {
      f0_given.reinit (N);
      for (int i = 0; i < N; i++)
	f0_given[i] = 1.;
      std::vector<double> x;
      get_x (x);
      for (int i = 0; i < N; i++)
	{
	  if (x[i] <= tau)
	    {
	      f0_given[i] = pow (
		  (exp (4 * tau * x[i] / (tau * tau - x[i] * x[i])) - 1), 2)
		  / pow (
		      ((exp (4 * tau * x[i] / (tau * tau - x[i] * x[i])) + 1)),
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
    void
    HeatEquation<dim>::output_results_for_yita_full_2D () const
    {
      DataOut<dim> data_out;
      data_out.attach_dof_handler (dof_handler);
      data_out.add_data_vector (yita_full_2D, "U");
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
    HeatEquation<dim>::print_and_save_yita_1D (std::vector<double> solution)
    {
      printf ("Solved ! \n middle_solution: N=%d \n", get_N ());
      // print and write  solution;
      std::vector<double> x;
      get_x (x);
      FILE * fp;
      char filename[64], num[64];
      sprintf (num, "%d", get_N ());
      strcpy (filename, "solution_yita_1D_N= ");
      strcat (filename, num);
      strcat (filename, ".txt");

      fp = fopen (filename, "w+");
      fprintf (fp, "N= %d \n", get_N ());
      for (unsigned int i = 0; i < solution.size (); i++)
	{
	  printf ("solution[%d]=%2.15f \n", i, solution[i]);
	  fprintf (fp, "%d,%2.15f,%2.15f\n", i, x[i], solution[i]);
	}
      fclose (fp);
      printf ("%s is written. \n", filename);
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
      Xnp1.reinit (dof_handler.n_dofs ());
      Xn.reinit (dof_handler.n_dofs ());
      system_rhs.reinit (dof_handler.n_dofs ());
      N = triangulation.n_active_cells () + 1;
      get_f0_given ();
      update_internal_data ();
      build_solution_table ();
    }

  template<int dim>
    void
    HeatEquation<dim>::assemble_system () /* assemble for A,B,C
     where A:(fi[i],fi[j]) B:(d(fi[i]),d(fi[j])) C:(yita[i]*fi[i],fi[j])
     /*/
    {
      A = 0.;
      B = 0.;
      C = 0.;

//      MatrixCreator::create_mass_matrix (dof_handler,
//					 QGauss<dim> (fe.degree + 1), A);
      const QGauss<dim> quadrature_formula (fe.degree + 1);

      FEValues<dim> fe_values (
	  fe, quadrature_formula,
	  update_values | update_gradients | update_JxW_values);

      const unsigned int dofs_per_cell = fe.dofs_per_cell;
      const unsigned int n_q_points = quadrature_formula.size ();

      FullMatrix<double> cell_A (dofs_per_cell, dofs_per_cell);
      FullMatrix<double> cell_B (dofs_per_cell, dofs_per_cell);
      FullMatrix<double> cell_C (dofs_per_cell, dofs_per_cell);

      std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

      DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active (),
	  endc = dof_handler.end ();

      for (; cell != endc; ++cell)
	{
	  cell_A = 0.;
	  cell_B = 0.;
	  fe_values.reinit (cell);

	  for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
	    for (unsigned int i = 0; i < dofs_per_cell; ++i)
	      for (unsigned int j = 0; j < dofs_per_cell; ++j)
		{
		  cell_A (i, j) += fe_values.shape_value (i, q_point)
		      * fe_values.shape_value (j, q_point)
		      * fe_values.JxW (q_point);
		  cell_B (i, j) += fe_values.shape_grad (i, q_point)
		      * fe_values.shape_grad (j, q_point)
		      * fe_values.JxW (q_point);
		}

	  cell->get_dof_indices (local_dof_indices);
	  cell->get_dof_indices (local_dof_indices);
	  for (unsigned int i = 0; i < dofs_per_cell; ++i)
	    {
	      for (unsigned int j = 0; j < dofs_per_cell; ++j)
		{
		  A.add (local_dof_indices[i], local_dof_indices[j],
			 cell_A (i, j));
		  B.add (local_dof_indices[i], local_dof_indices[j],
			 cell_B (i, j));
		}
	    }
	}

      std::vector<double> yita_full_1D (N, 0.);
      for (int i = 0; i < N - 2; i++)
	yita_full_1D[i + 1] = yita_middle_1D[i];
      for (int i = 0; i < 2 * N; i++)
	{
	  int j = solution_table_2D_to_1D.find (i)->second;
	  yita_full_2D[i] = yita_full_1D[j];
	}
      C.copy_from (A);
      for (unsigned int i = 0; i < C.m (); i++)
	{
	  SparseMatrix<double>::iterator begin = C.begin (i), end = C.end (i);
	  for (; begin != end; ++begin)
	    {
	      begin->value () *= yita_full_2D[i];
	    }
	}
      system_matrix.copy_from (A);
      system_matrix.add (time_step, B);
      system_matrix.add (time_step, C);

    }

  template<int dim>
    void
    HeatEquation<dim>::solve_time_step ()
    {
      SolverControl solver_control (80000, 1e-17);
      SolverCG<> solver (solver_control);
      solver.solve (system_matrix, Xnp1, system_rhs, PreconditionIdentity ());
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

      Point<dim> P;
      std::vector<types::global_dof_index> loc_dof_indices (fe.dofs_per_cell);
      typename DoFHandler<dim>::active_cell_iterator cell =
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

  template class HeatEquation<2> ;
  template class HeatEquation<3> ;
}
