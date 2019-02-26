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
#include <deal.II/numerics/fe_field_function.h>
#include "nrutil.h"
#include "nr.h"

template class std::vector<dealii::Point<2>>;
//explicit instantiation for debug purpose

namespace SCFT
{
  template<int dim>
    HeatEquation<dim>::HeatEquation (double tau, int N, int total_time_step,
				     double L) :
	fe (1), dof_handler (triangulation), tau (tau), time (0.0), time_step (
	    1. / (total_time_step - 1)), timestep_number (0), N (N), total_time_step (
	    total_time_step), n_dof (dof_handler.n_dofs ()), L (L), refine_times (
	    0), iteration (0), mean_field_free_energy (0.0)
    {
      time_step = 1. / (total_time_step - 1);
      solution_store = Matcreate (N + 1, total_time_step);
      f0.reinit (N);
      yita_middle_1D = NULL;
      yita_full_1D.reinit (N);
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
      solution_store = Matcreate (N + 1, total_time_step);
      f0 = other.f0;
      yita_full_1D = other.yita_full_1D; //TODO: make sure it is notdeep copy
      yita_full_2D = other.yita_full_2D;
      local_iteration = other.local_iteration;
      out = other.out;
      iteration = other.iteration;
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
      Point < dim > p;
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

      KellyErrorEstimator < dim
	  > ::estimate (dof_handler, QGauss < dim - 1 > (fe.degree + 1),
			typename FunctionMap<dim>::type (), yita_full_2D,
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
//	  if (cell->refine_flag_set ())
	  cell->set_refine_flag (RefinementCase < dim > ::cut_axis (0));
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
      printf ("N=%d  \n", N);
      Matfree (solution_store);
      solution_store = Matcreate (N + 1, total_time_step + 1);
      f0.reinit (N);
      yita_full_1D.reinit (N);
      out.reinit (N);
      lookup_table_1D_to_2D.clear ();
      lookup_table_2D_to_1D.clear ();
      iteration = 0;
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
      DataOut < dim > data_out;
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
    HeatEquation<dim>::print_and_save_yita_1D ()
    {
      printf ("Solved ! \n Full_1D_solution: N=%d \n", get_N ());

      // get function value and nPlot points
      int nPlot = pow (2, 17);
      std::vector<Point<dim> > vP (nPlot); // stores location that I am interested
      for (unsigned int i = 0; i < nPlot; i++)
	{
	  vP[i][0] = L * i / (nPlot - 1);
	  vP[i][1] = 0.0;
	}

      std::vector<double> detailedSolutionYita1D (nPlot);
      Functions::FEFieldFunction < dim > fefunction (dof_handler, yita_full_2D);
      fefunction.value_list (vP, detailedSolutionYita1D);

      std::vector<double> xp (nPlot);

      for (int i = 0; i < nPlot; i++)
	{
	  xp[i] = vP[i][0];
	}

      //      set_mean_field_free_energy ();
      set_mean_field_free_energy_romint (nPlot, xp, detailedSolutionYita1D);

      FILE * fp;
      // write nPlot points solution
      const std::string filename = "detailedsolution_yita_1D_N="
	  + Utilities::int_to_string (get_N (), 3) + ".txt";

      fp = fopen (filename.c_str (), "w+");
      if (fp == NULL)
	{
	  printf ("cannot create file %s \n;", filename.c_str ());
	  exit (-1);
	}
      fprintf (fp, "N= %d \n", get_N ());

      for (int i = 0; i < nPlot; i++)
	{
	  fprintf (fp, "%i,%2.15f,%2.15f\n", i, xp[i],
		   detailedSolutionYita1D[i]);
	}

      fprintf (fp, "mean_field_free_energy, %2.15f \n", mean_field_free_energy);
      fclose (fp);


      // print and write  solution;
      std::vector<double> x;
      get_x (x);

      const std::string filename2 = "solution_yita_1D_N="
	  + Utilities::int_to_string (get_N (), 3) + ".txt";

      fp = fopen (filename2.c_str (), "w+");
      if (fp == NULL)
	{
	  printf ("cannot create file %s \n;", filename.c_str ());
	  exit (-1);
	}
      fprintf (fp, "N= %d \n", get_N ());
      for (int i = 0; i < N; i++)
	{
	  printf ("solution[%d]=%2.15f \n", i,
		  yita_full_2D[lookup_table_1D_to_2D.find (i)->second]);
	  fprintf (fp, "%d,%2.15f,%2.15f\n", i, x[i],
		   yita_full_2D[lookup_table_1D_to_2D.find (i)->second]);
	}

      printf ("mean_field_free_energy=%2.15f \n", mean_field_free_energy);
      fprintf (fp, "mean_field_free_energy, %2.15f \n", mean_field_free_energy);

      fclose (fp);
      printf ("%s is written. \n", filename.c_str ());
    }

  template<int dim>
    double
    HeatEquation<dim>::get_L ()
    {
      return L;
    }

  template<int dim>
    void
    HeatEquation<dim>::set_mean_field_free_energy ()
    {
      const QGauss<dim> quadrature_formula (fe.degree + 1);
      mean_field_free_energy = 0.;
      FEValues < dim
	  > fe_values (
	      fe,
	      quadrature_formula,
	      update_values | update_gradients | update_JxW_values
		  | update_quadrature_points);
      DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active (),
	  endc = dof_handler.end ();
      const unsigned int n_q_points = quadrature_formula.size (); // which is 4.
      Point < dim > p;
      std::vector<Point<dim> > quadrature (n_q_points);
      std::vector<double> solution_values (quadrature_formula.size ());
      double x, fi;
      double y = L / 33;
      double f0Bar = 0.892581217773656; // calculated from romint, use testFiBar.cc tp calculate it.

      MappingQ1 < dim > mapping;
      for (; cell != endc; ++cell)
	{
	  fe_values.reinit (cell);
	  quadrature = fe_values.get_quadrature_points ();
	  fe_values.get_function_values (yita_full_2D, solution_values);
	  for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
	    {
	      p = quadrature[q_point];
	      x = p[0];
	      if (x >= L / 2)
		x = L - x;
	      if (x <= tau)
		{
		  fi = pow ((exp (4 * tau * x / (tau * tau - x * x)) - 1), 2)
		      / pow (((exp (4 * tau * x / (tau * tau - x * x)) + 1)),
			     2);
		  if (std::isnan (fi))
		    fi = 1.;
		}
	      else
		fi = 1.;
//			printf("fi=%2.15f when x = %2.15f ,p[0]=%2.15f\n", fi, x, p[0]);

	      mean_field_free_energy += solution_values[q_point] * fi
		  * fe_values.JxW (q_point);
	    }
	}

      mean_field_free_energy = mean_field_free_energy / y; // Do I need it, yes
      mean_field_free_energy =
	  (mean_field_free_energy / f0Bar / L + log (f0Bar)) / (-1000.);
    }

  template<int dim>
    void
    HeatEquation<dim>::set_mean_field_free_energy_romint (
	int nPlot, std::vector<double> x,
	std::vector<double> detailedSolutionYita1D)
    {

      double f0Bar = 0.892581217773656;

      std::vector<double> f0_given;
      f0_given.resize (nPlot);
      f0_given.assign (nPlot, 1.0);

      for (int i = 0; i < nPlot; i++)
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

	      f0_given[nPlot - i - 1] = f0_given[i];

	    }
	  else
	    break;
	}

      for (int i = 0; i < nPlot; i++)
	{
	  detailedSolutionYita1D[i] = detailedSolutionYita1D[i] * f0_given[i];
	}

      mean_field_free_energy = romint (&detailedSolutionYita1D[0], nPlot - 1,
				       L / (nPlot - 1));

      mean_field_free_energy =
	  (mean_field_free_energy / f0Bar / L + log (f0Bar)) / (-1000.);
      printf ("integration=%2.15f\n", mean_field_free_energy);

    }

  template<int dim>
    void
    HeatEquation<dim>::set_yita_full_2D ()
    {
      std::vector<double> x;
      get_x (x);
      x.pop_back (); // remove the last and the first elements
      x.erase (x.begin ()); // yita_middle_1D si a double*

      MappingQ1 < dim > mapping;
      std::vector<Point<dim> > support_points (dof_handler.n_dofs ());
      DoFTools::map_dofs_to_support_points < dim
	  > (mapping, dof_handler, support_points);
      double* xp = (double*) malloc (sizeof(double) * support_points.size ());
      double* yp = (double*) malloc (sizeof(double) * support_points.size ());
      for (unsigned int i = 0; i < support_points.size (); i++)
	{
	  xp[i] = support_points[i] (0);
	}

      double m = 0.;
      spline_chen (&x[0], yita_middle_1D, xp, yp, x.size (),
		   support_points.size (), &m);

//      for (int i = 0; i < support_points.size (); i++)
//	printf ("%2.15f\n", yp[i]);

      for (unsigned int i = 0; i < support_points.size (); i++)
	{
	  yita_full_2D[i] = yp[i];
	}

      free (xp);
      free (yp);

//      int de;
//      scanf ("%d", &de);
    }

  template<int dim>
    void
    HeatEquation<dim>::build_lookup_table ()
    {
      std::vector<double> x;
      get_x (x);
      std::vector<double> yita_full_1D_temp (N, 0.);
      for (int i = 1; i < N - 1; i++)
	yita_full_1D_temp[i] = yita_middle_1D[i - 1];

      MappingQ1 < dim > mapping;
      std::vector<Point<dim> > support_points (dof_handler.n_dofs ());
      DoFTools::map_dofs_to_support_points < dim
	  > (mapping, dof_handler, support_points);
      double px;
      for (unsigned int i = 0; i < support_points.size (); i++)
	{
	  px = support_points[i] (0);
	  int index = std::lower_bound (x.begin (), x.end (), px) - x.begin ();
	  lookup_table_1D_to_2D.insert (std::pair<int, int> (index, i));
	  lookup_table_2D_to_1D.insert (std::pair<int, int> (i, index));
	}

      /* print lookup tables */
//      printf("lookup_table_1D_to_2D:\n");
//      for (auto itr = lookup_table_1D_to_2D.begin ();
//	  itr != lookup_table_1D_to_2D.end (); ++itr)
//	{
//	  std::cout << '\t' << itr->first << "---" << itr->second << '\n';
//	}
//      std::cout << std::endl;
//      printf("lookup_table_2D_to_1D:\n");
//      for (auto itr = lookup_table_2D_to_1D.begin ();
//	  itr != lookup_table_2D_to_1D.end (); ++itr)
//	{
//	  std::cout << '\t' << itr->first << "---" << itr->second << '\n';
//	}
//      std::cout << std::endl;
//      int de;
//      scanf("%d",&de);
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
      DoFTools::make_sparsity_pattern (dof_handler, dsp, constraint_matrix,
      /*keep_constrained_dofs = */true);
      sparsity_pattern.copy_from (dsp);
      A.reinit (sparsity_pattern);
      B.reinit (sparsity_pattern);
      system_matrix.reinit (sparsity_pattern);
      C.reinit (sparsity_pattern);
      D.reinit (sparsity_pattern);
      Xnp1.reinit (dof_handler.n_dofs ());
      Xn.reinit (dof_handler.n_dofs ());
      yita_full_2D.reinit (dof_handler.n_dofs ());
      system_rhs.reinit (dof_handler.n_dofs ());
      N = triangulation.n_active_cells () + 1;
      n_dof = dof_handler.n_dofs ();
      solutionBlock.reinit (2 * n_dof);
      systemRhsBlock.reinit (2 * n_dof);
      tmp.reinit (n_dof);
      constraint_matrix.clear ();

      // Manually set SparsityPattern
      DynamicSparsityPattern dspBlock (dof_handler.n_dofs () * 2);
      for (unsigned int i = 0; i < A.m (); i++)
	{
	  SparseMatrix<double>::iterator begin = A.begin (i), end = A.end (i);
	  for (; begin != end; ++begin)
	    {
	      const dealii::SparseMatrixIterators::Accessor<double, false> acc =
		  *begin;
	      unsigned int row = acc.row (), col = acc.column ();
	      dspBlock.add (row, col);
	      dspBlock.add (row, col + n_dof);
	      dspBlock.add (row + n_dof, col);
	      dspBlock.add (row + n_dof, col + n_dof);
	    }
	}
      sparsityPatternBlock.copy_from (dspBlock);
      systemMatrixBlock.reinit (sparsityPatternBlock);

      get_f0_given ();
      update_internal_data ();
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

      VectorTools::interpolate_boundary_values (dof_handler, 0,
						ZeroFunction<dim> (),
						constraint_matrix);

      VectorTools::interpolate_boundary_values (dof_handler, 1,
						ZeroFunction<dim> (),
						constraint_matrix);
      constraint_matrix.close ();

      //      MatrixCreator::create_mass_matrix (dof_handler,
      //					 QGauss<dim> (fe.degree + 1), A);
      const QGauss<dim> quadrature_formula (fe.degree + 1);

      FEValues < dim
	  > fe_values (
	      fe,
	      quadrature_formula,
	      update_values | update_gradients | update_JxW_values
		  | update_quadrature_points);

      const unsigned int dofs_per_cell = fe.dofs_per_cell;
      const unsigned int n_q_points = quadrature_formula.size ();

      FullMatrix<double> cell_A (dofs_per_cell, dofs_per_cell);
      FullMatrix<double> cell_B (dofs_per_cell, dofs_per_cell);
      FullMatrix<double> cell_C (dofs_per_cell, dofs_per_cell);

      std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

      DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active (),
	  endc = dof_handler.end ();

      std::vector<double> old_solution_values (quadrature_formula.size ());
      for (; cell != endc; ++cell)
	{
	  cell_A = 0.;
	  cell_B = 0.;
	  cell_C = 0.;
	  fe_values.reinit (cell);
	  fe_values.get_function_values (yita_full_2D, old_solution_values);

//	std::vector<Point<dim> > qpositions(10);
//	qpositions=fe_values.get_quadrature_points();

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
		  cell_C (i, j) += fe_values.shape_value (i, q_point)
		      * fe_values.shape_value (j, q_point)
		      * old_solution_values[q_point] * fe_values.JxW (q_point);
		}

	  cell->get_dof_indices (local_dof_indices);

	  constraint_matrix.distribute_local_to_global (cell_A,
							local_dof_indices, A);
	  constraint_matrix.distribute_local_to_global (cell_B,
							local_dof_indices, B);
	  constraint_matrix.distribute_local_to_global (cell_C,
							local_dof_indices, C);
	}

      D.copy_from (B);
      D.add (1., C);

      double c01 = (1. / 4. - sqrt (3) / 6.) * time_step, c10 = (1. / 4
	  + sqrt (3) / 6.) * time_step;
      for (unsigned int i = 0; i < A.m (); i++)
	{
	  SparseMatrix<double>::iterator begin = A.begin (i), end = A.end (i);
	  for (; begin != end; ++begin)
	    {
	      const dealii::SparseMatrixIterators::Accessor<double, false> acc =
		  *begin;
	      int row = acc.row (), col = acc.column ();
	      // block(0,0)
	      systemMatrixBlock.set (
		  row, col, A (row, col) + time_step / 4 * D (row, col));
	      // block(0,1)
	      systemMatrixBlock.set (row, col + n_dof, c01 * D (row, col));
	      // block(1,0)
	      systemMatrixBlock.set (row + n_dof, col, c10 * D (row, col));
	      // block(1,1)
	      systemMatrixBlock.set (
		  row + n_dof, col + n_dof,
		  A (row, col) + time_step / 4 * D (row, col));
	    }
	}

      A_direct.initialize (systemMatrixBlock);
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
      DataOut < dim > data_out;

      data_out.attach_dof_handler (dof_handler);
      data_out.add_data_vector (Xnp1, "U");
      data_out.build_patches ();
      const std::string filename = "solution-"
	  + Utilities::int_to_string (timestep_number, 3) + ".vtk";
      std::ofstream output (filename.c_str ());
      data_out.write_vtk (output);
      printf ("%s is written \n", filename.c_str ());
    }

  template class HeatEquation<2> ;
  template class HeatEquation<3> ;

}
