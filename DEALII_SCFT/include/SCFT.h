/*
 * SCFT.h
 *
 *  Created on: Oct 8, 2018
 *      Author: chen
 */

#ifndef INCLUDE_SCFT_H_
#define INCLUDE_SCFT_H_
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
#include "NR_chen.h"

namespace SCFT
{

  using namespace dealii;

  template<int dim>
    class HeatEquation
    {
    public:
      HeatEquation (double tau, int N, int total_time_step, double L);
      HeatEquation ();
      HeatEquation (HeatEquation&& other);
      HeatEquation&
      operator = (HeatEquation &other);  // assignment constructor
      ~HeatEquation () = default;
      int
      get_refine_times ();
      void
      get_x (std::vector<double> & in); // for a fake 1D problem, it returns the x for the 1D problem.
      void
      set_refine_times (int a);
      void
      set_N (int in);
      int
      get_N () const;
      void
      set_local_iteration (int in);
      int
      get_local_iteration () const;
      void
      refine_mesh (std::vector<double> & in);
      void
      update_internal_data ();
      double *
      run (double* yita_full_1D_in);
      double *
      run_experiemnt ();
      void
      get_f0_given ();
      void
      output_results_for_yita_full_2D () const;
      void
      output_mesh ();
      void
      print_and_save_yita_1D ();
    private:
      void
      setup_system ();
      void
      solve_time_step ();
      void
      output_results () const;
      void
      build_solution_table ();

    private:
      // private data
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
      int N, total_time_step, local_iteration;
      double** solution_store;
      double L;
      Vector<double> f0;
      Vector<double> yita_full_1D;
      Vector<double> yita_full_2D;
      Vector<double> out;
      Vector<double> f0_given;
      int refine_times;
      std::map<int, int> solution_table_1D_to_2D;
      std::map<int, int> solution_table_2D_to_1D;
      std::map<double, int> solution_table_x_to_2D;
      std::map<int, double> solution_table_2D_to_x;
    };

}

#endif /* INCLUDE_SCFT_H_ */
