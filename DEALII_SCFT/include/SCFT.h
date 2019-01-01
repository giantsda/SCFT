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
#include <vector>

#include <deal.II/base/time_stepping.h>
#include <deal.II/lac/sparse_direct.h>

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
      refine_mesh (std::vector<double> oldSolution,
		   std::vector<double> & newSolution);
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
      set_yita_full_2D (); // use cubic spline interpotation
      void
      build_lookup_table ();
      void
      print_and_save_yita_1D (std::vector<double> solution);

    private:
      void
      setup_system ();
      void
      assemble_system ();
      void
      solve_time_step ();
      void
      output_results () const;

    private:
      // private data
      Triangulation<dim> triangulation;
      FE_Q<dim> fe;
      DoFHandler<dim> dof_handler;
      ConstraintMatrix constraint_matrix;

      SparsityPattern sparsity_pattern;
      SparsityPattern sparsityPatternBlock;
      SparseMatrix<double> A; // A: mass matrix (fi[i],fi[j])
      SparseMatrix<double> B; // B: laplace_matrix (d(fi[i]),d(fi[j]))
      SparseMatrix<double> system_matrix;
      SparseMatrix<double> C; // C: (yita[i]*fi[i],fi[j])
      SparseMatrix<double> D; // D: B+C
      SparseMatrix<double> systemMatrixBlock;

      Vector<double> Xnp1;
      Vector<double> Xn;
      Vector<double> system_rhs;
      Vector<double> solutionBlock, systemRhsBlock, tmp;

      double tau;
      double time;
      double time_step;
      int timestep_number;
      int N, total_time_step, local_iteration;
      int n_dof;
      double** solution_store;
      double L;
      Vector<double> f0;
      double* yita_middle_1D;
      Vector<double> yita_full_1D;
      Vector<double> yita_full_2D;
      Vector<double> out;
      Vector<double> f0_given;
      int refine_times;
      std::map<int, int> lookup_table_1D_to_2D;
      std::map<int, int> lookup_table_2D_to_1D;
      SparseDirectUMFPACK A_direct;
      int iteration;
    };

}

#endif /* INCLUDE_SCFT_H_ */
