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
#include "SCFT.h"
#include "scft_urtil.h"
#include <stdio.h>
#include <math.h>
#include "nr.h"
#include "nrutil.h"
#include <stdlib.h>
#include <fstream>
#include <vector>

#define BROYDN

int de; // My debug varaibe

SCFT::HeatEquation<2> heat_equation_solver;
namespace dealii
{
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
}

template<int dim>
  double *
  SCFT::HeatEquation<dim>::run (double* yita_middle_1D_in)
  {
    std::vector<double> yita_full_1D (N, 0.);
    for (int i = 0; i < N - 2; i++)
      yita_full_1D[i + 1] = yita_middle_1D_in[i];
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

    // convert yita_full_1D to yita_full_2D;
    for (int i = 0; i < 2 * N; i++)
      {
	int j = solution_table_2D_to_1D.find (i)->second;
	yita_full_2D[i] = yita_full_1D[j];
      }

//      for (int i = 0; i < 2 * N; i++)
//	printf ("yita_full_2D[%d]=%2.15f \n", i, yita_full_2D[i]);
//      scanf ("%d", &de);

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
    MatrixTools::apply_boundary_values (boundary_values_l, system_matrix, Xnp1,
					system_rhs);
    VectorTools::interpolate_boundary_values (dof_handler, 1,
					      ConstantFunction<2> (0.),
					      boundary_values_r);
    MatrixTools::apply_boundary_values (boundary_values_r, system_matrix, Xnp1,
					system_rhs);
    time = 0.;
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

//      {
//	FILE * fp;
//	fp = fopen ("solution_store.txt", "w+");
//	for (int i = 0; i < N + 1; i++)
//	  {
//	    for (int j = 0; j < total_time_step; j++)
//	      fprintf (fp, "%2.15f,", solution_store[i][j]);
//	    fprintf (fp, "\n");
//	  }
//
//	fclose (fp);
//      }
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

    for (int i = 0; i < N; i++)
      {
	out[i] = f0_given[i] - f0[i];// +yita_full_1D[i];  // for adm // so f0 and f0_given are full sized and so out is full sized.
      }
    return &out[1];
  }

void
SCFT_wrapper (int N, double * in, double * out)
{
#ifdef BROYDN
  double* res = heat_equation_solver.run (&in[1]);
#else
  double* res = heat_equation_solver.run (&in[0]);
#endif
  N = heat_equation_solver.get_N ();

#ifdef BROYDN
  for (int i = 0; i < N - 2; i++)
    out[i + 1] = res[i];
#else
  for (int i = 0; i < N - 2; i++)
  out[i] = res[i];
#endif
  int local_interation = heat_equation_solver.get_local_iteration ();
  for (int i = 0; i < N; i++)
    printf (
	"in[%d]=%2.15f ; out[%d]=%2.15f ; local_iteration:%d refine_times:%d \n",
	i, in[i], i, out[i], local_interation,
	heat_equation_solver.get_refine_times ());
  fflush (stdout);
  heat_equation_solver.set_local_iteration (local_interation + 1);
}

template class std::vector<double>;
// enable std::vector::size()

#ifdef BROYDN
// stuff for broydn
double* f0_given;
double **qt, **r, *d, /* qt[1:n,1:n], r[1:n,1:n] and d[1:n] must be allocated in the calling program */
err; /* Passed in as the convergence criterion, and returned with the actual residual error */
int funcerr, /* Flag for error in evaluating the function in Broyden method */
jc; /* For re-use of Jacobian. jc denotes the status of Jacobian calculation: 0 for not calculated,
 1 for previously calculated, 2 for currently calculated */
int PRINT;
#endif

int
main ()
{
  try
    {
      using namespace dealii;
      using namespace SCFT;

      int N;
      std::vector<double> x_old; // initial guess, the ends are bounded   // this is the middle of yita_1D, because the boundary are fixed.
      double tau = 0.5302, L = 3.72374; // tau is for calculating f0_given, L is the length.
      read_yita_middle_1D (x_old, "inputFiles/N=33_for_read.txt", N); // read data from file, also set N;
      HeatEquation<2> other (tau, N, 2049, L); /* 2049 are points, 2048 intervals */
      heat_equation_solver = other; // I need this global class to do stuffs
      std::vector<double> interpolated_solution_yita_1D;
//      double* out=heat_equation_solver.run (&x_old[0]);

//      for(int i=0;i<N;i++)
//	printf("out[%d]=%2.15f\n",i,out[i]);
//      /*--------------------------------------------------------------*/

#ifdef BROYDN
      int check = 1;
      qt = dmatrix (1, N - 2, 1, N - 2);
      r = dmatrix (1, N - 2, 1, N - 2);
      d = dvector (1, N - 2);
      jc = 0;
      err = 0.00000001;
#endif

      for (int i = 0; i < 10; i++)
	{
//	  adm_chen (&SCFT_wrapper, &x_old[1], 1e-7, 3000, N - 2);
	  broydn (&x_old[0], N - 2, &check, SCFT_wrapper);
	  heat_equation_solver.print_and_save_yita_1D (x_old);
	  heat_equation_solver.output_results_for_yita_full_2D ();
	  heat_equation_solver.output_mesh ();
	  heat_equation_solver.refine_mesh (interpolated_solution_yita_1D);
	  N = heat_equation_solver.get_N ();
#ifdef BROYDN
	  free_dmatrix (qt, 1, N - 2, 1, N - 2);
	  qt = dmatrix (1, N - 2, 1, N - 2);
	  free_dmatrix (r, 1, N - 2, 1, N - 2);
	  r = dmatrix (1, N - 2, 1, N - 2);
	  free_dvector (d, 1, N - 2);
	  d = dvector (1, N - 2);
	  jc = 0;
#endif
	  x_old = interpolated_solution_yita_1D;
	  heat_equation_solver.set_local_iteration (0);
	}
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
