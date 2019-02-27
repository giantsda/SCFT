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
#include <deal.II/fe/fe_system.h>
#include <deal.II/dofs/dof_renumbering.h>

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

//#define CHECK;
int de; // My debug varaibe

std::string fileNameMesh;

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
      Assert (component == 0, ExcIndexRange (component, 0, 1));
      Assert (dim == 2, ExcNotImplemented ());
      if (p[0] == 0. || p[0] == heat_equation_solver.get_L ())
	return 0.;
      else
	return 1.;
    }
}

template<int dim>
  double *
  SCFT::HeatEquation<dim>::run (double* yita_middle_1D_in)
  {
    yita_middle_1D = yita_middle_1D_in; //set yita_middle_1D

    if (refine_times == 0)
      {
	if (0)
	  {
	    std::vector<unsigned int> repetitions;
	    repetitions.push_back (N - 1);
	    repetitions.push_back (1);
	    GridGenerator::subdivided_hyper_rectangle (triangulation,
						       repetitions,
						       Point < 2 > (0.0, 0.0),
						       Point < 2 > (L, L / N),
						       true);
	  }
	else
	  {
	    GridIn < dim > grid_in;
	    grid_in.attach_triangulation (triangulation);
	    std::ifstream input_file (fileNameMesh.c_str ());
	    Assert (dim == 2, ExcInternalError ());
	    grid_in.read_msh (input_file);
	  }
	setup_system (); // The first time, the triangulation is generated and system is set up. The
	// Next time, it is setup in the refine();
	refine_times++;
      }
    if (iteration == 0)
      {
	build_lookup_table ();
      }

    set_yita_full_2D ();
    assemble_system ();

    VectorTools::interpolate (dof_handler, Initial_condition<dim> (), Xn);
    Xnp1 = Xn;

    for (int i = 2; i < N; i++)
      solution_store[i][0] = 1.;
    solution_store[0][0] = 0.;
    solution_store[1][0] = 0.; // have to initialize it or it may be a NaN stored in it.
    solution_store[N][0] = 0.;

    time = 0.;
    for (timestep_number = 1; timestep_number < total_time_step;
	timestep_number++)
      {
	time += time_step;
	D.vmult (tmp, Xn);
	tmp *= -1.;
	for (int i = 0; i < n_dof; i++)
	  {
	    systemRhsBlock[i] = tmp[i];
	    systemRhsBlock[i + n_dof] = tmp[i];
	  }
	A_direct.vmult (solutionBlock, systemRhsBlock);
	for (unsigned int i = 0; i < Xnp1.size (); i++)
	  Xnp1[i] = Xn[i]
	      + 0.5 * time_step * (solutionBlock[i] + solutionBlock[i + n_dof]);

	Xn = Xnp1;
	solution_store[0][timestep_number] = time;

	for (int i = 0; i < N; i++)
	  {
	    solution_store[i + 1][timestep_number] =
		Xnp1[lookup_table_1D_to_2D.find (i)->second];
#ifdef CHECK
	    if (std::isnan (solution_store[i + 1][timestep_number]))
	      {
		printf ("Got an NAN from solution_store:\n");
		printf ("solution_store[%d][%d]=%2.15f\n", i + 1,
		    timestep_number,
		    solution_store[i + 1][timestep_number]);
		exit (-1);
	      }
#endif
	  }

      }

    /* write solution; */
    int write = 0;

    if (write)
      {
	FILE * fp;
	fp = fopen ("solution_store.txt", "w+");
	for (int i = 0; i < N + 1; i++)
	  {
	    for (int j = 0; j < total_time_step; j++)
	      fprintf (fp, "%2.15f,", solution_store[i][j]);
	    fprintf (fp, "\n");
	  }
	fclose (fp);
      }
//    scanf ("%d", &de);
    /*   integrate for f0 use romint   */
    double v_for_romint[total_time_step];
    for (int i = 0; i < N; i++)
      {
	for (int j = 0; j < total_time_step; j++)
	  {
	    v_for_romint[j] = solution_store[i + 1][j]
		* solution_store[i + 1][total_time_step - 1 - j];
	  }
	f0[i] = romint (v_for_romint, total_time_step - 1,
			1. / (total_time_step - 1));
#ifdef CHECK
	if (std::isnan (f0[i]))
	  {
	    printf ("Got an NAN from romint(): f0[%d]=%2.15f\n", i, f0[i]);
	    for (int k = 0; k < total_time_step; k++)
	    printf ("v_for_romint[%d]=%2.15f\n", k, v_for_romint[k]);
	    printf ("solution_store[%d][0]=%f\n", i + 1,
		solution_store[i + 1][0]);
	    printf ("solution_store[%d][%d]=%f\n", i + 1, total_time_step - 0,
		solution_store[i + 1][total_time_step - 0]);
	    exit (-1);
	  }
#endif
      }
    //      scanf ("%d", &de);

    for (int i = 0; i < N; i++)
      {
	out[i] = f0_given[i] - f0[i]; // +yita_full_1D[i];  // for adm // so f0 and f0_given are full sized and so out is full sized.
      }
    iteration++;
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
  //  for (int i = 0; i < N; i++)
  //    printf (
  //	"in[%d]=%2.15f ; out[%d]=%2.15f ; local_iteration:%d refine_times:%d \n",
  //	i, in[i], i, out[i], local_interation,
  //	heat_equation_solver.get_refine_times ());
  //  fflush (stdout);
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
main (int argc, char *argv[])
{
  try
    {
      using namespace dealii;
      using namespace SCFT;

      int N;
      std::vector<double> x_old; // initial guess, the ends are bounded   // this is the middle of yita_1D, because the boundary are fixed.
      double tau = 5.30252230020752e-01, L = 3.72374357332160; // tau is for calculating f0_given, L is the length.
      std::string fileNameResults;
      std::string NStore = argv[1];
      fileNameResults = "solution_yita_1D_N=";
      fileNameResults += NStore;
      fileNameResults += ".txt";
      fileNameMesh = "yita_full_2D_N=";
      fileNameMesh += NStore;
      fileNameMesh += ".msh";
      printf ("%s and %s \n", fileNameResults.c_str (), fileNameMesh.c_str ());
      read_yita_middle_1D (x_old, fileNameResults.c_str (), N);
      HeatEquation<2> other (tau, N, 2049, L); /* 2049 are points, 2048 intervals */
      heat_equation_solver = other; // I need this global class to do stuffs
      double* out = heat_equation_solver.run (&x_old[1]);
      double max = -1.;
      for (int i = 0; i < N; i++)
	{
	  if (fabs (out[i]) >= max)
	    max = fabs (out[i]);
	}

      heat_equation_solver.print_and_save_yita_1D ();
      heat_equation_solver.output_results_for_yita_full_2D ();
      heat_equation_solver.output_mesh ();

//      FILE * fp;
//      fp = fopen (fileNameResults.c_str (), "a");
//      if (fp == NULL)
//	{
//	  printf ("cannot create file %s \n;", fileNameResults.c_str ());
//	  exit (-1);
//	}
//      fprintf (fp, "solverError=%2.10E \n", max);
//      fclose (fp);
//
//      const std::string filename = "detailedsolution_yita_1D_N=" + NStore
//	  + ".txt";
//      fp = fopen (filename.c_str (), "a");
//      if (fp == NULL)
//	{
//	  printf ("cannot create file %s \n;", filename.c_str ());
//	  exit (-1);
//	}
//      fprintf (fp, "solverError=%2.10E \n", max);
//      fclose (fp);

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
