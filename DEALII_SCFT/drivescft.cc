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

//#define BROYDN

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

namespace SCFT
{
  class SCFTsolve : public Subscriptor
  {
  public:
    SCFTsolve (const BlockSparseMatrix<double> &A);
    void
    vmult (Vector<double> &dst, const Vector<double> &src) const;
  private:
    const SmartPointer<const BlockSparseMatrix<double> > system_matrix;
    mutable Vector<double> srcUp, srcLow, resUp, resLow;
  };

  SCFTsolve::SCFTsolve (const BlockSparseMatrix<double> &A) :
      system_matrix (&A), srcUp (A.block (0, 0).m ()), srcLow (
	  A.block (0, 0).m ()), resUp (A.block (0, 0).m ()), resLow (
	  A.block (0, 0).m ())
  {
  }

  void
  SCFTsolve::vmult (Vector<double> &dst, const Vector<double> &src) const
  {
    int n_dof = src.size () / 2;
    for (int i = 0; i < n_dof; i++)
      {
	srcUp[i] = src[i];
	srcLow[i] = src[i + n_dof];
      }
    Vector<double> tmp1, tmp2;
    tmp1.reinit (n_dof);
    tmp2.reinit (n_dof);

    system_matrix->block (0, 0).vmult (tmp1, srcUp);
    system_matrix->block (0, 1).vmult (tmp2, srcLow);
    resUp = tmp1;
    resUp += tmp2;
    system_matrix->block (1, 0).vmult (tmp1, srcUp);
    system_matrix->block (1, 1).vmult (tmp2, srcLow);
    resLow = tmp1;
    resLow += tmp2;

    for (int i = 0; i < n_dof; i++)
      {
	dst[i] = resUp[i];
	dst[i + n_dof] = resLow[i];
      }
  }

}

template<int dim>
  double *
  SCFT::HeatEquation<dim>::run (double* yita_middle_1D_in)
  {
    yita_middle_1D = yita_middle_1D_in; //set yita_middle_1D

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
    assemble_system ();

    VectorTools::interpolate (dof_handler, Initial_condition<dim> (), Xn);
    Xnp1 = Xn;

    for (int i = 2; i < N; i++)
      solution_store[i][0] = 1.;
    solution_store[0][0] = 0.;
    solution_store[N][0] = 0.;

    FESystem<dim> fe_s (FE_Q<dim> (1), 1, FE_Q<dim> (1), 1);

//    dof_handler.distribute_dofs (fe_s);
//    DynamicSparsityPattern dsp (dof_handler.n_dofs (), dof_handler.n_dofs ());
//    DoFTools::make_sparsity_pattern (dof_handler, dsp);
//    sparsity_pattern.copy_from (dsp);
//    std::ofstream outf ("block.svg");
//    sparsity_pattern.print_svg (outf);
//    scanf("%d",&de);

    BlockSparseMatrix<double> system_matrix_s;
    BlockDynamicSparsityPattern dsp_s (2, 2);
    const unsigned int n_d = dof_handler.n_dofs ();
    dsp_s.block (0, 0).reinit (n_d, n_d);
    dsp_s.block (1, 0).reinit (n_d, n_d);
    dsp_s.block (0, 1).reinit (n_d, n_d);
    dsp_s.block (1, 1).reinit (n_d, n_d);
    dsp_s.collect_sizes ();
    DoFHandler<dim> dof_handler_s (triangulation);
    dof_handler_s.distribute_dofs (fe_s);
    DoFRenumbering::component_wise (dof_handler_s);
    DoFTools::make_sparsity_pattern (dof_handler_s, dsp_s);
    BlockSparsityPattern sparsity_pattern_s;
    sparsity_pattern_s.copy_from (dsp_s);

    system_matrix_s.reinit (sparsity_pattern_s);
    Vector<double> solution_s, system_rhs_s;
    solution_s.reinit (2 * n_d);
    system_rhs_s.reinit (2 * n_d);

    SparseMatrix<double> D; // D=B+C;
    D.reinit (sparsity_pattern);
    D.copy_from (B);
    D.add (1., C);

    double c01 = 1. / 4. - sqrt (3) / 6., c10 = 1. / 4 + sqrt (3) / 6.;
    for (unsigned int i = 0; i < A.m (); i++)
      {
	SparseMatrix<double>::iterator begin = A.begin (i), end = A.end (i);
	for (; begin != end; ++begin)
	  {
	    const dealii::SparseMatrixIterators::Accessor<double, false> acc =
		*begin;
	    int row = acc.row (), col = acc.column ();
	    // block(0,0)
	    system_matrix_s.set (row, col,
				 A (row, col) + time_step / 4 * D (row, col));
	    // block(0,1)
	    system_matrix_s.set (row, col + n_d, c01 * D (row, col));
	    // block(1,0)
	    system_matrix_s.set (row + n_d, col, c10 * D (row, col));
	    // block(1,1)
	    system_matrix_s.set (row + n_d, col + n_d,
				 A (row, col) + time_step / 4 * D (row, col));
	  }
      }
    printf ("B:\n");
    B.print (std::cout);
    printf ("C:\n");
    C.print (std::cout);
    printf ("D:\n");
    D.print (std::cout);
    printf ("system_matrix_s:\n");
    system_matrix_s.print (std::cout);
    scanf ("%d", &de);
    // assmble rhs
    Vector<double> tmp;
    tmp.reinit (n_d);
    D.vmult (tmp, Xn);
    tmp *= -1.;
    for (unsigned int i = 0; i < n_d; i++)
      {
	system_rhs_s[i] = tmp[i];
	system_rhs_s[i + n_d] = tmp[i];
      }

    SCFT::SCFTsolve SCFT_solve (system_matrix_s);
    SolverControl solver_control (80000, 1e-17);

    SolverCG<> cg (solver_control);

    cg.solve (SCFT_solve, solution_s, system_rhs_s, PreconditionIdentity ());

    // solve;

//    SolverControl solver_control (80000, 1e-17);
//    SolverCG<> solver (solver_control);
//    solver.solve (system_matrix_s, solution_s, solution_s,
//		  PreconditionIdentity ());

    for (int i = 0; i < N; i++)
      {
	out[i] = f0_given[i] - f0[i]; // +yita_full_1D[i];  // for adm // so f0 and f0_given are full sized and so out is full sized.
      }
    return &out[1];
  }

template<int dim>
  unsigned int
  SCFT::HeatEquation<dim>::embedded_explicit_method (
      const TimeStepping::runge_kutta_method method,
      const unsigned int n_time_steps, const double initial_time,
      const double final_time)
  {
    double time_step = (final_time - initial_time)
	/ static_cast<double> (n_time_steps);
    double time = initial_time;
    const double coarsen_param = 1.2;
    const double refine_param = 0.8;
    const double min_delta = 1e-8;
    const double max_delta = 1 * time_step;
    const double refine_tol = 1e-1;
    const double coarsen_tol = 1e-5;

    TimeStepping::EmbeddedExplicitRungeKutta<Vector<double> > embedded_explicit_runge_kutta (
	method, coarsen_param, refine_param, min_delta, max_delta, refine_tol,
	coarsen_tol);

    Vector<double> X_1D (N);

    unsigned int n_steps = 0;
    while (time < final_time)
      {
	if (time + time_step > final_time)
	  time_step = final_time - time;

	printf ("time=%f:Before solving, X=:\n\n", time);
	for (int i = 0; i < N; i++)
	  X_1D[i] = Xnp1[solution_table_1D_to_2D.find (i)->second];

	X_1D.print (std::cout);
	time = embedded_explicit_runge_kutta.evolve_one_time_step (
	    std::bind (&HeatEquation<dim>::evaluate_diffusion, this,
		       std::placeholders::_1, std::placeholders::_2),
	    time, time_step, Xnp1);

	printf ("time=%f:After solving, X=:\n\n", time);
	for (int i = 0; i < N; i++)
	  X_1D[i] = Xnp1[solution_table_1D_to_2D.find (i)->second];

	X_1D.print (std::cout);

	scanf ("%d", &de);

	solution_store[0][timestep_number] = time;
	for (int i = 0; i < N; i++)
	  {
	    solution_store[i + 1][timestep_number] =
		Xnp1[solution_table_1D_to_2D.find (i)->second];
	  }

	time_step = embedded_explicit_runge_kutta.get_status ().delta_t_guess;
	++n_steps;
	timestep_number++;
      }

    return n_steps;
  }

template<int dim>
  dealii::Vector<double>
  SCFT::HeatEquation<dim>::evaluate_diffusion (const double time,
					       const Vector<double> &y) const // evaluate inv(A)*(-B*y-C*y)
  {
    Vector<double> tmp1 (dof_handler.n_dofs ());
    tmp1 = 0.;
    Vector<double> tmp2 (dof_handler.n_dofs ());
    tmp2 = 0.;

    const double factor = -1.;

    B.vmult (tmp1, y);
    C.vmult (tmp2, y);
    tmp1 *= -1.;
    tmp1 -= tmp2;

    Vector<double> value (dof_handler.n_dofs ());

    inverse_mass_matrix.vmult (value, tmp1);

    return value;
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
      N = 33;
      HeatEquation<2> other (tau, N, 2049, L); /* 2049 are points, 2048 intervals */
      heat_equation_solver = other; // I need this global class to do stuffs
      std::vector<double> interpolated_solution_yita_1D;

      x_old.clear ();
      x_old.resize (N, 0.);
      double* out = heat_equation_solver.run (&x_old[0]);

      for (int i = 0; i < N; i++)
	printf ("out[%d]=%2.15f\n", i, out[i]);
      /*--------------------------------------------------------------*/

      return 0;

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
	  adm_chen (&SCFT_wrapper, &x_old[1], 1e-1, 200, N - 2, 0.99, 2);
	  adm_chen (&SCFT_wrapper, &x_old[1], 1e-4, 400, N - 2, 0.9, 5);
	  adm_chen (&SCFT_wrapper, &x_old[1], 1e-7, 800, N - 2, 0.9, 15);
	  adm_chen (&SCFT_wrapper, &x_old[1], 1e-7, 1000, N - 2, 0.9, 30);
//	  broydn (&x_old[0], N - 2, &check, SCFT_wrapper);
	  heat_equation_solver.print_and_save_yita_1D (x_old);
	  heat_equation_solver.output_results_for_yita_full_2D ();
	  heat_equation_solver.output_mesh ();
	  heat_equation_solver.refine_mesh (x_old,
					    interpolated_solution_yita_1D);
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
