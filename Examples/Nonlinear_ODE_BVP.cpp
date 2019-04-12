/// \file  Nonlinear_ODE_BVP.cpp
/// \ingroup Examples
/// The equation \f[ (f(x)-1)f''(x) + \left( \frac{2}{x^2} - 1 \right)f'(x)^2
/// = 1 \f] is solved subject to the boundary conditions \f$ f(0) = 2 \f$ and
/// \f$ f \left(\frac{\sqrt{3}}{2} \right) = \frac{3}{2}. \f$ The system is
/// solved as a set of first-order equations in the form \f[ \begin{pmatrix}
/// 1 & 0 \\ 0 & f(x) - 1 \end{pmatrix} \begin{bmatrix} f(x) \\ f'(x)
/// \end{bmatrix}' = \begin{bmatrix} f'(x) \\ 1 + \left(1 - \frac{2}{x^2}
/// \right)f'(x)^2 \end{bmatrix}, \f] 
/// The exact solution is given by \f$ f(x) = 1 + \sqrt{1-x^2}. \f$

#include "Luna/Core"
#include "Luna/ODE"

// ODE enumeration
enum{ f, fd };

namespace Luna
{
  class nonlinear_equation : public Equation_1matrix<double>
	{

		public:
			// The nonlinear equation is 2nd order
			nonlinear_equation() : Equation_1matrix<double> ( 2 ) {}

			// Define the equation
			void residual_fn( const Vector<double>& u, Vector<double>& F  ) const
			{
        double x( coord( 0 ) );
				F[ f ]   = u[ fd ];
				F[ fd ]  = 1. + ( 1. - ( 2. / ( x * x ) ) ) * u[ fd ] * u[ fd ];
			}

      void matrix0( const Vector<double>&u, Matrix<double> &m ) const
      {
        m( 0, 0 ) = 1.0;
        m( 1, 1 ) = u[ f ] - 1.0;
      }
	};

  class left_BC : public Residual<double>
  {
      public:
        // f(0) = 2
        left_BC() : Residual<double> ( 1, 2 ) {}

        void residual_fn( const Vector<double> &z, Vector<double> &B ) const
        {
            B[ 0 ] = z[ f ] - 2.0;
        }
  };

  class right_BC : public Residual<double>
  {
      public:
        // f(sqrt(3)/2) = 3/2
        right_BC() : Residual<double> ( 1, 2 ) {}

        void residual_fn( const Vector<double> &z, Vector<double> &B ) const
        {
            B[ 0 ] = z[ f ] - 1.5;
        }
  };

}

using namespace std;
using namespace Luna;

int main()
{
  cout << "------------------- Nonlinear_ODE_BVP ------------------" << endl;
  // Setup the domain
  double left( 0.0 );                // Left boundary
  double right( sqrt(3.0) / 2.0 );   // Right boundary
	std::size_t N( 50 );               // Number of nodes
	Vector<double> nodes;						   // Vector of nodes ( uniformly spaced )
	nodes.linspace( left, right, N);

  cout << " * The nonlinear boundary problem is solved on a uniform" << endl
       << " * grid in the domain " << left << " < x < " << right
       << " containing N = " << N << endl << " * nodal points." << endl;

  // Create instances of the equation and BCs
  nonlinear_equation equation;
  left_BC left_bc;
  right_BC right_bc;

  // Create boundary value problem
  ODE_BVP<double> ode( &equation, nodes, &left_bc, &right_bc );

  /* ----- Set the initial guess ----- */
  // f = 2 - (1/sqrt(3))x (linear fit to boundary conditions)
  double gradient( - 1.0 / sqrt(3) );
	for (std::size_t j=0; j < N; ++j )
	{
		double x = nodes[ j ];
		ode.solution()( j, f )  = 2.0 + gradient * x;
    ode.solution()( j, fd ) = gradient;
	}

  ode.tolerance() = 1e-8;
  ode.max_iterations() = 50;
  ode.set_monitor_det( false );

  // Solve
  ode.solve_bvp();

  nodes = ode.solution().nodes();
  Mesh1D<double, double> exact_solution( nodes, 2 );
  for (std::size_t j=0; j < nodes.size(); ++j )
	{
		double x = nodes[ j ];
		exact_solution( j, f )  = 1.0 + sqrt( 1.0 - x * x );
    exact_solution( j, fd ) = - x / sqrt( 1.0 - x * x );
	}

  // Output the solution data
  ode.solution().output( "./DATA/ODE_BVP_numerical_solution.dat" );
  exact_solution.output( "./DATA/ODE_BVP_exact_solution.dat" );

  cout << " * To plot the numerical solution, the exact solution " << endl
       << " * and the error between the two run:" << endl
       << "python Plotting/Nonlinear_ODE_BVP_plot.py" << endl;
  cout << "--- FINISHED ---" << endl;
}
