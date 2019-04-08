/// \file  ODE_BVP_test.cpp
/// \ingroup Examples
/// TODO description / rename file

#include "Luna/Core"
#include "Luna/ODE"

// ODE enumeration
enum{ y, yd };

namespace Luna
{
  class test_equation : public Equation_1matrix<double>
	{

		public:
			// The test equation is 2nd order ( y'' + 4y = 0 )
			test_equation() : Equation_1matrix<double> ( 2 ) {}

			// Define the equation
			void residual_fn( const Vector<double>& u, Vector<double>& F  ) const
			{
        double x( coord( 0 ) );
				F[ y ]   = u[ yd ];
				F[ yd ]  = - 4 * u[ y ];
			}

      void matrix0( const Vector<double>&x, Matrix<double> &m ) const
      {
        m.eye();
      }
	};

  class left_BC : public Residual<double>
  {
      public:
        // y(0) = -2
        left_BC() : Residual<double> ( 1, 2 ) {}

        void residual_fn( const Vector<double> &z, Vector<double> &B ) const
        {
            B[ 0 ] = z[ y ] + 2.0;
        }
  };

  class right_BC : public Residual<double>
  {
      public:
        // y(pi/4) = 10
        right_BC() : Residual<double> ( 1, 2 ) {}

        void residual_fn( const Vector<double> &z, Vector<double> &B ) const
        {
            B[ 0 ] = z[ y ] - 10.0;
        }
  };

}

using namespace std;
using namespace Luna;

int main()
{
  cout << "--------------------- ODE_BVP --------------------" << endl;

  // Setup the domain

  double left( 0.0 );               // Left boundary
  double right( M_PI / 4.0 );			  // Right boundary
	std::size_t N( 1000 );            // Number of nodes
	Vector<double> nodes;						  // Vector of nodes ( uniformly spaced )
	nodes.linspace( left, right, N);

  // Create instances of the equation and BCs
  test_equation equation;
  left_BC left_bc;
  right_BC right_bc;

  // Create boundary value problem
  ODE_BVP<double> ode( &equation, nodes, &left_bc, &right_bc );

  /* ----- Set the initial guess ----- */
  // y = (48/pi)x - 2 (linear fit to boundary conditions)
  double gradient( 48.0 / M_PI );
	for (std::size_t j=0; j < N; ++j )
	{
		double x = nodes[ j ];
		ode.solution()( j, y )  = gradient * x - 2.0;
    ode.solution()( j, yd ) = gradient;
	}

  ode.tolerance() = 1e-8;
  ode.max_iterations() = 50;
  ode.set_monitor_det( false );

  // Adaptively solve the BVP
  //ode.adapt_until( 1e-7 );

  cout << " Number of nodes in the solution = " << ode.solution().nnodes()
       << endl;

  nodes = ode.solution().nodes();
  Mesh1D<double, double> exact_solution( nodes, 2 );
  for (std::size_t j=0; j < nodes.size(); ++j )
	{
		double x = nodes[ j ];
		exact_solution( j, y )  = - 2.0 * cos( 2 * x ) + 10.0 * sin( 2 * x );
    exact_solution( j, yd ) =   4.0 * sin( 2 * x ) + 20.0 * cos( 2 * x );
	}
  // Solve
  ode.solve_bvp();

  // Output the solution data
  ode.solution().output( "./DATA/ODE_BVP_numerical_solution.dat" );
  exact_solution.output( "./DATA/ODE_BVP_exact_solution.dat" );

  //TODO plotting instructions

  cout << "--- FINISHED ---" << endl;
}
