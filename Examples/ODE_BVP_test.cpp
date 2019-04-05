/// \file  ODE_BVP_test.cpp
/// \ingroup Examples
/// TODO description / rename file

#include "Luna/Core"
#include "Luna/ODE"

// ODE enumeration
enum{ y, yd };

namespace Luna
{
  class test_equation : public Equation<double>
	{

		public:
			// The test equation is 2nd order ( y'' + 4y = 0 )
			test_equation() : Equation<double> ( 2 ) {}

			// Define the equation
			void residual_fn( const Vector<double>& u, Vector<double>& F  ) const
			{
        double x( coord( 0 ) );
				F[ y ]   = u[ yd ];
				F[ yd ]  = - 4 * u[ y ];
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
	std::size_t N( 20 );              // Number of nodes
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
		double x = nodes[ j ];					// eta value at node j
		ode.solution()( j , y )  = gradient * x - 2.0;
    ode.solution()( j , yd ) = gradient;
	}

  // Solve
  ode.solve_bvp();

  // Output the solution data
  ode.solution().output( "./DATA/ODE_BVP_test.dat" );

  //TODO compare to exact solution y'(0) = 20

  cout << "--- FINISHED ---" << endl;
}
