/// \file Root_finding.cpp
/// \ingroup Examples
/// Examples of the numerical solution of polynomials

#include "Luna/Core"

using namespace std;
using namespace Luna;

namespace Luna
{
  class test_residual : public Residual<double>
  {
    public:
			// The test equation is 2nd order
			test_residual() : Residual<double> ( 2 ) {}

      // Define the equation
			void residual_fn( const Vector<double>& x, Vector<double>& F ) const
			{
				F[ 0 ] = std::pow( x[ 0 ], 3.0 ) + x[ 1 ] - 1;
				F[ 1 ] = std::pow( x[ 1 ], 3.0 ) - x[ 0 ] + 1;
				/*
					x^3 + y - 1 = 0,
					y^3 - x + 1 = 0,
					(x,y) = (1,0) is the only (real) solution
				*/
			}
  }; // End of class test_residual

  // Define the residual for arc-length continuation of a circle
  class Arc_problem : public Residual<double>
  {
    public:
      double p;

      Arc_problem() : Residual<double>( 1 ) {}

      void residual_fn( const Vector<double>& x, Vector<double>& F ) const
      {
        F[ 0 ] = x[ 0 ] * x[ 0 ] + p * p - 2.0;
        /*
          x^2 + p^2 - 2 = 0
        */
      }

  }; // End of class Arc_problem
}

int main()
{
  cout << "----------------- Root finding -----------------" << endl;

  Luna::test_residual Res;                           // Create residual
  Vector<double> x_0(2, 0.5);                        // Initial guess
  x_0[1] = 0.25;
  cout << " * x_0 = " << x_0 << endl;
  Luna::Newton<double> newton( &Res );               // Create Newton object
  newton.solve( x_0 );                               // Solve the system

  cout << " * Solution = " << x_0 << endl;

  cout << endl << "----- Arc continuation of x^2 + p^2 = 2.0 -----" << endl;

  // Instantiate the problem
  Luna::Arc_problem arc_prob;
  arc_prob.p = 1.0;

  const double tol = 1.e-6;
  // Scalar newton iteration problem
  Luna::Newton<double> newton_arc( &arc_prob, 8, tol );

  // Initial guess
  Vector<double> state( 1, 1.0 );
  double arc_step( 0.1 );
  double max_arc_step( 0.5 );
  newton_arc.init_arc( state, &arc_prob.p, arc_step, max_arc_step );

  do
  {
    newton_arc.arclength_solve( state );
    cout << "x = " << state[ 0 ] << ", p = " << arc_prob.p << endl;
  }while( state[ 0 ] < 1.0 );



  cout << "--- FINISHED ---" << endl;

}
