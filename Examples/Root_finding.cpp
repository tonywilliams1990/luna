/// \file Root_finding.cpp
/// \ingroup Examples
/// Examples of the numerical solution of equations using Newton's method and
/// arc-length continuation and the solution of polynomial equations.

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

  cout << " * Solving the set of equations " << endl
       << " * x^3 + y - 1 = 0, " << endl
       << " * y^3 - x + 1 = 0, " << endl
       << " * using Newton's method and the initial guess " << endl;

  Luna::test_residual Res;                           // Create residual
  Vector<double> x_0(2, 0.5);                        // Initial guess
  x_0[1] = 0.25;
  cout << " * x_0 = " << x_0 << endl;
  Luna::Newton<double> newton( &Res );               // Create Newton object
  newton.solve( x_0 );                               // Solve the system

  cout << " * gives the solution (x,y)= " << x_0 << endl;

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

  cout << endl << "--------- Polynomial equations ---------" << endl;
  cout << " * Define the function f(x) = -1 + 2x + x^2" << endl;
  Luna::Polynomial<double> f( { -1.0, 2.0, 1.0, 1.0 } );
  cout << " * Access the x coefficient f[1] = " << f[1] << endl;
  cout << " * Evaluate the polynomial at x = 0.1 " << endl
       << " * f( 0.1 ) = " << f( 0.1 ) << endl << endl;

  cout << " * Solve the complex quadratic equation " << endl
       << " * ( 1 + i )x^2 - ( 2 + i )x - 1 = 0 " << endl;
  Luna::Polynomial< std::complex<double> > p;
  std::complex<double> a ( 1.0, 1.0 );
  std::complex<double> b ( - 2.0, - 1.0 );
  std::complex<double> c ( - 1.0, 0.0 );
  Vector< std::complex<double> > roots;
  roots = p.quadratic_solve( a, b, c );
  cout << " * The two roots of the equation are" << endl;
  cout << " * x_0 = " << roots[ 0 ] << endl;
  cout << " * x_1 = " << roots[ 1 ] << endl << endl;

  cout << " * Solve the cubic equation" << endl
       << " * 2x^3 + x^2 + 2x - 1 = 0 " << endl;
  Vector<std::complex<double>> f_roots;
  f_roots = f.cubic_solve( 2.0, 1.0, 2.0, - 1.0 );
  cout << " * The three roots of the equation are" << endl;
  cout << " * x_0 = " << f_roots[ 0 ] << endl;
  cout << " * x_1 = " << f_roots[ 1 ] << endl;
  cout << " * x_2 = " << f_roots[ 2 ] << endl << endl;

  cout << " * Dividing the polynomial u(x) = x^3 - 2x^2 - 4" << endl
       << " * by v(x) = x - 3 gives the quotient " << endl;
  Polynomial<double> u( std::vector<double> { -4.0, 0.0, -2.0, 1.0 } );
  Polynomial<double> v( std::vector<double> { -3.0, 1.0 } );
  Polynomial<double> q, r;
  u.polydiv( v, q, r );
  cout << " * q(x) = " << q.coeffs()[0] << " + " << q.coeffs()[1] << "x + "
       << q.coeffs()[2] << "x^2" << endl
       << " * and the remainder " << endl
       << " * r(x) = " << r.coeffs() << endl
       << " * such that u(x) = q(x) * v(x) + r(x)." << endl << endl;

  cout << " * Find the roots of the quartic equation " << endl
       << " * 3 + x + x^2 + x^4 = 0 using Laguerre's method. " << endl;
  Polynomial<double> poly( std::vector<double> { 3.0, 1.0, 1.0, 0.0, 1.0 } );
  Vector<std::complex<double>> poly_roots;
  poly_roots = poly.roots( true );
  cout << " * The roots of the equation are" << endl
       << " * x_0 = " << poly_roots[0] << endl
       << " * x_1 = " << poly_roots[1] << endl
       << " * x_2 = " << poly_roots[2] << endl
       << " * x_3 = " << poly_roots[3] << endl << endl;

  cout << "--- FINISHED ---" << endl;

}
