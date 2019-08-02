/// \file Falkner.cpp
/// \ingroup Examples
/// \ingroup Spectral
/// The Falkner-Skan equation \f[ f'''(\eta) + f(\eta)f''(\eta) + \beta \left[1
/// - \left(f'(\eta) \right)^2 \right] = 0 \f] is solved on \f$ \eta \in
/// [0,\infty) \f$ subject to the boundary conditions \f$ f(0) = f'(0) =  0\f$
/// and \f$ f'(\eta) \to 1 \f$ as \f$ \eta \to \infty\f$.
/// The system is solved using rational Chebyshev functions on a semi-infinite
/// interval \f$ TL_i(\eta; L) \f$. The rational chebyshev functions all have zero
/// derivative as \f$ \eta \to \infty \f$ so it is necessary to use the change
/// of variables \f$ u(\eta) = f(\eta) - \eta \f$. The Falkner-Skan system is
/// transformed to \f[ u''' + (\eta + u)u'' - \beta u'\left(2 + u' \right) = 0
/// \f] subject to the boundary conditions \f$ u(0) = 0\f$, \f$ u'(0) = -1 \f$
/// and \f$ u' \to 0 \f$ as \f$ \eta \to \infty \f$. This nonlinear equation is
/// solver using Newton iteration on the spectral coefficients.

#include "Luna/Core"
#include "Luna/ODE"
#include "Luna/Spectral"

namespace Luna
{
  namespace Example
  {
    double initial_guess( const double& y )
    {
      return - y * exp( - y );
    }
  } // End of namespace Example
} // End of namespace Luna

using namespace std;
using namespace Luna;
using namespace Example;

int main()
{
  cout << "------------------------- Falkner -------------------------" << endl;

  cout << " * The Falkner-Skan equation, for a boundary layer on flat " << endl
       << " * plate, is solved using a rational Chebyshev spectral " << endl
       << " * method." << endl;

  double beta( 0.0 );               // Hartree parameter
  int n( 40 );                      // Number of spectral coefficients
  double L( 6.0 );                  // Map parameter
  double tol( 1e-10 );              // Tolerance correction coefficients norm
  Vector<double> y;                 // Vector for collocation grid
  y.rational_semi_grid( n, L );
  y.push_back( 0.0 );

  int N( 1000 );                        // Number of output points
  Vector<double> grid;                  // Output grid
  grid.powspace( 0.0, y[ 0 ], N, 2 );   // Squash nodes near to the boundary
  Mesh1D<double> solution( grid, 3 );   // Output mesh


  RationalSemi<double> rationalsemi( L );
  Vector<double> c( n, 0.0 );                 // Vector of coefficients

  // Approximate the guess as a semi-infinite rational Chebyshev polynomial
  c = rationalsemi.approximate( Example::initial_guess, n );
  for ( std::size_t i = 0; i < 2; i++ )
  {
    c.push_back( 0.0 ); // Extra coefficients for BCs
  }
  // Make this into a spectral solution
  Spectral<double> u_g( c, "RationalSemi", L );

  double norm;
  int max_iter( 100 );
  int iter( 0 );

  Timer timer;
  timer.start();

  do
  {
    // Setup the system of equations to find the correction
    Matrix<double> M( n + 2, n + 2, 0.0 );
    Vector<double> F( n + 2 ), a_c( n + 2 );

    // Don't need the BC at infinity as this is a "natural" BC

    // Internal nodes
    for ( std::size_t i = 0; i < n; i++ )
    {
      double yi = y[ i ];
      for ( std::size_t j = 0; j < n + 2; j++ )
      {
        // u_c'''
        M( i, j )  = rationalsemi( yi, j, 3 );
        // ( y + u_g ) * u_c''
        M( i, j ) += ( yi + u_g( yi ) ) * rationalsemi( yi, j, 2 );
        // - 2 * beta * ( 1 + u_g' ) * u_c'
        M( i, j ) -= 2 * beta * ( 1 + u_g( yi, 1 ) ) * rationalsemi( yi, j, 1 );
        // u_g'' * u_c
        M( i, j ) += u_g( yi, 2 ) * rationalsemi( yi, j );
      }
      // - u_g''' - ( y + u_g ) * u_g'' + beta * u_g' ( 2 + u_g' )
      F[ i ]  = - u_g( yi, 3 ) - ( yi + u_g( yi ) ) * u_g( yi, 2 );
      F[ i ] += beta * u_g( yi, 1 ) * ( 2 + u_g( yi, 1 ) );
    }

    // y = 0 boundary u_c = - u_g
    for ( std::size_t j = 0; j < n + 2; j++ )
    {
      double yi = 0.0;
      M( n, j ) = rationalsemi( yi, j );
    }
    F[ n ] = - u_g( 0.0 );

    // y = 0 boundary u_c' = -1 - u_g'
    for ( std::size_t j = 0; j < n + 2; j++ )
    {
      double yi = 0.0;
      M( n + 1, j ) = rationalsemi( yi, j, 1 );
    }
    F[ n + 1 ] = - 1.0 - u_g( 0.0, 1 );

    // Solve the system for the correction spectral coefficients
    a_c = M.solve( F );
    u_g.update_coefficients( a_c );
    norm = a_c.norm_2();
    ++iter;

  }while( norm > tol && iter < max_iter );

  cout << " * The final spectral coefficient for n = " << n << " is: "
       << u_g.get_coefficients()[ n-1 ] << endl;
  cout << " * f''(0) = " << u_g( 0.0, 2 ) << endl;


  // Put spectral solution into output mesh
  for ( std::size_t i = 0; i < N; i++ )
  {
    double y( grid[ i ] );
    solution( i, 0 ) = y + u_g( y );              // f   = eta + u
    solution( i, 1 ) = 1 + u_g( y, 1 );           // f'  = 1 + u'
    solution( i, 2 ) = u_g( y, 2 );               // f'' = u''
  }

  solution.output( "./DATA/Spectral_Falkner.dat" );

  timer.print();
  timer.stop();

  cout << " * To see the spectral solution and its derivatives run: " << endl;
  cout << "python Plotting/Spectral_Falkner_plot.py" << endl;

  cout << "--- FINISHED ---" << endl;
}
