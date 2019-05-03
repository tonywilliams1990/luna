/// \file  Spectral_ODE_Nonlinear.cpp
/// \ingroup Examples
/// TODO \todo description of the problem

#include "Luna/Core"
#include "Luna/ODE"
#include "Luna/Spectral"

namespace Luna
{
  namespace Example
  {
    double initial_guess( const double& x )
    {
      double e( std::exp( 1.0 ) );
      return 0.5 * ( e - ( 1 / e ) ) * x;
    }
  } // End of namespace Example
} // End of namespace Luna

using namespace std;
using namespace Luna;
using namespace Example;

int main()
{
  cout << "------------------ Spectral_ODE_Nonlinear -----------------" << endl;

  int n( 8 );                          // Number of spectral coefficients
  Vector<double> x;
  x.lobatto_grid( n );
  cout << " x (Chebyshev) = " << x << endl;

  // Output grid
  int N( 100 );         // Number of output points
  Vector<double> grid;
  grid.linspace( -1.0, 1.0, N );
  Mesh1D<double> solution( grid, 3 );

  // Exact solution
  for ( std::size_t i = 0; i < N; i++ )
  {
    solution( i, 0 ) = std::sinh( grid[ i ] );
  }

  Chebyshev<double> cheby;

  Vector<double> c( n, 0.0 );                 // Vector of coefficients
  Spectral<double> u_c( c, "Chebyshev" );  // correction spectral solution

  // Approximate the initial guess as a Chebyshev polynomial
  c = cheby.approximate( Example::initial_guess, n );
  // Make this into a spectral solution
  Spectral<double> u_g( c, "Chebyshev" );
  double norm;
  double m( 0.5 * ( std::exp( 1.0 ) - std::exp( - 1.0 ) ) );

  do
  {

    cout << " u_c = " << u_c.get_coefficients() << endl;
    cout << " u_g = " << u_g.get_coefficients() << endl;

    // Setup the system of equations to find the correction u_c
    Matrix<double> L( n, n, 0.0 );
    Vector<double> F( n ), a( n );

    // x = 1 boundary u(1) = 0.5 * ( e - e^(-1) )
    for ( std::size_t j = 0; j < n; j++ )
    {
      L( 0, j ) = cheby( 1.0, j );
    }
    F[ 0 ] = 0.5 * ( std::exp( 1.0 ) - std::exp( - 1.0 ) ) - u_g( 1.0 );

    // Internal nodes
    for ( std::size_t i = 1; i < n - 1; i++ )
    {
      double xi = x[ i ];
      for ( std::size_t j = 0; j < n; j ++ )
      {
        L( i, j ) = u_g( xi ) * cheby( xi, j, 2 );          // u_g * u_c''
        L( i, j ) += - 2.0 * u_g( xi, 1 ) * cheby( xi, j, 1 ); // -2 * u_g' * u_c'
        L( i, j ) += u_g( xi, 2 ) * cheby( xi, j );         // u_g'' * u_c
      }
      // u_g' * u_g' - u_g * u_g'' - 1
      F[ i ] = u_g( xi, 1 ) * u_g( xi, 1 ) - u_g( xi ) * u_g( xi, 2 ) - 1.0;
    }

    // x = -1 boundary u(-1) = 0.5 * ( e^(-1) - e )
    for ( std::size_t j = 0; j < n; j++ )
    {
      L( n - 1, j ) = cheby( - 1.0, j );
    }
    F[ n - 1 ] = 0.5 * ( std::exp( - 1.0 ) - std::exp( 1.0 ) ) - u_g( - 1.0 );


    //cout << " L = " << L << endl;
    //cout << " F = " << F << endl;

    // Solve the system for the spectral coefficients
    a = L.solve( F );

    //cout << " a = " << a << endl;

    u_c.set_coefficients( a );
    u_g.update_coefficients( a );

    norm = a.norm_2();
    cout << " a.norm_2 = " << norm << endl;

  }while( norm > 1e-8 );

  cout << " u_c = " << u_c.get_coefficients() << endl;
  cout << " u_g = " << u_g.get_coefficients() << endl;

  // Put spectral solution into output mesh
  for ( std::size_t i = 0; i < N; i++ )
  {
    solution( i, 1 ) = u_g( grid[ i ] );
    solution( i, 2 ) = solution( i, 1 ) - solution( i, 0 ); // difference
  }

  solution.output( "./DATA/Spectral_ODE_Nonlinear.dat" );

  cout << "--- FINISHED ---" << endl;
}
