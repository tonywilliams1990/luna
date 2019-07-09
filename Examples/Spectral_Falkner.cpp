/// \file  Spectral_Falkner.cpp
/// \ingroup Examples

/// \todo TODO Problems with f' -> 1 at infinity as rationsl chebyshev functions
/// all tend to zero at infinity.


#include "Luna/Core"
#include "Luna/ODE"
#include "Luna/Spectral"

namespace Luna
{
  namespace Example
  {
    double initial_guess( const double& y )
    {
      //return y * ( 1.0 - exp( - y ) );
      //return ( y + 1.0 ) * exp( 1.0 - y );
      return - y * exp( - y );
    }
  } // End of namespace Example
} // End of namespace Luna

using namespace std;
using namespace Luna;
using namespace Example;

int main()
{
  cout << "--------------------- Spectral_Falkner --------------------" << endl;

  double beta( 0.0 );               // Hartree parameter
  int n( 25 );                      // Number of spectral coefficients
  int bcs( 2 );                     // Number of boundary conditions
  double L( 5.0 );                  // Map parameter
  double tol( 1e-10 );              // Tolerance correction coefficients norm
  Vector<double> y;                 // Vector for collocation grid
  y.rational_semi_grid( n, L );
  y.push_back( 0.0 );

  int N( 1000 );                        // Number of output points
  Vector<double> grid;                  // Output grid
  grid.powspace( 0.0, y[ 0 ], N, 2 );   // Squash nodes near to the boundary
  Mesh1D<double> solution( grid, 3 );   // Output mesh


  RationalSemi<double> rationalsemi( L );
  Vector<double> c( n , 0.0 );                 // Vector of coefficients

  // Approximate the guess as a semi-infinite rational Chebyshev polynomial
  c = rationalsemi.approximate( Example::initial_guess, n  );
  for ( std::size_t i = 0; i < bcs; i++ )
  {
    c.push_back( 0.0 );
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
    // Setup the system of equations to find the correction u_c
    Matrix<double> M( n + bcs, n + bcs, 0.0 );
    Vector<double> F( n + bcs ), a_c( n + bcs );

    // Don't need the BC at infinity as this is a "natural" BC

    // Internal nodes
    for ( std::size_t i = 0; i < n; i++ )
    {
      double yi = y[ i ];
      for ( std::size_t j = 0; j < n + bcs; j++ )
      {
        // u_c'''
        M( i, j )  = rationalsemi( yi, j, 3 );
        // (y + u_g) * u_c''
        M( i, j ) += ( yi + u_g( yi ) ) * rationalsemi( yi, j, 2 );
        // - 2 * beta * (1 + u_g') * u_c'
        M( i, j ) -= 2 * beta * ( 1 + u_g( yi, 1 ) ) * rationalsemi( yi, j, 1 );
        // u_g'' * u_c
        M( i, j ) += u_g( yi, 2 ) * rationalsemi( yi, j );
      }
      // - u_g''' - (y + u_g) * u_g'' + beta * u_g' ( 2 + u_g' )
      F[ i ]  = - u_g( yi, 3 ) - ( yi + u_g( yi ) ) * u_g( yi, 2 );
      F[ i ] += beta * u_g( yi, 1 ) * ( 2 + u_g( yi, 1 ) );
    }

    // y = 0 boundary f = 0
    for ( std::size_t j = 0; j < n + bcs; j++ )
    {
      double yi = 0.0;
      M( n, j ) = rationalsemi( yi, j );
    }
    F[ n ] = - u_g( 0.0 );

    // y = 0 boundary f' = -1
    for ( std::size_t j = 0; j < n + bcs; j++ )
    {
      double yi = 0.0;
      M( n + 1, j ) = rationalsemi( yi, j, 1 );
    }
    F[ n + 1 ] = - 1.0 - u_g( 0.0, 1 );

    // Solve the system for the correction spectral coefficients
    a_c = M.solve( F );
    u_g.update_coefficients( a_c );
    norm = a_c.norm_2();

    cout << " * norm = " << std::scientific << norm << endl;
    timer.print();

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

  cout << "--- FINISHED ---" << endl;
}
