/// \file  Spectral_ODE_Rational.cpp
/// \ingroup Examples


#include "Luna/Core"
#include "Luna/ODE"
#include "Luna/Spectral"

namespace Luna
{
  namespace Example
  {

  } // End of namespace Example
} // End of namespace Luna

using namespace std;
using namespace Luna;
using namespace Example;

int main()
{
  cout << "------------------- Spectral_ODE_Rational ------------------" << endl;

  // Solve BVP u_yy - u = - 2*e^{1-y}, u(0)=e, u(inf)=0.
  // exact solution = (y+1)*e^{1-y}

  int n( 5 );                       // Number of spectral coefficients
  double L( 5.0 );                  // Map parameter
  Vector<double> y;                 // Vector for collocation grid
  y.rational_semi_grid( n, L );
  //cout << " y (Rational Chebyshev Semi) = " << y << endl;
  y.push_back( 0.0 );
  cout << " y (Rational Chebyshev Semi) = " << y << endl;

  int N( 10000 );                       // Number of output points
  Vector<double> grid;                  // Output grid
  grid.linspace( 0.0, y[ 0 ], N );
  Mesh1D<double> solution( grid, 3 );   // Output mesh
  // Exact solution
  for ( std::size_t i = 0; i < N; i++ )
  {
    double yi = grid[ i ];
    solution( i, 0 ) = ( yi + 1.0 ) * exp( 1.0 - yi );
  }

  Matrix<double> M( n + 1, n + 1, 0.0 );
  Vector<double> F( n + 1 ), a( n + 1 );
  RationalSemi<double> rationalsemi( L );         // Rational Chebyshev basis

  // y = inf boundary u = 0
  for ( std::size_t j = 0; j < n + 1; j++ )
  {
    double yi = y[ 0 ];
    M( 0, j ) = rationalsemi( yi, j );
  }
  F[ 0 ] = 0.0;

  // Internal nodes
  for ( std::size_t i = 1; i < n; i++ )
  {
    double yi = y[ i ];
    for ( std::size_t j = 0; j < n + 1; j++ )
    {
      M( i, j )  = rationalsemi( yi, j, 2 );                     // u_yy
      M( i, j ) -= rationalsemi( yi, j );                        // -u
    }
    F[ i ] = - 2 * exp( 1.0 - yi );                              // = -2*e^{1-y}
  }

  // y = 0 boundary u = e
  for ( std::size_t j = 0; j < n + 1; j++ )
  {
    double yi = y[ n ];
    M( n, j ) = rationalsemi( yi, j );
  }
  F[ n ] = exp( 1.0 );

  cout << " M = " << M << endl;
  cout << " F = " << F << endl;

  // Solve the system for the spectral coefficients
  a = M.solve( F );

  cout << " a = " << a << endl;

  // Put spectral solution into output mesh
  for ( std::size_t i = 0; i < N; i++ )
  {
    double yi( grid[ i ] );
    for ( std::size_t k = 0; k < n + 1; k++ )
    {
      solution( i, 1 ) += a[ k ] * rationalsemi( yi, k );
    }
    solution( i, 2 ) = solution( i, 1 ) - solution( i, 0 ); // difference
  }


  solution.output( "./DATA/Spectral_ODE_Rational.dat" );

  //TODO describe how to show solution

  cout << "--- FINISHED ---" << endl;
}
