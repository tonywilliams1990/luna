/// \file ODE_BVP_System.cpp
/// \ingroup Examples
/// \ingroup Spectral
/// The system of equations \f[ u''(x) - u(x) + 2 v(x) = 2 \left(1 - 2 x e^{-x}
/// \right), \f] \f[ v''(x) + x^2 v(x) + u(x) = x^2 - e^{-x}, \f] is solved
/// subject to the boundary conditions \f$ u(0) = v(0) = 0 \f$ and \f$ u(x) \to
/// 0 \f$, \f$ v(x) \to 1 \f$ as \f$ x \to \infty \f$. The system is solved
/// using a rational Chebyshev spectral method such that \f[ u(x) \approx u_n(x)
/// = \sum_{i=0}^n a_i TL_{i}(x; L), \f] \f[ v(x) \approx v_n(x) = \sum_{i=0}^n
/// a_{n + 1 + i} TL_{i}(x; L). \f]
/// Here \f$ L \f$ is the map parameter associated with the rational Chebyshev
/// functions \f$ TL_i(y; L) \f$. The exact solution is given by \f$ u(x) = x^2
/// e^{-x} \f$ and \f$ v(x) = 1 - e^{-x} \f$.


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
  cout << "---------------------- ODE_BVP_System ---------------------" << endl;

  cout << " * The system of equations u'' - u + 2v = 2( 1 - 2xe^(-x) ), " <<
  endl << " * v'' + x^2 v + u = x^2 - e^(-x), is solved using a rational " <<
  endl << " * Chebyshev spectral method." << endl;

  int n( 30 );                      // Number of spectral coefficients
  double L( 5.0 );                  // Map parameter
  double tol( 1e-10 );              // Tolerance correction coefficients norm
  Vector<double> x;                 // Vector for collocation grid
  x.rational_semi_grid( n, L );
  x.push_back( 0.0 );

  int N( 1000 );                        // Number of output points
  Vector<double> grid;                  // Output grid
  grid.powspace( 0.0, x[ 0 ], N, 3 );   // Squash nodes near to the boundary
  Mesh1D<double> solution( grid, 4 );   // Output mesh


  RationalSemi<double> rationalsemi( L );

  // Setup the system of equations
  Matrix<double> M( 2 * n + 2, 2 * n + 2, 0.0 );
  Vector<double> F( 2 * n + 2 ), a( 2 * n + 2 );
  int row( 0 );
  int v_offset( n + 1 );    // n coefficients + 1 bc

  // Don't need the BCs at infinity as these are "natural" BC

  // u equation
  for ( std::size_t i = 0; i < n; i++ )
  {
    double xi = x[ i ];
    for ( std::size_t j = 0; j < n + 1; j++ )
    {
      // u''
      M( row, j )  = rationalsemi( xi, j, 2 );
      // - u
      M( row, j ) -= rationalsemi( xi, j );
      // + 2v
      M( row, v_offset + j ) = 2 * rationalsemi( xi, j );
    }
    // = 2( 1 - 2 x e^(-x) )
    F[ row ] = 2 * ( 1. - 2 * xi * exp( - xi ) );
    ++row;
  }

  // x = 0 boundary u = 0
  for ( std::size_t j = 0; j < n + 1; j++ )
  {
    double xi = 0.0;
    M( row, j ) = rationalsemi( xi, j );
  }
  F[ row ] = 0.0;
  ++row;

  // v equation
  for ( std::size_t i = 0; i < n; i++ )
  {
    double xi = x[ i ];
    for ( std::size_t j = 0; j < n + 1; j++ )
    {
      // v''
      M( row, v_offset + j )  = rationalsemi( xi, j, 2 );
      // + x^2 v
      M( row, v_offset + j ) += xi * xi * rationalsemi( xi, j );
      // + u
      M( row, j ) = rationalsemi( xi, j );
    }
    // = x^2 - e^(-x)
    F[ row ]  = xi * xi - exp( -xi );
    ++row;
  }

  // x = 0 boundary v = 0
  for ( std::size_t j = 0; j < n + 1; j++ )
  {
    double xi = 0.0;
    M( row, v_offset + j ) = rationalsemi( xi, j );
  }
  F[ row ] = 0.0;
  ++row;

  // Solve the system for the correction spectral coefficients
  a = M.solve( F );

  // Split a into u and v parts and make Spectral solution for each
  Vector<double> a_u, a_v;
  for ( std::size_t i = 0; i < n + 1; i++ )
  {
    a_u.push_back( a[ i ] );
    a_v.push_back( a[ n + 1 + i ] );
  }

  Spectral<double> u( a_u, "RationalSemi", L );
  Spectral<double> v( a_v, "RationalSemi", L );

  // Put spectral solution into output mesh
  for ( std::size_t i = 0; i < N; i++ )
  {
    double x( grid[ i ] );
    solution( i, 0 ) = u( x );                        // u
    solution( i, 1 ) = v( x );                        // v
    solution( i, 2 ) = x * x * exp( - x );            // u_exact
    solution( i, 3 ) = 1.0 - exp( - x );              // v_exact
  }

  solution.output( "./DATA/ODE_BVP_System.dat" );
  cout << " * To see the spectral solution and the exact solution run: " << endl
       << "python Plotting/ODE_BVP_System_plot.py" << endl;

  cout << "--- FINISHED ---" << endl;
}
