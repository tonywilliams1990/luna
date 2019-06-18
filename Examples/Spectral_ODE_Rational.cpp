/// \file  Spectral_ODE_Rational.cpp
/// \ingroup Examples
/// The equation \f[ u''(y) - u(y) = -2e^{1-y} \f] is solved on
/// \f$ y \in [0,\infty)\f$ subject to the boundary conditions \f$ u(0) = e\f$
/// and \f$ u(y) \to 0 \f$ as \f$ y \to \infty\f$. The system is solved using
/// rational Chebyshev functions on a semi-infinite interval \f$ TL_i(y; L) \f$.
/// The solution \f$ u(y) \f$ is approximated by
/// \f[ u(y) \approx u_n(y) = \sum_{i=0}^n a_i TL_{i}(y; L) .\f]
/// Here \f$ n \f$ is the number of spectral coefficients \f$ a_i \f$ and
/// \f$ L \f$ is the map parameter associated with the rational Chebyshev
/// functions \f$ TL_i(y; L) \f$. The exact solution is given by
/// \f[ u(y) = (y+1) e^{1-y}. \f]
/// The spectral solution is found on \f$ y \in [0,y_{\infty}] \f$ where
/// \f$ y_{\infty} \f$ is determined by the number of spectral coefficients. For
/// more information on rational Chebyshev functions on a semi-infinite interval
/// see <b> Boyd, J.P., 1987. Orthogonal rational functions on a semi-infinite
/// interval. Journal of Computational Physics, 70(1), pp.63-88 </b>.


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

  cout << " * The equation u'' - u = - 2 * e^(1-y) is solved " << endl
       << " * subject to the boundary conditions u(0) = e and " << endl
       << " * u(y) -> 0 as y -> infinity." << endl;

  int n( 10 );                       // Number of spectral coefficients
  double L( 5.0 );                   // Map parameter
  Vector<double> y;                  // Vector for collocation grid
  y.rational_semi_grid( n, L );
  y.push_back( 0.0 );                // Extra point for imposing the zero BC

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

  // Create the pseudospectrally discretised matrix problem
  Matrix<double> M( n + 1, n + 1, 0.0 );
  Vector<double> F( n + 1 ), a( n + 1 );
  RationalSemi<double> rationalsemi( L );         // Rational Chebyshev basis

  // Don't need the BC at infinity as this is a "natural" BC for details see
  // Boyd (1987) Orthogonal rational functions on a semi-infinite interval.

  // Internal nodes
  for ( std::size_t i = 0; i < n; i++ )
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

  // Solve the system for the spectral coefficients
  a = M.solve( F );

  // Put spectral solution into the output mesh
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

  cout << " * For a comparison of the spectral/exact solutions run: " << endl;
  cout << "python Plotting/Spectral_ODE_Rational_plot.py" << endl;

  cout << "--- FINISHED ---" << endl;
}
