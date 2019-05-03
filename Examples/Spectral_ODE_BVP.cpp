/// \file  Spectral_ODE_BVP.cpp
/// \ingroup Examples
/// The equation \f[ u''(x) - (x^6 + 3x^2)u(x) = 0 \f] is solved subject to the
/// boundary conditions \f$ u(-1) = u(1) = 1. \f$ The system is solved using a
/// Chebyshev spectral method. The solution \f$ u(x) \f$ is symmetric so only
/// even Chebyshev polynomials are needed to approximate the solution such that
/// \f[ u(x) \approx u_N(x) = \sum_{i=0}^N a_i T_{2i}(x) .\f] The exact
/// solution is given by \f[ u(x) = \exp \left[ \frac{x^4 - 1}{4} \right]. \f]
/// The spectral solution is found on \f$ x \in [0,1] \f$ and then because the
/// solution is symmetric this gives the solution in the entire domain
/// \f$ x \in [-1,1]. \f$

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
  cout << "------------------- Spectral_ODE_BVP ------------------" << endl;

  cout << " * The equation u'' + (x^6 + 3x^2)u = 0 is solved " << endl
       << " * subject to the boundary conditions u(-1) = u(1) = 0." << endl
       << " * Only symmetric even Chebyshev polynomials are used in" << endl
       << " * the spectral solution. " << endl;

  int n( 1 );           // Number of spectral coefficients
  int N( 100 );         // Number of output points
  int num_out( 3 );     // Number of different n values
  int step( 2 );

  Vector<int> n_values;

  Vector<double> grid;
  grid.linspace( -1.0, 1.0, N );
  Mesh1D<double> solution( grid, 1 + num_out );

  Vector<double> exact( N ), numerical( N );

  for ( std::size_t i = 0; i < N; i++ )
  {
    solution( i, 0 ) = std::exp( 0.25 * ( std::pow( grid[ i ], 4 ) - 1.0 ) );
    exact[ i ] = solution( i, 0 );
  }

  int counter( 0 );

  do
  {

    Vector<double> x;
    x.half_lobatto_grid( n );

    Matrix<double> L( n, n, 0.0 );
    Vector<double> F( n ), a( n );
    Chebyshev<double> cheby;

    // x = 1 boundary u(1) = 1
    for ( std::size_t j = 0; j < n; j++ )
    {
      L( 0, j ) = cheby( 1.0, 2 * j );
    }
    F[ 0 ] = 1.0;

    // Internal nodes and x = 0 symmetric boundary
    for ( std::size_t i = 1; i < n; i++ )
    {
      double xi = x[ i ];
      double xi2 = xi * xi;
      double xi6 = std::pow( xi, 6 );
      for ( std::size_t j = 0; j < n; j ++ )
      {
        L( i, j ) = cheby( xi, 2 * j, 2 );                     // u_xx
        L( i, j ) += - ( xi6 + 3 * xi2 ) * cheby( xi, 2 * j ); // (x^6 + 3x^2)u
      }
      F[ i ] = 0.0;
    }

    // Solve the system for the spectral coefficients
    a = L.solve( F );

    Spectral<double> spectral( a, "EvenChebyshev" );
    numerical = spectral( grid );

    Vector<double> diff;
    diff = exact - numerical;
    cout << " * ||u_exact - u_spectral|| ( n = " << n << " ) = "
         << std::scientific << diff.norm_2() << endl;

    for ( std::size_t i = 0; i < N; i++ )
    {
      solution( i, 1 + counter ) = spectral( grid[ i ] );
    }
    n_values.push_back( n );
    n += step;
    counter++;

  }while( counter < num_out );

  solution.output( "./DATA/Spectral_ODE_BVP.dat" );

  cout << " * For a comparison of the spectral/exact solutions run: " << endl;
  cout << "python Plotting/Spectral_ODE_BVP_plot.py --n";
  for ( std::size_t i = 0; i < n_values.size(); i++ )
  {
    cout << " " << n_values[ i ];
  }
  cout << endl;

  cout << "--- FINISHED ---" << endl;
}
