/// \file  Spectral_ODE_BVP.cpp
/// \ingroup Examples
/// TODO full description including equations


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

  //TODO introduce the program - actually should only use even Chebyshev functions
  // as the solution is symmetric

  int n( 2 );         // Number of spectral coefficients
  int N( 100 );        // Number of output points
  int num_out( 4 );     // Number of different n values
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
    x.lobatto_grid( n );

    Matrix<double> L( n, n, 0.0 );
    Vector<double> F( n ), a( n );
    Chebyshev<double> cheby;

    // x = -1 boundary u(-1) = 1
    for ( std::size_t j = 0; j < n; j++ )
    {
      L( 0, j ) = cheby( - 1.0, j );
    }
    F[ 0 ] = 1.0;

    // Internal nodes
    for ( std::size_t i = 1; i < n - 1; i++ )
    {
      double xi = x[ i ];
      double xi2 = xi * xi;
      double xi6 = std::pow( xi, 6 );
      for ( std::size_t j = 0; j < n; j ++ )
      {
        L( i, j ) = cheby( xi, j, 2 );                       // u_xx
        L( i, j ) += - ( xi6 + 3 * xi2 ) * cheby( xi, j );   // (x^6 + 3x^2) * u
      }
      F[ i ] = 0.0;
    }

    // x = 1 boundary u(1) = 1
    for ( std::size_t j = 0; j < n; j++ )
    {
      L( n - 1, j ) = cheby( 1.0, j );
    }
    F[ n - 1 ] = 1.0;

    // Solve the system for the spectral coefficients
    a = L.solve( F );

    Spectral<double> spectral( a, "Chebyshev" );
    numerical = spectral( grid );

    Vector<double> diff;
    diff = exact - numerical;
    cout << " * ||u - u_spectral|| ( n = " << n << " ) = " << std::scientific
         << diff.norm_2() << endl;

    for ( std::size_t i = 0; i < N; i++ )
    {
      solution( i, 1 + counter ) = spectral( grid[ i ] );
    }
    n_values.push_back( n );

    n += step;
    counter++;

  }while( counter < num_out );

  solution.output( "./DATA/Spectral_ODE_BVP.dat" );

  cout << " * For a comparison of the spectral solutions run: " << endl;
  cout << "python Plotting/Spectral_ODE_BVP_plot.py --n";
  for ( std::size_t i = 0; i < n_values.size(); i++ )
  {
    cout << " " << n_values[ i ];
  }
  cout << endl;

  cout << "--- FINISHED ---" << endl;
}
