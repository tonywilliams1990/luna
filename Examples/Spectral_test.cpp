/// \file  Spectral_test.cpp
/// \ingroup Examples


#include "Luna/Core"
#include "Luna/ODE"
#include "Luna/Spectral"

namespace Luna
{
  namespace Example
  {
    double function( const double& x )
    {
      return std::sinh( x );
    }
  } // End of namespace Example
} // End of namespace Luna

using namespace std;
using namespace Luna;
using namespace Example;

int main()
{
  cout << "------------------- Spectral_test ------------------" << endl;

  Basis<double> basis_function;

  int N( 4 );
  Vector<double> x;
  x.chebyshev_grid( N );
  cout << " x (Chebyshev) = " << x << endl;


  Chebyshev<double> cheby;

  /*cout << " cheby( 0.5, 3 ) = " << cheby( 0.5, 3 ) << endl;
  cout << " cheby( x, 3 ) = " << cheby( x, 3 ) << endl;
  cout << " gegenbauer( 0.0, 4, 2 ) = " << cheby.gegenbauer( 0.0, 4, 2 ) << endl;
  cout << " cheby( 1.0, 3, 1 ) = " << cheby( 1.0, 3, 1 ) << endl;
  cout << " cheby( x, 3, 1 ) = " << cheby( x, 3, 1 ) << endl;*/

  Vector<double> c( N, 0.0 );

  for ( std::size_t j = 0; j < N; j++ )
  {
    double sum( 0.0 );
    for ( std::size_t k = 0; k < N; k++ )
    {
      sum += Example::function( x[ k ] ) * std::cos( ( j * M_PI * ( k + 0.5 ) ) / N );
    }
    c[ j ] = ( 2.0 / N ) * sum;
  }

  cout << " c = " << c << endl;

  c = cheby.approximate_odd( Example::function, N / 2 );

  cout << " c = " << c << endl;

  Spectral<double> spectral( c, "Chebyshev" );



  cout << "--- FINISHED ---" << endl;
}
