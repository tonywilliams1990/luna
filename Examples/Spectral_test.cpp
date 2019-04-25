/// \file  Spectral_test.cpp
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
  cout << "------------------- Spectral_test ------------------" << endl;

  Basis<double> basis_function;

  Vector<double> x;
  x.lobatto_grid( 7 );
  cout << " x (Lobatto) = " << x << endl;
  //x.chebyshev_grid( 7 );
  //cout << " x (Chebyshev) = " << x << endl;

  Chebyshev<double> cheby;

  cout << " cheby( 0.5, 3 ) = " << cheby( 0.5, 3 ) << endl;
  cout << " cheby( x, 3 ) = " << cheby( x, 3 ) << endl;
  cout << " gegenbauer( 0.0, 4, 2 ) = " << cheby.gegenbauer( 0.0, 4, 2 ) << endl;
  cout << " cheby( 1.0, 3, 1 ) = " << cheby( 1.0, 3, 1 ) << endl;
  cout << " cheby( x, 3, 1 ) = " << cheby( x, 3, 1 ) << endl;

  Vector<double> coefficients( 3, 1.0 );

  Spectral<double> spectral( coefficients, "Chebyshev" );

  cout << "spectral.get_basis() = " << spectral.get_basis() << endl;
  cout << "spectral( 0.2 ) = " << spectral( 0.2 ) << endl;
  cout << "spectral( x ) = " << spectral( x ) << endl;

  cout << "spectral( 0.2, 1 ) = " << spectral( 0.2, 1 ) << endl;
  cout << "spectral( x, 1 ) = " << spectral( x, 1 ) << endl;

  coefficients[ 0 ] = 2.0;
  spectral.set_coefficients( coefficients );
  spectral[ 1 ] = 3.14;
  spectral.push_front( 7.45 );

  cout << "spectral.get_coefficients() = " << spectral.get_coefficients() << endl;
  cout << "spectral[ 2 ] = " << spectral[ 2 ] << endl;

  cout << "--- FINISHED ---" << endl;
}
