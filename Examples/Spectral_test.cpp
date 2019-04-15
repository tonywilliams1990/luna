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

  cout << "--- FINISHED ---" << endl;
}
