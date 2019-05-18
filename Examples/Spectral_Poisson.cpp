/// \file  Spectral_Poisson.cpp
/// \ingroup Examples
/// The equation \f[ \nabla^2 u(x,y) = 2(x^2 + y^2) \f] is solved subject to the
/// boundary conditions \f$ u(x,\pm 1) = x^2 \f$ \f$ u(\pm 1, y ) = y^2. \f$
/// The system is solved using a Chebyshev spectral method.
/// TODO full description of the problem and method of solution

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
  cout << "------------------- Spectral_Poisson ------------------" << endl;

  int I( 3 );           // Number of x collocation points
  int J( 4 );           // Number of y collocation points
  int size( I * J );

  Vector<double> x, y;
  x.lobatto_grid( I );
  y.lobatto_grid( J );

  cout << " x = " << x << endl;
  cout << " y = " << y << endl;

  Matrix<double> L( size, size, 0.0 );
  Vector<double> F( size ), a( size );

  double xi, yj;

  // x = -1 boundary u = y^2 (left)
  xi = x[ I - 1 ];
  cout << xi << endl;


  //loop over internal nodes


  // y = -1 boundary u = x^2 (bottom)

  // Internal nodes

  // y = 1 boundary u = x^2 (top)


  // x = 1 boundary u = y^2 (right)


  // Solve the system for the spectral coefficients

  cout << " L = " << L << endl;
  cout << " F = " << F << endl;
  cout << " a = " << a << endl;

  cout << "--- FINISHED ---" << endl;
}
