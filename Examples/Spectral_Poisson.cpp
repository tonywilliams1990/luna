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
  Chebyshev<double> cheby;

  double xi, yj;
  int f, g, i, j;


  for ( int M = 0; M < size; M++ )
  {
    i = M / J;
    j = M % J;
    //cout << " i, j = " << i << ", " << j;

    xi = x[ I - ( i + 1 ) ];
    yj = y[ J - ( j + 1 ) ];

    //cout << "\t x, y = " << xi << ", " << yj << endl;

    // x = -1 boundary u = y^2 (left)
    if ( i == 0 )
    {
      for ( int N = 0; N < size; N++ )
      {
        f = N / J;
        g = N % J;
        L( M, N ) = cheby( -1.0, f ) * cheby( yj, g );
      }
      F[ M ] = yj * yj;
    }

    // y = -1 boundary u = x^2 (bottom)
    if ( j == 0 && i != 0 && i != I - 1 )
    {
      for ( int N = 0; N < size; N++ )
      {
        f = N / J;
        g = N % J;
        L( M, N ) = cheby( xi, f ) * cheby( -1.0, g );
      }
      F[ M ] = xi * xi;
    }

    // Internal nodes
    if ( i != 0 && i != I - 1 && j != 0 && j != J - 1 )
    {
      for ( int N = 0; N < size; N++ )
      {
        f = N / J;
        g = N % J;
        L( M, N )  = cheby( xi, f, 2 ) * cheby( yj, g ); // u_xx
        L( M, N ) += cheby( xi, f ) * cheby( yj, g, 2 ); // u_yy
      }
      F[ M ] = 2 * ( xi * xi + yj * yj );
    }

    // y = 1 boundary u = x^2 (top)
    if ( j == J - 1 && i != 0 && i != I - 1 )
    {
      for ( int N = 0; N < size; N++ )
      {
        f = N / J;
        g = N % J;
        L( M, N ) = cheby( xi, f ) * cheby( 1.0, g );
      }
      F[ M ] = xi * xi;
    }

    // x = 1 boundary u = y^2 (right)
    if ( i == I - 1 )
    {
      for ( int N = 0; N < size; N++ )
      {
        f = N / J;
        g = N % J;
        L( M, N ) = cheby( 1.0, f ) * cheby( yj, g );
      }
      F[ M ] = yj * yj;
    }

  }

  //TODO 1 inner N loop with if statements

  // Solve the system for the spectral coefficients
  a = L.solve( F );

  cout << " L = " << L << endl;
  cout << " F = " << F << endl;
  cout << " a = " << a << endl;

  // TODO output to a 2D mesh + compare with exact u = x^2 * y^2
  // Solution is very symmetrical = lots of 0 coefficients 

  cout << "--- FINISHED ---" << endl;
}
