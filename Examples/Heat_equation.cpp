/// \file  Heat_equation.cpp
/// \ingroup Examples
/// TODO Description + tex equations

#include "Luna/Core"
#include "Luna/Sparse"

using namespace std;
using namespace Luna;

int main()
{
  cout << "----- Heat equation -----" << endl;

  // Diagonals of the matrix
  Vector<double> a( 3 );
  a[ 0 ] = - 1.; a[ 1 ] = - 1.; a[ 2 ] = - 2.;

  Vector<double> b( 4, 3. );

  Vector<double> c( 3 );
  c[ 0 ] = - 2.; c[ 1 ] = - 1.; c[ 2 ] = - 1.;

  Tridiagonal<double> tridiag( a, b, c );

  //cout << " * a = " << tridiag.get_sub() << endl;
  //cout << " * b = " << tridiag.get_main() << endl;
  //cout << " * c = " << tridiag.get_super() << endl;
  cout << " * tridiag = " << tridiag << endl;
  cout << " * tridiag( 1, 0 ) = " << tridiag( 1, 0 ) << endl;

  // RHS vector
  Vector<double> d( 4 );
  d.linspace( 1., 4., 4 );
  cout << " * d = " << d << endl;

  // Solution vector
  Vector<double> x( 4 );
  cout << " * x = " << x << endl;

  cout << "--- FINISHED ---" << endl;

}
