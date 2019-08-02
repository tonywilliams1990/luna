/// \file SVD_test.cpp
/// \ingroup Examples
/// TODO description / rename file

#include "Luna/Core"
#include "Luna/Sparse"

using namespace std;
using namespace Luna;

int main()
{
  cout << "--------------------- SVD_test --------------------" << endl;

  Matrix<double> A( 3, 2, 0.0 );
  A( 0, 0 ) = 1.0; A( 0, 1 ) = 3.0;
  A( 1, 0 ) = 2.0; A( 1, 1 ) = 4.0;
  A( 2, 0 ) = 1.0; A( 2, 1 ) = 6.0;

  cout << " A = " << A << endl;

  SVD<double> svd( A );

  cout << " w = " << svd.singular_values() << endl;
  cout << " nullspace = " << svd.nullspace() << endl;
  cout << " u = " << svd.left_vectors() << endl;
  cout << " v = " << svd.right_vectors() << endl;

  Vector<double> b( 3 ), x( 2 );
  b[ 0 ] = 4.0; b[ 1 ] = 1.0; b[ 2 ] = 3.0;
  cout << " b = " << b << endl;
  svd.solve( b, x );
  cout << " x = " << x << endl;


  /// \todo TODO test complex matrix too

  cout << "--- FINISHED ---" << endl;
}
