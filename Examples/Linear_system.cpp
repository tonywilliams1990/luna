/// \file  Linear_system.cpp
/// \ingroup Examples
/// Solve some linear systems of equations

#include "Luna/Core"

using namespace std;
using namespace Luna;

int main()
{
  cout << "----- Linear system -----" << endl;

  Matrix<double> A( 3, 3, 1.0 );
  Matrix<double> B( A );
  Matrix<double> C;

  A( 0, 1 ) = 2.0; A( 0, 2 ) = -3.0; A( 2, 1 ) = 7.0;
  B( 1, 1 ) = 5.5; B( 2, 2 ) = 8.0; B( 2, 1 ) = 8.0;

  cout << "  * A = " << A << endl;
  cout << "  * B = " << B << endl;

  cout << "  * A.rows() = " << A.rows() << endl;
  cout << "  * A.cols() = " << A.cols() << endl;
  cout << "  * A.numel() = " << A.numel() << endl;

  C = A * 2;
  C -= A;
  C -= 2;
  cout << "  * C = " << C << endl;

  Vector<double> row( 3, 3.14 );
  C.set_row( 2, row );
  cout << "  * C = " << C << endl;

  C = A * B;
  cout << "  * C = A * B = " << C << endl;

  Vector<double> vec;
  vec = C.get_col( 0 );

  vec = A * row;
  cout << "  * x = " << row << endl;
  cout << "  * b = A * x = " << vec << endl;

  //C.transpose_in_place();
  cout << "  * C^T = " << C.transpose() << endl;

  Matrix< std::complex<double> > D( 4, 3, 1.0 );
  //D( 1, 2 ) = -7.0;
  D.random();
  cout << "  * D = " << D << endl;
  D.transpose_in_place();
  cout << "  * D^T = " << D << endl;

  cout << "--- FINISHED ---" << endl;
}
