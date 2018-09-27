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

  A( 0, 1 ) = 2.0; A( 0, 2 ) = 3.0; A( 2, 1 ) = 7.0;
  B( 1, 1 ) = 5.0; B( 2, 2 ) = 8.0; B( 2, 1 ) = 8.0;

  cout << "  * A = " << A << endl;
  cout << "  * B = " << B << endl;

  cout << "  * A.rows() = " << A.rows() << endl;
  cout << "  * A.cols() = " << A.cols() << endl;
  cout << "  * A.numel() = " << A.numel() << endl;

  C = A * 4;
  cout << "  * C = " << C << endl;

  cout << "--- FINISHED ---" << endl;
}
