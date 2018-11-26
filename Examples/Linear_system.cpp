/// \file  Linear_system.cpp
/// \ingroup Examples
/// Solve some linear systems of equations

#include "Luna/Core"

using namespace std;
using namespace Luna;

int main(int argc, char ** argv)
{
  cout << "---------------- Linear system ----------------" << endl;

  Matrix<double> A( 3, 3, 0.0 );
  Vector<double> b( 3, 0.0 );
  Vector<double> x;

  cout << "  * Solve the linear system Ax=b for x where" << endl;

  A( 0, 0 ) = 1.0; A( 0, 1 ) = 1.0; A( 0, 2 ) =   1.0;
  A( 1, 0 ) = 0.0; A( 1, 1 ) = 2.0; A( 1, 2 ) =   5.0;
  A( 2, 0 ) = 2.0; A( 2, 1 ) = 5.0; A( 2, 2 ) = - 1.0;
  cout << "  * A = " << A << endl;

  b[ 0 ] = 6.0; b[ 1 ] = - 4.0; b[ 2 ] = 27.0;
  cout << "  * b^T = " << b << endl;
  cout << "  * gives the solution vector " << endl;

  x = A.solve_basic( b );
  cout << "  * x^T = " << x << endl;

  cout << "-----------------------------------------------" << endl;

  cout << "  * Solve the complex linear system CX=B for X " << endl;
  cout << "  * where C, X and B are complex matrices. " << endl;

  Matrix< std::complex<double> > C( 2, 2, 0.0 );
  Matrix< std::complex<double> > B( 2, 3, 0.0 );
  Matrix< std::complex<double> > X( 2, 3, 0.0 );

  C.random();
  cout << "  * C = " << C << endl;

  B.random();
  cout << "  * B = " << B << endl;

  cout << "  * gives the solution matrix " << endl;
  X = C.solve_basic( B );
  cout << "  * X = " << X << endl;

  cout << "  * This solution may easily be checked i.e. " << endl;
  cout << "  * CX - B = " << C * X - B << endl;

  cout << "-----------------------------------------------" << endl;

  cout << "  * In theory we may solve systems of any size" << endl;
  cout << "  * but in practice this can be time consuming. " << endl;

  Timer timer;
  double time_in_ms;
  std::size_t N( 16 );
  Matrix<double> D( N, N, 0.0 );

  cout << "  * For example: " << endl;

  for ( std::size_t i = 0; i < 7; ++i )
  {
    D.resize( N, N );
    D.random();
    x.resize( N );
    Vector<double> y( N, 1.0 );

    timer.start();
    x = D.solve_basic( y );
    time_in_ms = timer.get_time();
    timer.stop();

    cout << "  * solving a " << N << "x" << N << " system takes "
         << time_in_ms << " ms. " << endl;
    N *= 2;
  }

  cout << "-----------------------------------------------" << endl;
  cout << "  * Calculate the LU decomposition of a 3x3  " << endl;
  cout << "  * matrix with partial pivoting ( PA = LU )." << endl;
  A( 0, 0 ) = 2.0; A( 0, 1 ) =   4.0; A( 0, 2 ) =  1.0;
  A( 1, 0 ) = 4.0; A( 1, 1 ) = -10.0; A( 1, 2 ) =  2.0;
  A( 2, 0 ) = 1.0; A( 2, 1 ) =   2.0; A( 2, 2 ) =  4.0;
  N = 3;
  cout << "  * A = " << A << endl;

  Matrix<double> L( N, N, 0.0 );
  Matrix<double> U( N, N, 0.0 );
  Matrix<double> P( N, N, 0.0 );

  A.LU_decomp( L, U, P );

  cout << "  * L = " << L << endl;
  cout << "  * U = " << U << endl;
  cout << "  * P = " << P << endl;
  cout << "  * Check that the decomposition has worked." << endl;
  cout << "  * PA - LU = " << P * A - L * U << endl;
  cout << "  * Now caluculate the determinant of the matrix." << endl;
  cout << "  * A.det() = " << A.det() << endl;
  Matrix<double> Inv( N, N, 0.0 );
  Inv = A.inverse();
  cout << "  * Also calculate the inverse and check it." << endl;
  cout << "  * A.inverse() = " << Inv << endl;
  cout << "  * AA^{-1} = " << A * Inv << endl;

  cout << "------------------ FINISHED ------------------" << endl;
}
