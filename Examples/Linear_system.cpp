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

  for ( std::size_t i = 0; i < 8; ++i )
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

  /*N /= 2;
  D.resize( N, N );
  D.random();
  Vector<double> y( N, 1.0 );
  timer.start();
  x = D.solve_parallel( y, 1 );
  time_in_ms = timer.get_time();
  timer.stop();

  cout << "  * Solving a " << N << "x" << N << " system in parallel takes "
       << time_in_ms << " ms. " << endl;*/


  cout << "------------------ FINISHED ------------------" << endl;
}
