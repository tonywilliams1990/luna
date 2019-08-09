/// \file Matrix_test.cpp
/// \ingroup Examples
/// TODO description / rename file

#include "Luna/Core"
#include "Luna/Sparse"

using namespace std;
using namespace Luna;

int main()
{
  cout << "--------------------- Matrix test --------------------" << endl;

  int tests( 100 );
  std::size_t n( 2500 );
  Matrix<double> dense( n, n, 0.0 );

  cout << " * size = " << n << " x " << n << endl;
  cout << " * number of tests = " << tests << endl;

  Vector<double> times( tests );

  Vector<double> b( n, 1.0 );
  Vector<double> x;

  for ( int i = 0; i < tests; i++ )
  {
    //cout << " * test " << i << endl;
    dense.random();

    Timer timer;
    timer.start();

    x = dense.solve( b );
    times[ i ] = timer.get_time(); // Get time in ms
    //timer.print();
    timer.stop();
  }


  cout << " * average time = " << times.sum() / tests << " ms." << endl;



  cout << "--- FINISHED ---" << endl;
}
