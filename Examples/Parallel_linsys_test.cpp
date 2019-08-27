/// \file  Parallel_linsys_test.cpp
/// \ingroup Examples
/// TODO description / rename file

#include "Luna/Core"
#include "Luna/Sparse"
//#include "Luna/Eigen/Dense"

using namespace std;
using namespace Luna;
using namespace Eigen;

int main()
{
  cout << "----------------- Eigen linear system test ----------------" << endl;

  int n = 512;
  Timer timer;
  /*MatrixXd A = MatrixXd::Random(n,n);
  //cout << " A = " << endl << A << endl;
  MatrixXd b = MatrixXd::Random(n,50);
  timer.start();
  MatrixXd x = A.partialPivLu().solve(b);
  timer.print();
  timer.stop();
  double relative_error = (A*x - b).norm() / b.norm(); // norm() is L2 norm
  cout << "The relative error is: " << relative_error << endl;*/

  Luna::Matrix<double> mat( n, n, 0.0 );
  mat.random();
  //cout << " mat = " << mat << endl;
  Luna::Vector<double> b_vec( n, 1.0 );
  Luna::Vector<double> x_vec( n );

  //cout << " b_vec = " << b_vec << endl;
  double b_norm( b_vec.norm_2() );

  timer.start();
  //x_vec = mat.solve_parallel( b_vec );
  x_vec = mat.solve( b_vec );
  timer.print();
  timer.stop();

  b_vec -= mat * x_vec;
  cout << " relative error = " << std::scientific << b_vec.norm_2() / b_norm << endl;

  cout << "--- FINISHED ---" << endl;
}
