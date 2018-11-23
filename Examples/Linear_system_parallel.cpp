/// \file  Linear_system.cpp
/// \ingroup Examples
/// Solve some linear systems of equations

#include "Luna/Core"

//#include "omp.h"

using namespace std;
using namespace Luna;

int main(int argc, char ** argv)
{
  /*int mynode, totalnodes;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
  MPI_Comm_rank(MPI_COMM_WORLD, &mynode);

  if ( mynode == 0 ) // Run this part on a single core
  {

    cout << "----- Linear system parallel -----" << endl;


  }
  // This part is run on every core

  MPI_Barrier( MPI_COMM_WORLD ); // Try to make processes wait but this is not guaranteed

  cout << "Hello world from processor " << mynode << " of " << totalnodes << endl;

  // Delay all processes except 0
  if ( mynode != 0 )
  {
    std::this_thread::sleep_until( std::chrono::system_clock::now()
                                 + std::chrono::seconds(1));
  }

  if ( mynode == 0 )
  {
    cout << "--- FINISHED ---" << endl;
  }
  MPI_Finalize();*/

  /*int num_proc, max_thread;
  num_proc = omp_get_num_procs();
  max_thread = omp_get_max_threads();

  cout << "  * Number of processors available = " << num_proc << endl;
  cout << "  * Number of threads available = " << max_thread << endl;

  // Start of parallel region
  #pragma omp parallel num_threads( max_thread )
  {
    int ID = omp_get_thread_num();
    printf("hello world %d \n", ID);
  }
  // End of parallel region

  std::size_t len( 10 );
  Vector<double> a( len, 1.0 );
  Vector<double> b( len, 2.0 );

  double sum( 0.0 );

  #pragma omp parallel for reduction(+:sum)
  for ( std::size_t i = 0; i < len; ++i )
  {
    sum = sum + a[ i ] + b[ i ];
  }

  cout << "  * parallel sum = " << sum << endl;
  cout << "  * a.sum() + b.sum() = " << a.sum() + b.sum() << endl;

  // Calculate pi by integration

  std::size_t N( 1e8 );
  std::size_t threads( 4 );
  double step( 0.0 );
  double pi( 0.0 );
  sum = 0.0;

  std::size_t i = 0;
  double X;

  step = 1.0 / N;
  omp_set_num_threads( threads );
  double timer_start = omp_get_wtime();

  #pragma omp parallel for private(X) reduction(+:sum)
  for ( i = 0; i < N; ++i )
  {
    X = ( i + 0.5 ) * step;
    sum += 4.0 / ( 1.0 + X * X );
  }

  pi = sum * step;
  double timer_took = omp_get_wtime() - timer_start;
  cout << "  * pi = " << pi << ", " << threads << " threads took "
       << timer_took << " seconds." << endl;
  */

  // Solve a small matrix system in parallel

  Matrix<double> A( 3, 3, 0.0 );
  Vector<double> rhs( 3, 0.0 );
  Vector<double> x;

  cout << "  * Solve the linear system Ax=b for x where" << endl;

  A( 0, 0 ) = 1.0; A( 0, 1 ) = 1.0; A( 0, 2 ) =   1.0;
  A( 1, 0 ) = 0.0; A( 1, 1 ) = 2.0; A( 1, 2 ) =   5.0;
  A( 2, 0 ) = 2.0; A( 2, 1 ) = 5.0; A( 2, 2 ) = - 1.0;
  cout << "  * A = " << A << endl;

  rhs[ 0 ] = 6.0; rhs[ 1 ] = - 4.0; rhs[ 2 ] = 27.0;
  cout << "  * rhs^T = " << rhs << endl;
  cout << "  * gives the solution vector " << endl;

  x = A.solve_parallel( rhs, 4 );
  cout << "  * x^T = " << x << endl; // Should be 5 3 -2

  //TODO solve a large random linear system in parallel

  //Timer timer;
  //double time_in_ms;

  std::size_t N( 4096 );
  Matrix<double> D( N, N, 0.0 );
  D.random();
  Vector<double> y( N, 1.0 );
  //timer.start();
  double time_start = omp_get_wtime();
  x = D.solve_parallel( y, 4 );
  double time_total = omp_get_wtime() - time_start;

  //time_in_ms = timer.get_time();
  //timer.stop();

  cout << "  * Solving a " << N << "x" << N << " system in parallel takes "
       << time_total << " seconds. " << endl;


  cout << "--- FINISHED ---" << endl;
}
