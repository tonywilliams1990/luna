/// \file Laplace_equation.cpp
/// \ingroup Examples
/// TODO description

#include "Luna/Core"
#include "Luna/Sparse"

using namespace std;
using namespace Luna;

int main()
{
  cout << "--------------------- Laplace's equation --------------------" << endl;

  int N( 256 );                              // Number of intervals
  double dx( 1.0 / N );                      // Step size
  std::size_t size( ( N + 1 ) * ( N + 1 ) ); // Size of the linear system
  Vector<double> rhs( size, 0.0 );           // Right hand side vector
  Vector<double> u( size, 0.0 );             // Solution vector
  Vector<double> u_exact( size, 0.0 );       // Exact solution
  Vector<double> u_guess( size, 0.0 );       // Guess solution

  typedef Triplet<double> Td;
  Vector<Td> triplets;

  double x( 0.0 );
  double y( 0.0 );
  std::size_t row( 0 );

  // x = 0 boundary ( u = y / ( 1 + y^2 ) )
  std::size_t i( 0 );
  std::size_t j( 0 );
  for ( j = 0; j < N + 1; j++ )
  {
    y = j * dx;
    triplets( Td( row, i * ( N + 1 ) + j, 1.0 ) );
    rhs[ row ] = y / ( 1 + y * y );
    ++row;
  }


  for ( i = 1; i < N; i++ )
  {
    x = i * dx;
    // y = 0 boundary ( u = 0 )
    j = 0;
    triplets( Td( row, i * ( N + 1 ) + j, 1.0 ) );
    rhs[ row ] = 0.0;
    ++row;
    // Interior nodes
    for ( j = 1; j < N; j++ )
    {
      triplets( Td( row, ( i - 1 ) * ( N + 1 ) + j, 1.0 ) );
      triplets( Td( row, ( i + 1 ) * ( N + 1 ) + j, 1.0 ) );
      triplets( Td( row, i * ( N + 1 ) + j - 1, 1.0 ) );
      triplets( Td( row, i * ( N + 1 ) + j + 1, 1.0 ) );
      triplets( Td( row, i * ( N + 1 ) + j, - 4.0 ) );
      rhs[ row ] = 0.0;
      ++row;
    }
    // y = 1 boundary ( u = 1 / ( ( 1 + x )^2 + 1 ) )
    j = N;
    triplets( Td( row, i * ( N + 1 ) + j, 1.0 ) );
    rhs[ row ] = 1. / ( ( 1. + x ) * ( 1. + x ) + 1. );
    ++row;
  }

  // x = 1 boundary ( u = y / ( 4 + y^2 ) )
  i = N;
  for ( j = 0; j < N + 1; j++ )
  {
    y = j * dx;
    triplets( Td( row, i * ( N + 1 ) + j, 1.0 ) );
    rhs[ row ] = y / ( 4 + y * y );
    ++row;
  }

  // Exact solution u = y / ( ( 1 + x )^2 + y^2 )
  row = 0;
  for ( i = 0; i < N + 1; i++ )
  {
    x = i * dx;
    for ( j = 0; j < N + 1; j++ )
    {
      y = j * dx;
      u_exact[ row ] = y / ( ( 1. + x ) * ( 1. + x ) + y * y );
      u_guess[ row ] = y / ( 1. + x * x + y * y );
      ++row;
    }
  }

  // Create the sparse matrix from the triplets
  SparseMatrix<double> sparse( size, size, triplets );

  // Solve with Biconjugate gradient method
  //u.random();
  u = u_guess;
  double tol( 1e-8 );
  int code;
  int max_iter( 10000 );

  double time_start = omp_get_wtime();
  code = sparse.solve_BiCG( rhs, u, max_iter, tol );
  double time_total = omp_get_wtime() - time_start;

  cout << " * time = " << time_total << endl;
  cout << " * iter = " << max_iter << endl;
  cout << " * error estimate = " << scientific << tol << endl;
  u = u - u_exact;
  cout << " * BiCG solution error = " << u.norm_2() << endl;
  cout << " * code = " << code << endl;
  cout << endl;

  // Solve with QMR
  max_iter = 10000;
  //u.random();
  u = u_guess;
  tol = 1e-8;

  time_start = omp_get_wtime();
  code = sparse.solve_QMR( rhs, u, max_iter, tol );
  time_total = omp_get_wtime() - time_start;

  cout << " * time = " << time_total << endl;
  cout << " * iter = " << max_iter << endl;
  u = u - u_exact;
  cout << " * QMR solution error = " << u.norm_2() << endl;
  cout << " * code = " << code << endl;
  cout << endl;

  cout << "--- FINISHED ---" << endl;
}
