/// \file Laplace_equation.cpp
/// \ingroup Examples
/// Solve Laplace's equation  \f[ \nabla^2 u = 0, \f] on the unit square
/// \f$ (x,y) \in [0,1] \times [0,1] \f$. The equation is subject to the
/// boundary conditions \f[ u(x,0) = 0, \hspace{0.5cm}
/// u(x,1) = \frac{1}{(1+x)^2 + 1}, \hspace{0.5cm} 0 \leq x \leq 1, \f]
/// \f[ u(0,y) = \frac{y}{1+y^2}, \hspace{0.5cm} u(1,y) = \frac{y}{4+y^2},
/// \hspace{0.5cm} 0 \leq y \leq 1. \f] The exact solution is given by \f[
/// u(x,y) = \frac{y}{(1+x)^2 + y^2}, \hspace{0.5cm} 0 \leq x,y \leq 1. \f]
/// The equation is discretised using a second order finite-difference scheme
/// with step size \f$ \Delta x = \Delta y = 1 / N \f$ where \f$ N \f$ is the
/// number of intervals. The internal nodes in the discretisation are governed
/// by the finite-difference equation \f[ u_{i-1,j} + u_{i,j-1} - 4u_{i,j} +
/// u_{i,j+1} + u_{i+1,j} = 0, \f] where \f$ u(x,y) \approx u(x_i,y_j) =
/// u_{i,j} \f$ and \f$ x_i = i \Delta x, \hspace{0.2cm} y_j = j \Delta y \f$
/// with \f$ i,j = 0,1,\ldots,N. \f$

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
      ++row;
    }
  }

  // Create the sparse matrix from the triplets
  SparseMatrix<double> sparse( size, size, triplets );

  // Solve using the Biconjugate gradient method
  u.random();                                    // Random initial guess
  double tol( 1e-8 );
  int code;
  int max_iter( 10000 );

  Timer timer;
  double time_in_ms;
  timer.start();
  code = sparse.solve_BiCG( rhs, u, max_iter, tol );
  time_in_ms = timer.get_time();
  timer.stop();

  cout << endl << " * time = " << time_in_ms / 1000 << " s" << endl;
  cout << " * iter = " << max_iter << endl;
  cout << " * system error estimate = " << scientific << tol << endl;
  u = u - u_exact;
  cout << " * solution error = " << u.norm_2() << endl;
  cout << " * code ( 0 = success ) = " << code << endl << endl;
  cout << "--- FINISHED ---" << endl;
}
