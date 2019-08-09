/// \file Rational_2D.cpp
/// \ingroup Examples
/// \ingroup Spectral
/// \todo TODO description

#include "Luna/Core"
#include "Luna/ODE"
#include "Luna/Spectral"

namespace Luna
{
  namespace Example
  {
    double EPS( 0.1 );

    double initial_guess_x( const double& x )
    {
      return exp( - EPS * x );
    }

    double initial_guess_y( const double& y )
    {
      return exp( - EPS * y );
    }

    double exact( const double& x, const double& y )
    {
      return exp( - EPS * x * y );
    }
  } // End of namespace Example
} // End of namespace Luna

using namespace std;
using namespace Luna;
using namespace Example;

int main()
{
  cout << "----------------------- Rational_2D -----------------------" << endl;
  cout << " * The equation u_xx + u_yy + u_x * u_y = " << endl
       << " * eps^2 (x^2 + y^2) * u + eps^2 * x * y * u^2 " << endl
       << " * is solved subject to the boundary conditions " << endl
       << " * u(x,0) = u(0,y) = 1 and u -> 0 as x,y -> infinity." << endl
       << " * epsilon is a given small parameter. The exact solution" << endl
       << " * is u(x,y) = e^(-eps*x*y)." << endl;

  int I( 60 );                       // Number of x collocation points
  int J( I );                       // Number of y collocation points ( J = I )
  double L( 7.0 );                  // Map parameter
  double eps2( Example::EPS * Example::EPS );
  double tol( 1e-10 );              // Tolerance correction coefficients norm
  int size( I * J );
  int n( 250 );                     // Number of output points (x and y)
  double out_max( 2.0 );

  Timer total_timer;
  total_timer.start();

  cout << " * I = J = " << I << endl;
  cout << " * L = " << L << endl;
  cout << " * epsilon = " << Example::EPS << endl;

  Vector<double> x, y, x_grid, y_grid;
  x.rational_semi_grid( I - 1, L );
  y.rational_semi_grid( J - 1, L );
  x.push_back( 0.0 );
  y.push_back( 0.0 );

  x_grid.linspace( 0.0, out_max, n );
  y_grid.linspace( 0.0, out_max, n );
  Mesh2D<double> solution( x_grid, y_grid, 3 );

  // Fill in the exact solution (store in first variable)
  solution.apply( Example::exact, 0 );

  Matrix<double> mat( size, size, 0.0 );
  Vector<double> F( size ), a_c( size );
  RationalSemi<double> rationalsemi( L );
  Vector<double> c( size, 0.0 );                 // Vector of coefficients
  int i, j;

  // Approximate solution as product of 1D spectral functions
  Vector<double> c_x( I, 0.0 );
  c_x = rationalsemi.approximate( Example::initial_guess_x, I );

  Vector<double> c_y( J, 0.0 );
  c_y = rationalsemi.approximate( Example::initial_guess_y, J );

  for ( int k = 0; k < size; k++ )
  {
    i = k / J;
    j = k % J;
    c[ k ] = c_x[ i ] * c_y[ j ];
  }

  // Set the spectral solution guess
  Spectral2D<double> u_g( c, I, J, "RationalSemi", L );
  u_g( solution, 2 ); // output initial guess to mesh

  double xi, yj;
  int f, g, ix, jy;

  double norm;
  int max_iter( 100 );
  int iter( 0 );

  Vector<double> phix, phiy;

  Timer timer;

  do{

    cout << "---------------------------------------------------------" << endl;
    timer.start();

    // Make mesh for storing u_g, u_g_x, u_g_y, u_g_xx, u_g_yy, u_g_xy
    Mesh2D<double> u_g_mesh = u_g.mesh_derivatives( x, y );
    timer.print( "Solution mesh construction time" );

    for ( int M = 0; M < size; M++ )
    {
      i = M / J;
      j = M % J;

      ix = I - ( i + 1 );
      jy = J - ( j + 1 );

      xi = x[ ix ];
      yj = y[ jy ];

      // Matrix
      for ( int N = 0; N < size; N++ )
      {
        f = N / J;
        g = N % J;

        phix = rationalsemi.eval_2( xi, f );  // phix, phix' and phix''
        phiy = rationalsemi.eval_2( yj, g );  // phiy, phiy' and phiy''

        // x = 0 boundary u_c = 1 - u_g (left)
        if ( i == 0 )
        {
          mat( M, N ) = rationalsemi( 0.0, f ) * phiy[ 0 ];
        }
        // y = 0 boundary u_c = 1 - u_g (bottom)
        if ( j == 0 && i != 0 )
        {
          mat( M, N ) = phix[ 0 ] * rationalsemi( 0.0, g );
        }
        // Internal nodes
        if ( i != 0 && j != 0 )
        {
          // u_c_xx
          mat( M, N )  = phix[ 2 ] * phiy[ 0 ];
          // u_c_yy
          mat( M, N ) += phix[ 0 ] * phiy[ 2 ];
          // u_g_y * u_c_x
          /*mat( M, N ) += u_g( xi, yj, 0, 1 ) * rationalsemi( xi, f, 1 )
                                             * rationalsemi( yj, g );*/
          mat( M, N ) += u_g_mesh( ix, jy, 2 )
                         * phix[ 1 ] * phiy[ 0 ];
          // u_g_x * u_c_y
          /*mat( M, N ) += u_g( xi, yj, 1, 0 ) * rationalsemi( xi, f )
                                             * rationalsemi( yj, g, 1 );*/
          mat( M, N ) += u_g_mesh( ix, jy, 1 )
                         * phix[ 0 ] * phiy[ 1 ];
          // - eps^2 * a * u_c
          mat( M, N ) += - eps2 * ( xi * xi + yj * yj )
                         * phix[ 0 ] * phiy[ 0 ];
          // - 2 * eps^2 * b * u_g * u_c
          /*mat( M, N ) += - 2 * eps * eps * Example::b( xi, yj ) * u_g( xi, yj )
                         * rationalsemi( xi, f ) * rationalsemi( yj, g );*/
          mat( M, N ) += - 2 * eps2 * xi * yj * u_g_mesh( ix, jy, 0 )
                         * phix[ 0 ] * phiy[ 0 ];
        }
      }

      // Residuals
      // x = 0 boundary u_c = 1 - u_g (left)
      if ( i == 0 )
      {
        F[ M ] = 1.0 - u_g_mesh( I - 1, jy, 0 );
      }
      // y = 0 boundary u_c = 1 - u_g (bottom)
      if ( j == 0 && i != 0 )
      {
        F[ M ] = 1.0 - u_g_mesh( ix, J - 1, 0 );
      }
      // Internal nodes
      if ( i != 0 && j != 0 )
      {
        // - u_g_xx - u_g_yy - u_g_x * u_g_y
        // + eps^2 * u_g * ( x^2 + y^2 + x * y * u_g )
        F[ M ] = - u_g_mesh( ix, jy, 3 ) - u_g_mesh( ix, jy, 4 )
                 - u_g_mesh( ix, jy, 1 ) * u_g_mesh( ix, jy, 2 )
                 + eps2 * u_g_mesh( ix, jy, 0 ) * ( xi * xi + yj * yj
                 + xi * yj * u_g_mesh( ix, jy, 0 ) );
      }

      // BCs at infinity are natural boundary conditions

    }
    timer.print( "Matrix construction time" );
    // Solve the system for the spectral coefficients
    //a_c = mat.solve( F );
    a_c = mat.solve_parallel( F );
    u_g.update_coefficients( a_c );
    norm = a_c.norm_2();
    cout << "  * iter = " << iter << ", norm = " << std::scientific << norm
         << std::fixed << endl;
    ++iter;
    timer.print( "Total time" );
    timer.stop();

  }while( norm > tol && iter < max_iter );

  cout << "---------------------------------------------------------" << endl;
  // Output the spectral solution to the 2D mesh ( store in 2nd variable )
  u_g( solution, 1 );

  // Output the mesh to a file
  solution.output( "./DATA/Rational_2D.dat" );

  total_timer.print();
  total_timer.stop();

  cout << " * For a comparison of the spectral/exact solutions run: " << endl;
  cout << "python Plotting/Rational_2D_plot.py" << endl;

  cout << "--- FINISHED ---" << endl;
}
