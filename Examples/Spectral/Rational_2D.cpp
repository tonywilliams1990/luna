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

    double a( const double& x, const double& y )
    {
      return x * x + y * y;
    }

    double b( const double& x, const double& y )
    {
      return x * y;
    }

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

  int I( 15 );                       // Number of x collocation points
  int J( I );                       // Number of y collocation points ( J = I )
  double L( 7.0 );                  // Map parameter
  double eps( Example::EPS );
  double tol( 1e-10 );              // Tolerance correction coefficients norm
  int size( I * J );
  int n( 250 );                     // Number of output points (x and y)
  double out_max( 2.0 );

  cout << " * I = J = " << I << endl;
  cout << " * L = " << L << endl;
  cout << " * epsilon = " << eps << endl;

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
  int f, g;

  double norm;
  int max_iter( 100 );
  int iter( 0 );

  Timer timer;

  do{

    timer.start();
    for ( int M = 0; M < size; M++ )
    {
      i = M / J;
      j = M % J;

      xi = x[ I - ( i + 1 ) ];
      yj = y[ J - ( j + 1 ) ];

      // x = 0 boundary u_c = 1 - u_g (left)
      if ( i == 0 )
      {
        for ( int N = 0; N < size; N++ )
        {
          f = N / J;
          g = N % J;
          mat( M, N ) = rationalsemi( 0.0, f ) * rationalsemi( yj, g );
        }
        F[ M ] = 1.0 - u_g( 0.0, yj );
      }

      // y = 0 boundary u_c = 1 - u_g (bottom)
      if ( j == 0 && i != 0 )
      {
        for ( int N = 0; N < size; N++ )
        {
          f = N / J;
          g = N % J;
          mat( M, N ) = rationalsemi( xi, f ) * rationalsemi( 0.0, g );
        }
        F[ M ] = 1.0 - u_g( xi, 0.0 );
      }

      // Internal nodes
      if ( i != 0 && j != 0 )
      {
        for ( int N = 0; N < size; N++ )
        {
          f = N / J;
          g = N % J;
          // u_c_xx
          mat( M, N )  = rationalsemi( xi, f, 2 ) * rationalsemi( yj, g );
          // u_c_yy
          mat( M, N ) += rationalsemi( xi, f ) * rationalsemi( yj, g, 2 );
          // u_g_y * u_c_x
          mat( M, N ) += u_g( xi, yj, 0, 1 ) * rationalsemi( xi, f, 1 )
                                             * rationalsemi( yj, g );
          // u_g_x * u_c_y
          mat( M, N ) += u_g( xi, yj, 1, 0 ) * rationalsemi( xi, f )
                                             * rationalsemi( yj, g, 1 );
          // - eps^2 * a * u_c
          mat( M, N ) += - eps * eps * Example::a( xi, yj )
                         * rationalsemi( xi, f ) * rationalsemi( yj, g );
          // - 2 * eps^2 * b * u_g * u_c
          mat( M, N ) += - 2 * eps * eps * Example::b( xi, yj ) * u_g( xi, yj )
                         * rationalsemi( xi, f ) * rationalsemi( yj, g );
        }
        // - u_g_xx - u_g_yy - u_g_x * u_g_y + eps^2 * u_g * [ a +  b * u_g ]
        F[ M ]  = - u_g( xi, yj, 2, 0 ) - u_g( xi, yj, 0, 2 )
            - u_g( xi, yj, 1, 0 ) * u_g( xi, yj, 0, 1 )
            + eps * eps * u_g( xi, yj ) * ( Example::a( xi, yj )
            + Example::b( xi, yj ) * u_g( xi, yj ) );

      }

      // BCs at infinity are natural boundary conditions

    }

    // Solve the system for the spectral coefficients
    a_c = mat.solve( F );
    u_g.update_coefficients( a_c );
    norm = a_c.norm_2();
    cout << " * iter = " << iter << ", norm = " << std::scientific << norm
         << endl;
    ++iter;
    timer.print();
    timer.stop();

  }while( norm > tol && iter < max_iter );

  // Output the spectral solution to the 2D mesh ( store in 2nd variable )
  u_g( solution, 1 );

  // Output the mesh to a file
  solution.output( "./DATA/Rational_2D.dat" );
  cout << " * For a comparison of the spectral/exact solutions run: " << endl;
  cout << "python Plotting/Rational_2D_plot.py" << endl;

  cout << "--- FINISHED ---" << endl;
}
