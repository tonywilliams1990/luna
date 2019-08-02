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
    double a( const double& x )
    {
      return 2. / ( x + 1. );
    }

    double b( const double& y )
    {
      return 2. / ( y + 1. );
    }

    double c( const double& x, const double& y )
    {
      double sx( sin( x ) );
      double sy( sin( y ) );
      return sx * sx * sy * sy * pow( x + 1., - 2. ) * pow( y + 1., - 2. );
    }

    double initial_guess( const double& x, const double& y )
    {
      return exp( - x ) / ( ( x + 1. ) * ( y + 1. ) );
    }

    double initial_guess_x( const double& x )
    {
      //return 1. - exp( - x );
      //return x * pow( x + 1., - 2. );
      return sin( x ) / ( x + 1. );
    }

    double initial_guess_y( const double& y )
    {
      //return 1. - exp( - y );
      //return y * pow( y + 1., - 2. );
      return sin( y ) / ( y + 1. );
    }

    double exact( const double& x, const double& y )
    {
      //return x * y / ( ( x + 1. ) * ( y + 1. ) );
      return sin( x ) * sin( y ) / ( ( x + 1. ) * ( y + 1. ) );
    }
  } // End of namespace Example
} // End of namespace Luna

using namespace std;
using namespace Luna;
using namespace Example;

int main()
{
  cout << "----------------------- Rational_2D -----------------------" << endl;

  /// \todo TODO output equation being solved and boundary conditions

  int I( 10 );                       // Number of x collocation points
  int J( I );                       // Number of y collocation points ( J = I )
  double L( 1.0 );                  // Map parameter
  double tol( 1e-10 );              // Tolerance correction coefficients norm
  int size( I * J );
  int n( 300 );                     // Number of output points (x and y)
  double out_max( 10.0 );

  Vector<double> x, y, x_grid, y_grid;
  x.rational_semi_grid( I - 1, L );
  y.rational_semi_grid( J - 1, L );
  x.push_back( 0.0 );
  y.push_back( 0.0 );
  cout << " x = " << x << endl;
  cout << " y = " << y << endl;

  //x_grid.powspace( 0.0, x[ 0 ], n, 2 );
  //y_grid.powspace( 0.0, y[ 0 ], n, 2 );
  x_grid.linspace( 0.0, out_max, n );
  y_grid.linspace( 0.0, out_max, n );
  Mesh2D<double> solution( x_grid, y_grid, 3 );

  // Fill in the exact solution (store in first variable)
  solution.apply( Example::exact, 0 );

  Matrix<double> mat( size, size, 0.0 );
  Vector<double> F( size ), a_c( size );
  RationalSemi<double> rationalsemi( L );
  Vector<double> c( size, 0.0 );                 // Vector of coefficients

  Vector<double> c_x( I, 0.0 );
  c_x = rationalsemi.approximate( Example::initial_guess_x, I );
  cout << " c_x = " << c_x << endl;

  Vector<double> c_y( J, 0.0 );
  c_y = rationalsemi.approximate( Example::initial_guess_y, J );
  cout << " c_y = " << c_y << endl;

  int i, j;
  for ( int k = 0; k < size; k++ )
  {
    i = k / J;
    j = k % J;
    c[ k ] = c_x[ i ] * c_y[ j ];
  }
  cout << " c = " << c << endl;

  // Approximate the guess as a semi-infinite rational Chebyshev polynomial
  //c = rationalsemi.approximate2D( Example::initial_guess, I, J );

  // Set the spectral solution guess
  Spectral2D<double> u_g( c, I, J, "RationalSemi" );
  u_g( solution, 2 ); // output initial guess to mesh

  // Test 1D approximation
  /*Spectral<double> u_g_x( c_x, "RationalSemi", L );
  Spectral<double> u_g_y( c_y, "RationalSemi", L );
  cout << " u_g_x( 0.0 ) = " << std::scientific << u_g_x( 0.0 ) << endl;
  cout << " u_g_y( 0.0 ) = " << std::scientific << u_g_y( 0.0 ) << endl;
  Mesh1D<double> solution1D( x_grid, 3 );   // Output mesh
  for ( std::size_t i = 0; i < n; i++ )
  {
    double x( x_grid[ i ] );
    double y( y_grid[ i ] );
    solution1D( i, 0 ) = u_g_x( x );                     // u_g_x
    solution1D( i, 1 ) = u_g_y( y );                     // u_g_y
    //solution1D( i, 2 ) = sqrt( x ) / ( sqrt( x ) + 1. ); // exact ( initial guess )
    solution1D( i, 2 ) = 1 - exp( - x ); // exact ( initial guess )
  }
  solution1D.output( "./DATA/Initial_guess_xy.dat" );*/

  double xi, yj;
  int f, g;

  double norm;
  int max_iter( 100 );
  int iter( 0 );

  do{

    for ( int M = 0; M < size; M++ )
    {
      i = M / J;
      j = M % J;

      xi = x[ I - ( i + 1 ) ];
      yj = y[ J - ( j + 1 ) ];

      // x = 0 boundary u_c = - u_g (left)
      if ( i == 0 )
      {
        for ( int N = 0; N < size; N++ )
        {
          f = N / J;
          g = N % J;
          mat( M, N ) = rationalsemi( 0.0, f ) * rationalsemi( yj, g );
        }
        F[ M ] = - u_g( 0.0, yj );
        //cout << " u_g( 0.0, " << yj << " ) = " << u_g( 0.0, yj ) << endl;
      }

      // y = 0 boundary u_c = - u_g (bottom)
      if ( j == 0 && i != 0 )
      {
        for ( int N = 0; N < size; N++ )
        {
          f = N / J;
          g = N % J;
          mat( M, N ) = rationalsemi( xi, f ) * rationalsemi( 0.0, g );
        }
        F[ M ] = - u_g( xi, 0.0 );
        //cout << " u_g( " << xi << ", 0.0 ) = " << u_g( xi, 0.0 ) << endl;
      }

      // Internal nodes
      if ( i != 0 && j != 0 )
      {
        for ( int N = 0; N < size; N++ )
        {
          f = N / J;
          g = N % J;
          mat( M, N )  = rationalsemi( xi, f, 2 ) * rationalsemi( yj, g ); // u_c_xx
          mat( M, N ) += rationalsemi( xi, f ) * rationalsemi( yj, g, 2 ); // u_c_yy
          mat( M, N ) += Example::a( xi ) * rationalsemi( xi, f, 1 ) * rationalsemi( yj, g ); // a * u_c_x
          mat( M, N ) += Example::b( yj ) * rationalsemi( xi, f ) * rationalsemi( yj, g, 1 ); // b * u_c_y
          //mat( M, N ) += 2 * u_g( xi, yj ) * rationalsemi( xi, f ) * rationalsemi( yj, g ); // 2 * u_g * u_c
          mat( M, N ) += 2 * ( 1. + u_g( xi, yj ) ) * rationalsemi( xi, f ) * rationalsemi( yj, g ); // 2 * ( 1 + u_g ) * u_c
        }
        // c - u_g_xx - u_g_yy - a * u_g_x - b * u_g_y - u_g^2
        F[ M ]  = Example::c( xi, yj )
                - u_g( xi, yj, 2, 0 ) - u_g( xi, yj, 0, 2 )
                - Example::a( xi ) * u_g( xi, yj, 1, 0 )
                - Example::b( yj ) * u_g( xi, yj, 0, 1 )
                //- u_g( xi, yj ) * u_g( xi, yj ) ;
                - u_g( xi, yj ) * ( 2. + u_g( xi, yj ) ) ;

      }

      // BCs at infinity are natural boundary conditions

    }

    // Solve the system for the spectral coefficients
    a_c = mat.solve( F );
    //cout << " a_c = " << a_c << endl;
    u_g.update_coefficients( a_c );
    norm = a_c.norm_2();
    cout << iter << " norm = " << std::scientific << norm << endl;
    cout << u_g.get_coefficients() << endl;
    ++iter;

  }while( norm > tol && iter < max_iter );

  // Output the spectral solution to the 2D mesh ( store in 2nd variable )
  u_g( solution, 1 );

  // Output the mesh to a file
  solution.output( "./DATA/Rational_2D.dat" );
  //cout << " * For a comparison of the spectral/exact solutions run: " << endl;
  //cout << "python Plotting/Spectral_Poisson_plot.py" << endl;

  cout << "--- FINISHED ---" << endl;
}
