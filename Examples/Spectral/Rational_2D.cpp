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
    double a( const double& x, const double& y )
    {
      //return - x * y * exp( - x * y );
      //return 2 * ( pow( x + 1, - 2.0) + pow( y + 1, - 2.0 ) );
      //return 2 * ( pow( y + 1, - 2.0 ) - ( pow( x + 1, - 2.0) / x ) );
      return 2 * pow( y + 1, - 2.0 ) + ( x * x + 4. * x + 5. ) * pow( x + 1, - 2.0);
    }

    double b( const double& x, const double& y )
    {
      return ( x * y - 2.0 ) * ( x * x + y * y ) * exp( - x * y );
    }

    double initial_guess( const double& x, const double& y )
    {
      //return x * y * exp( - x * y );
      return exp( - x ) / ( ( x + 1. ) * ( y + 1. ) );
    }

    double exact( const double& x, const double& y )
    {
      //return x * y * exp( - x * y );
      //return exp( - x * y );
      //return 1. / ( ( x + 1. ) * ( y + 1. ) );
      //return x / ( ( x + 1. ) * ( y + 1. ) );
      return exp( - x ) / ( ( x + 1. ) * ( y + 1. ) );
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
  int J( 10 );                       // Number of y collocation points
  double L( 1.0 );                  // Map parameter
  double tol( 1e-10 );              // Tolerance correction coefficients norm
  int size( I * J );
  int n( 200 );                     // Number of output points (x and y)
  double out_max( 10.0 );

  Vector<double> x, y, x_grid, y_grid;
  x.rational_semi_grid( I - 1, L );
  y.rational_semi_grid( J - 1, L );
  x.push_back( 0.0 );
  y.push_back( 0.0 );
  //x.reverse();
  //y.reverse();
  cout << " x = " << x << endl;
  cout << " y = " << y << endl;

  //x_grid.powspace( 0.0, x[ 0 ], n, 2 );
  //y_grid.powspace( 0.0, y[ 0 ], n, 2 );
  x_grid.linspace( 0.0, out_max, n );
  y_grid.linspace( 0.0, out_max, n );
  Mesh2D<double> solution( x_grid, y_grid, 2 );

  // Fill in the exact solution (store in first variable)
  solution.apply( Example::exact, 0 );

  Matrix<double> mat( size, size, 0.0 );
  Vector<double> F( size ), a_c( size );
  RationalSemi<double> rationalsemi( L );
  Vector<double> c( size, 0.0 );                 // Vector of coefficients

  // Approximate the guess as a semi-infinite rational Chebyshev polynomial
  //c = rationalsemi.approximate2D( Example::initial_guess, I, J );

  Spectral2D<double> u_g( c, I, J, "RationalSemi" ); // Initially zero everywhere



  double xi, yj;
  int f, g, i, j;

  double norm;
  int max_iter( 1 );
  int iter( 0 );

  do{

    for ( int M = 0; M < size; M++ )
    {
      i = M / J;
      j = M % J;

      xi = x[ I - ( i + 1 ) ];
      yj = y[ J - ( j + 1 ) ];
      //cout << "(i,j) = " << i << ", " << j << "\t (x,y) = " << xi << ", " << yj << endl;

      // x = 0 boundary u = 0 (left)
      if ( i == 0 )
      {
        for ( int N = 0; N < size; N++ )
        {
          f = N / J;
          g = N % J;
          mat( M, N ) = rationalsemi( 0.0, f ) * rationalsemi( yj, g );
          //mat( M, N ) = 1.23;
        }
        //F[ M ] = - u_g( 0.0, yj );
        //F[ M ] = 1. / ( yj + 1. );
        F[ M ] = 1. / ( yj + 1. );
      }

      // y = 0 boundary u = 0 (bottom)
      if ( j == 0 && i != 0 )
      {
        for ( int N = 0; N < size; N++ )
        {
          f = N / J;
          g = N % J;
          mat( M, N ) = rationalsemi( xi, f ) * rationalsemi( 0.0, g );
        }
        //F[ M ] = - u_g( xi, 0.0 );
        //F[ M ] = 1. / ( xi + 1. );
        //F[ M ] = xi / ( xi + 1. );
        F[ M ] = exp( - xi ) / ( xi + 1. );
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
          //mat( M, N ) += - 2. * u_g( xi, yj ) * rationalsemi( xi, f ) * rationalsemi( yj, g ); // -2 * u_g * u_c
          //mat( M, N ) += - Example::a( xi, yj ) * rationalsemi( xi, f ) * rationalsemi( yj, g ); // - a * u_c
          mat( M, N ) += - Example::a( xi, yj ) * rationalsemi( xi, f ) * rationalsemi( yj, g );
        }
        // - u_g_xx - u_g_yy + u_g^2 + a * u_g + b
        F[ M ]  = //- u_g( xi, yj, 2, 0 ) - u_g( xi, yj, 0, 2 )
                  //+ u_g( xi, yj ) * u_g( xi, yj )
                  //+ Example::a( xi, yj ) * u_g( xi, yj ) +
                  //Example::b( xi, yj );
                  0.0;

      }

      // BCs at infinity are natural boundary conditions

    }

    // Solve the system for the spectral coefficients
    a_c = mat.solve( F );
    //cout << " a_c = " << a_c << endl;
    u_g.update_coefficients( a_c );
    norm = a_c.norm_2();
    cout << " norm = " << norm << endl;
    ++iter;

  }while( norm > tol && iter < max_iter );



  // Create the 2D spectral solution
  //Spectral2D<double> u( a, I, J, "Chebyshev" );

  // Output the spectral solution to the 2D mesh ( store in 2nd variable )
  u_g( solution, 1 );

  // Output the mesh to a file
  solution.output( "./DATA/Rational_2D.dat" );
  //cout << " * For a comparison of the spectral/exact solutions run: " << endl;
  //cout << "python Plotting/Spectral_Poisson_plot.py" << endl;

  cout << "--- FINISHED ---" << endl;
}
