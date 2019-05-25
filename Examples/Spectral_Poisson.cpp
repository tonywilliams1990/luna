/// \file  Spectral_Poisson.cpp
/// \ingroup Examples
/// The equation \f[ \nabla^2 u(x,y) = 2(x^2 + y^2) \f] is solved on the domain
/// \f$ (x,y) \in [-1,1] \times [-1,1] \f$ subject to the boundary conditions
/// \f$ u(-1,y) = e^{-1} \sin(y) + y^2 \f$, \f$ u(1,y) = e \sin(y) + y^2 \f$,
/// \f$ u(x,-1) = e^x \sin(-1) + x^2 \f$ and \f$ u(x,1) = e^x \sin(1) + x^2 \f$.
/// The system is solved using a Chebyshev spectral method such that
/// \f[ u(x,y) \approx u_{IJ}(x,y) = \sum_{N=0}^{IJ-1} a_N \Phi_N(x,y) .\f]
/// Here the basis functions \f$ \Phi_N(x,y) = T_i(x) T_j(y) \f$ are a product
/// of one-dimensional Chebyshev basis functions with \f$i = 0,1,\ldots,I-1 \f$
/// and \f$ j=0,1,\ldots, J-1\f$ where \f$ I \f$ and \f$ J \f$ are the number of
/// collocation points in the \f$ x \f$ and \f$ y \f$ direction respectively.
/// The exact solution is given by \f[ u(x,y) = e^x \sin(y) + x^2 y^2. \f]


#include "Luna/Core"
#include "Luna/ODE"
#include "Luna/Spectral"

namespace Luna
{
  namespace Example
  {

  } // End of namespace Example
} // End of namespace Luna

using namespace std;
using namespace Luna;
using namespace Example;

int main()
{
  cout << "------------------- Spectral_Poisson ------------------" << endl;

  cout << " * The equation u_xx + u_yy = 2(x^2 + y^2) is solved " << endl
       << " * subject to the boundary conditions" << endl
       << " * u(-1,y) = e^(-1) * sin(y) + y^2," << endl
       << " * u(1,y)  = e * sin(y) + y^2," << endl
       << " * u(x,-1) = e^x * sin(-1) + x^2," << endl
       << " * u(x,1)  = e^x * sin(1) + x^2." << endl;

  int I( 6 );           // Number of x collocation points
  int J( 6 );           // Number of y collocation points
  int size( I * J );
  int n( 100 );         // Number of output points (x and y)

  Vector<double> x, y, x_grid, y_grid;
  x.lobatto_grid( I );
  y.lobatto_grid( J );
  x_grid.linspace( -1.0, 1.0, n );
  y_grid.linspace( -1.0, 1.0, n );

  Mesh2D<double> solution( x_grid, y_grid, 2 );

  // Fill in the exact solution
  for ( std::size_t i = 0; i < n; i++ )
  {
    double x( x_grid[ i ] );
    for ( std::size_t j = 0; j < n; j++ )
    {
      double y( y_grid[ j ] );
      solution( i, j, 0 ) = exp( x ) * sin( y ) + x * x * y * y;
    }
  }

  Matrix<double> L( size, size, 0.0 );
  Vector<double> F( size ), a( size );
  Chebyshev<double> cheby;

  double xi, yj;
  int f, g, i, j;

  for ( int M = 0; M < size; M++ )
  {
    i = M / J;
    j = M % J;

    xi = x[ I - ( i + 1 ) ];
    yj = y[ J - ( j + 1 ) ];

    // x = -1 boundary u = e^-1 * sin(y) + y^2 (left)
    if ( i == 0 )
    {
      for ( int N = 0; N < size; N++ )
      {
        f = N / J;
        g = N % J;
        L( M, N ) = cheby( -1.0, f ) * cheby( yj, g );
      }
      F[ M ] = exp( -1.0 ) * sin( yj ) + yj * yj;
    }

    // y = -1 boundary u = e^x * sin(-1) + x^2 (bottom)
    if ( j == 0 && i != 0 && i != I - 1 )
    {
      for ( int N = 0; N < size; N++ )
      {
        f = N / J;
        g = N % J;
        L( M, N ) = cheby( xi, f ) * cheby( -1.0, g );
      }
      F[ M ] = exp( xi ) * sin( -1.0 ) + xi * xi;
    }

    // Internal nodes
    if ( i != 0 && i != I - 1 && j != 0 && j != J - 1 )
    {
      for ( int N = 0; N < size; N++ )
      {
        f = N / J;
        g = N % J;
        L( M, N )  = cheby( xi, f, 2 ) * cheby( yj, g ); // u_xx
        L( M, N ) += cheby( xi, f ) * cheby( yj, g, 2 ); // u_yy
      }
      F[ M ] = 2 * ( xi * xi + yj * yj );
    }

    // y = 1 boundary u = e^x * sin(1) + x^2 (top)
    if ( j == J - 1 && i != 0 && i != I - 1 )
    {
      for ( int N = 0; N < size; N++ )
      {
        f = N / J;
        g = N % J;
        L( M, N ) = cheby( xi, f ) * cheby( 1.0, g );
      }
      F[ M ] = exp( xi ) * sin( 1.0 ) + xi * xi;
    }

    // x = 1 boundary u = e * sin(y) +  y^2 (right)
    if ( i == I - 1 )
    {
      for ( int N = 0; N < size; N++ )
      {
        f = N / J;
        g = N % J;
        L( M, N ) = cheby( 1.0, f ) * cheby( yj, g );
      }
      F[ M ] = exp( 1.0 ) * sin( yj ) + yj * yj;
    }

  }

  // Solve the system for the spectral coefficients
  a = L.solve( F );

  // Output the spectral solution to the 2D mesh
  for ( std::size_t i = 0; i < n; i++ )
  {
    double x( x_grid[ i ] );
    for ( std::size_t j = 0; j < n; j++ )
    {
      double y( y_grid[ j ] );
      for ( std::size_t N = 0; N < size; N++ )
      {
        f = N / J;
        g = N % J;
        solution( i, j, 1 ) += a[ N ] * cheby( x, f ) * cheby( y, g );
      }
    }
  }

  solution.output( "./DATA/Spectral_Poisson.dat" );
  cout << " * For a comparison of the spectral/exact solutions run: " << endl;
  cout << "python Plotting/Spectral_Poisson_plot.py" << endl; 

  cout << "--- FINISHED ---" << endl;
}
