/// \file  ODE_BVP_Nonlinear.cpp
/// \ingroup Examples
/// \ingroup Spectral
/// The equation \f[ u(x) u''(x) - \left( u'(x) \right)^2 = -1 \f] is solved
/// subject to the boundary conditions \f$ u(1) = \frac{1}{2} \left(e - e^{-1}
/// \right) = - u(-1). \f$ The system is solved using a Chebyshev spectral
/// method. The solution \f$ u(x) \f$ is anti-symmetric so only odd Chebyshev
/// polynomials are needed to approximate the solution such that \f[ u(x)
/// \approx u_n(x) = \sum_{i=0}^n a_i T_{2i+1}(x) .\f] Here \f$ n \f$ is the
/// number of spectral coefficients \f$ a_i \f$ and \f$ T_{2i+1}(x) \f$ are the
/// corresponding Chebyshev polynomials. The solution to this nonlinear equation
/// is found by splitting it into a known part \f$ u^g(x) \f$ and a correction
/// \f$ u^c(x) \f$, such that \f$ u(x) = u^g(x) + u^c(x)\f$ , and using Newton
/// iteration to update the known part until convergence. The correction is
/// found from the resulting linearised equations. The exact solution is given by
/// \f$ u(x) = \sinh(x). \f$ The spectral solution is found on \f$ x \in [0,1]
/// \f$ and then because the solution is anti-symmetric this gives the solution
/// in the entire domain \f$ x \in [-1,1]. \f$

#include "Luna/Core"
#include "Luna/ODE"
#include "Luna/Spectral"

namespace Luna
{
  namespace Example
  {
    double initial_guess( const double& x )
    {
      double e( std::exp( 1.0 ) );
      return 0.5 * ( e - ( 1 / e ) ) * x;
    }
  } // End of namespace Example
} // End of namespace Luna

using namespace std;
using namespace Luna;
using namespace Example;

int main()
{
  cout << "------------------ Spectral_ODE_Nonlinear -----------------" << endl;

  cout << " * The equation u * u'' + (u')^2 = -1 is solved subject to" << endl
       << " * the boundary conditions u(1) = (e - e^(-1)) / 2 = -u(-1)." << endl
       << " * Only anti-symmetric odd Chebyshev polynomials are used in" << endl
       << " * the spectral solution. " << endl;

  int n( 6 );                       // Number of spectral coefficients
  double tol( 1e-10 );              // Tolerance correction coefficients norm
  Vector<double> x;                 // Vector for collocation grid
  x.half_lobatto_grid( n );

  int N( 100 );                         // Number of output points
  Vector<double> grid;                  // Output grid
  grid.linspace( -1.0, 1.0, N );
  Mesh1D<double> solution( grid, 3 );   // Output mesh

  // Exact solution
  for ( std::size_t i = 0; i < N; i++ )
  {
    solution( i, 0 ) = std::sinh( grid[ i ] );
  }

  Chebyshev<double> cheby;

  Vector<double> c( n, 0.0 );                 // Vector of coefficients

  // Approximate the initial guess as a Chebyshev polynomial
  c = cheby.approximate_odd( Example::initial_guess, n );
  // Make this into a spectral solution
  Spectral<double> u_g( c, "OddChebyshev" );
  double norm;

  do
  {
    // Setup the system of equations to find the correction u_c
    Matrix<double> L( n, n, 0.0 );
    Vector<double> F( n ), a_c( n );

    // x = 1 boundary u(1) = 0.5 * ( e - e^(-1) )
    for ( std::size_t j = 0; j < n; j++ )
    {
      L( 0, j ) = cheby( 1.0, 2 * j + 1 );
    }
    F[ 0 ] = 0.5 * ( std::exp( 1.0 ) - std::exp( - 1.0 ) ) - u_g( 1.0 );

    // Internal nodes and x = 0 anti-symmetric boundary
    for ( std::size_t i = 1; i < n; i++ )
    {
      double xi = x[ i ];
      for ( std::size_t j = 0; j < n; j ++ )
      {
        // u_g * u_c''
        L( i, j ) = u_g( xi ) * cheby( xi, 2 * j + 1, 2 );
        // -2 * u_g' * u_c'
        L( i, j ) += - 2.0 * u_g( xi, 1 ) * cheby( xi, 2 * j + 1, 1 );
        // u_g'' * u_c
        L( i, j ) += u_g( xi, 2 ) * cheby( xi, 2 * j + 1 );
      }
      // u_g' * u_g' - u_g * u_g'' - 1
      F[ i ] = u_g( xi, 1 ) * u_g( xi, 1 ) - u_g( xi ) * u_g( xi, 2 ) - 1.0;
    }

    // Solve the system for the correction spectral coefficients
    a_c = L.solve( F );
    u_g.update_coefficients( a_c );
    norm = a_c.norm_2();

  }while( norm > tol );

  cout << " * The spectral coefficients for n = " << n << " are: " << endl;

  for ( std::size_t i = 0; i < n; i++ )
  {
    cout << " * a_" << i << " = " << std::scientific
         << u_g.get_coefficients()[ i ] << endl;
  }

  // Put spectral solution into output mesh
  for ( std::size_t i = 0; i < N; i++ )
  {
    solution( i, 1 ) = u_g( grid[ i ] );
    solution( i, 2 ) = solution( i, 1 ) - solution( i, 0 ); // difference
  }

  solution.output( "./DATA/Spectral_ODE_Nonlinear.dat" );

  cout << " * For a comparison of the spectral/exact solutions run: " << endl;
  cout << "python Plotting/Spectral_ODE_Nonlinear_plot.py" << endl;

  cout << "--- FINISHED ---" << endl;
}
