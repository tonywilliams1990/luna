/// \file Falkner.cpp
/// \ingroup Examples
/// \ingroup Spectral
/// The Falkner-Skan equation \f[ f'''(\eta) + f(\eta)f''(\eta) + \beta \left[1
/// - \left(f'(\eta) \right)^2 \right] = 0 \f] is solved on \f$ \eta \in
/// [0,\infty) \f$ subject to the boundary conditions \f$ f(0) = f'(0) =  0\f$
/// and \f$ f'(\eta) \to 1 \f$ as \f$ \eta \to \infty\f$.
/// The system is solved using rational Chebyshev functions on a semi-infinite
/// interval \f$ TL_i(\eta; L) \f$. The rational chebyshev functions all have zero
/// derivative as \f$ \eta \to \infty \f$ so it is necessary to use the change
/// of variables \f$ u(\eta) = f(\eta) - \eta \f$. The Falkner-Skan system is
/// transformed to \f[ u''' + (\eta + u)u'' - \beta u'\left(2 + u' \right) = 0
/// \f] subject to the boundary conditions \f$ u(0) = 0\f$, \f$ u'(0) = -1 \f$
/// and \f$ u' \to 0 \f$ as \f$ \eta \to \infty \f$. This nonlinear equation is
/// solver using Newton iteration on the spectral coefficients.

#include "Luna/Core"
#include "Luna/ODE"
#include "Luna/Spectral"
#include "Luna/Falkner.h"

namespace Luna
{
  namespace Example
  {
    double initial_guess( const double& y )
    {
      return - y * exp( - y );
    }
  } // End of namespace Example
} // End of namespace Luna

using namespace std;
using namespace Luna;
using namespace Example;

int main()
{
  cout << "------------------------- Falkner -------------------------" << endl;

  cout << " * The Falkner-Skan equation, for a boundary layer on flat " << endl
       << " * plate, is solved using a rational Chebyshev spectral " << endl
       << " * method." << endl;

  double beta( 0.0 );               // Hartree parameter
  int n( 40 );                      // Number of spectral coefficients
  double L( 6.0 );                  // Map parameter
  double tol( 1e-10 );              // Tolerance correction coefficients norm
  int max_iter( 100 );              // Maximum number of iterations when solving

  int N( 1000 );                        // Number of output points
  double eta_max( 20.0 );               // Maximum eta value to output
  Vector<double> grid;                  // Output grid
  grid.linspace( 0.0, eta_max, N );
  Mesh1D<double> solution( grid, 3 );   // Output mesh

  Falkner<double> falkner( beta );
  Spectral<double> u_g;

  Timer timer;
  timer.start();

  // Solve for the spectral solution
  u_g = falkner.solve( n, L, tol, Example::initial_guess, max_iter );

  cout << " * The final spectral coefficient for n = " << n << " is: "
       << u_g.get_coefficients()[ n-1 ] << endl;
  cout << " * f''(0) = " << u_g( 0.0, 2 ) << endl;

  // Put spectral solution into output mesh
  for ( std::size_t i = 0; i < N; i++ )
  {
    double y( grid[ i ] );
    solution( i, 0 ) = y + u_g( y );              // f   = eta + u
    solution( i, 1 ) = 1 + u_g( y, 1 );           // f'  = 1 + u'
    solution( i, 2 ) = u_g( y, 2 );               // f'' = u''
  }

  solution.output( "./DATA/Spectral_Falkner.dat" );

  timer.print();
  timer.stop();

  cout << " * To see the spectral solution and its derivatives run: " << endl;
  cout << "python Plotting/Spectral_Falkner_plot.py" << endl;

  cout << "--- FINISHED ---" << endl;
}
