/// \file Boundary_region.cpp
/// \ingroup Examples
/// \ingroup Spectral
/// \todo TODO description

#include "Luna/Core"
#include "Luna/ODE"
#include "Luna/Spectral"
#include "Luna/Falkner.h"

enum{ UB, UBd, PhiB, ThetaB, ThetaBd, PsiB };                 // Base ODE
enum{ Phi, Psi, U, Theta };                                   // Correction PDE

namespace Luna
{
  namespace Example
  {
    double EPS( 0.1 );
    double GAMMA( 20.0 );
    double K( 0.0 );

    double base_guess( const double& y )
    {
      return - y * exp( - y );
    }

    double Phi_w_func( const double& hzeta ){
      // Top-hat injection
      return - K * 0.5 * ( tanh( GAMMA * ( hzeta - 1. ) )
             - tanh( GAMMA * ( hzeta - 2. ) ) );
      // Gaussian
      //return - K * exp( - hzeta * hzeta );
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
  cout << "--------------------- Boundary region ---------------------" << endl;

  double beta( 0.0 );             // Hartree parameter
  double zeta0( 1.0 );            // Slot width
  double K( Example::K );         // Injection rate
  int I( 10 );                    // Number of x collocation points
  int J( I );                     // Number of y collocation points ( J = I )
  int N_BASE( 60 );               // Number of coefficients in the base solution
  double L( 6.0 );                  // Map parameter
  double eps2( Example::EPS * Example::EPS );
  double tol( 1e-10 );              // Tolerance correction coefficients norm
  int max_iter( 20 );             // Maximum number of iterations
  int size( I * J );
  int n( 250 );                     // Number of output points (x and y)
  double out_max( 10.0 );

  Timer total_timer;
  total_timer.start();

  cout << " * I = J = " << I << endl;
  cout << " * L = " << L << endl;

  Vector<double> x, y, x_grid, y_grid;
  x.rational_semi_grid( I - 1, L );
  y.rational_semi_grid( J - 1, L );
  x.push_back( 0.0 );
  y.push_back( 0.0 );

  x_grid.linspace( 0.0, out_max, n );
  y_grid.linspace( 0.0, out_max, n );
  Mesh2D<double> solution( x_grid, y_grid, 1 );

  Matrix<double> mat( 4 * size, 4 * size, 0.0 );
  Vector<double> F( 4 * size ), a_c( 4 * size );
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
  Spectral2D<double> Phi_g( c, I, J, "RationalSemi", L );
  Spectral2D<double> Psi_g( c, I, J, "RationalSemi", L );
  Spectral2D<double> U_g( c, I, J, "RationalSemi", L );
  Spectral2D<double> Theta_g( c, I, J, "RationalSemi", L );

  // Solve the Falkner-Skan equation for the base flow
  cout << " * Solving the Falkner-Skan equation using " << N_BASE << endl
       << " * spectral coefficients." << endl;
  //TODO should this be timed too?
  Spectral<double> u_base;
  Falkner<double> falkner( beta );
  u_base = falkner.solve( N_BASE, L, tol, Example::base_guess, max_iter );
  cout << " * f''(0) = " << u_base( 0.0, 2 ) << endl;

  // Put base solution into a mesh at the collocation points for later use
  Mesh1D<double> base( y, 6 );
  for ( std::size_t j = 0; j < J; j++ )
  {
    double eta( y[ j ] );
    base( j, UB ) = u_base( eta, 1 ) + 1;           // UB = f' = u' + 1
    base( j, UBd ) = u_base( eta, 2 );              // UB' = f'' = u''
    base( j, PhiB ) = u_base( eta ) + eta;          // PhiB = f = u + eta
    base( j, ThetaB ) = ( 1. - beta ) * u_base( eta, 2 );  // ThetaB = (1-beta)u''
    base( j, ThetaBd ) = ( 1. - beta ) * u_base( eta, 3 ); // ThetaB = (1-beta)u'''
    base( j, PsiB ) = ( 1. - beta ) * ( u_base( eta, 1 ) + 1. ); // PsiB = (1-beta)*(u' + 1)
  }

  double xi, yj;
  int f, g, ix, jy;

  double norm;
  max_iter = 10;
  int iter( 0 );

  Vector<double> TLx, TLy;

  Timer timer;

  do{

    cout << "---------------------------------------------------------" << endl;
    timer.start();

    // Make mesh for storing U_g, U_g_x, U_g_y, U_g_xx, U_g_yy, U_g_xy
    Mesh2D<double> U_g_mesh = U_g.mesh_derivatives( x, y );
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

        TLx = rationalsemi.eval_2( xi, f );  // TLx, TLx' and TLx''
        TLy = rationalsemi.eval_2( yj, g );  // TLy, TLy' and TLy''

        // x = 0 boundary u_c = 1 - U_g (left)
        if ( i == 0 )
        {
          mat( Phi * size + M, Phi * size + N ) = rationalsemi( 0.0, f ) * TLy[ 0 ];
          mat( Psi * size + M, Psi * size + N ) = rationalsemi( 0.0, f ) * TLy[ 0 ];
          mat( U * size + M, U * size + N ) = rationalsemi( 0.0, f ) * TLy[ 0 ];
          mat( Theta * size + M, Theta * size + N ) = rationalsemi( 0.0, f ) * TLy[ 0 ];
        }
        // y = 0 boundary u_c = 1 - U_g (bottom)
        if ( j == 0 && i != 0 )
        {
          mat( Phi * size + M, Phi * size + N ) = TLx[ 0 ] * rationalsemi( 0.0, g );
          mat( Psi * size + M, Psi * size + N ) = TLx[ 0 ] * rationalsemi( 0.0, g );
          mat( U * size + M, U * size + N ) = TLx[ 0 ] * rationalsemi( 0.0, g );
          mat( Theta * size + M, Theta * size + N ) = TLx[ 0 ] * rationalsemi( 0.0, g );
        }
        // Internal nodes
        if ( i != 0 && j != 0 )
        {
          // u_c_xx
          mat( Phi * size + M, Phi * size + N )  = TLx[ 2 ] * TLy[ 0 ];
          mat( Psi * size + M, Psi * size + N )  = TLx[ 2 ] * TLy[ 0 ];
          mat( U * size + M, U * size + N )  = TLx[ 2 ] * TLy[ 0 ];
          mat( Theta * size + M, Theta * size + N )  = TLx[ 2 ] * TLy[ 0 ];
          // u_c_yy
          mat( Phi * size + M, Phi * size + N ) += TLx[ 0 ] * TLy[ 2 ];
          mat( Psi * size + M, Psi * size + N ) += TLx[ 0 ] * TLy[ 2 ];
          mat( U * size + M, U * size + N ) += TLx[ 0 ] * TLy[ 2 ];
          mat( Theta * size + M, Theta * size + N ) += TLx[ 0 ] * TLy[ 2 ];
          // U_g_y * u_c_x
          mat( Phi * size + M, Phi * size + N )     += U_g_mesh( ix, jy, 2 )
                                                    * TLx[ 1 ] * TLy[ 0 ];
          mat( Psi * size + M, Psi * size + N )     += U_g_mesh( ix, jy, 2 )
                                                    * TLx[ 1 ] * TLy[ 0 ];
          mat( U * size + M, U * size + N )         += U_g_mesh( ix, jy, 2 )
                                                    * TLx[ 1 ] * TLy[ 0 ];
          mat( Theta * size + M, Theta * size + N ) += U_g_mesh( ix, jy, 2 )
                                                    * TLx[ 1 ] * TLy[ 0 ];

          // U_g_x * u_c_y
          mat( Phi * size + M, Phi * size + N ) += U_g_mesh( ix, jy, 1 )
                         * TLx[ 0 ] * TLy[ 1 ];
          mat( Psi * size + M, Psi * size + N ) += U_g_mesh( ix, jy, 1 )
                         * TLx[ 0 ] * TLy[ 1 ];
          mat( U * size + M, U * size + N )         += U_g_mesh( ix, jy, 1 )
                         * TLx[ 0 ] * TLy[ 1 ];
          mat( Theta * size + M, Theta * size + N ) += U_g_mesh( ix, jy, 1 )
                         * TLx[ 0 ] * TLy[ 1 ];
          // - eps^2 * a * u_c
          mat( Phi * size + M, Phi * size + N ) += - eps2 * ( xi * xi + yj * yj )
                         * TLx[ 0 ] * TLy[ 0 ];
          mat( Psi * size + M, Psi * size + N ) += - eps2 * ( xi * xi + yj * yj )
                         * TLx[ 0 ] * TLy[ 0 ];
          mat( U * size + M, U * size + N ) += - eps2 * ( xi * xi + yj * yj )
                         * TLx[ 0 ] * TLy[ 0 ];
          mat( Theta * size + M, Theta * size + N ) += - eps2 * ( xi * xi + yj * yj )
                         * TLx[ 0 ] * TLy[ 0 ];
          // - 2 * eps^2 * b * U_g * u_c
          mat( Phi * size + M, Phi * size + N ) += - 2 * eps2 * xi * yj * U_g_mesh( ix, jy, 0 )
                         * TLx[ 0 ] * TLy[ 0 ];
          mat( Psi * size + M, Psi * size + N ) += - 2 * eps2 * xi * yj * U_g_mesh( ix, jy, 0 )
                         * TLx[ 0 ] * TLy[ 0 ];
          mat( U * size + M, U * size + N ) += - 2 * eps2 * xi * yj * U_g_mesh( ix, jy, 0 )
                         * TLx[ 0 ] * TLy[ 0 ];
          mat( Theta * size + M, Theta * size + N ) += - 2 * eps2 * xi * yj * U_g_mesh( ix, jy, 0 )
                         * TLx[ 0 ] * TLy[ 0 ];
        }
      }

      // Residuals
      // x = 0 boundary u_c = 1 - U_g (left)
      if ( i == 0 )
      {
        F[ Phi * size + M ] = 1.0 - U_g_mesh( I - 1, jy, 0 );
        F[ Psi * size + M ] = 1.0 - U_g_mesh( I - 1, jy, 0 );
        F[ U * size + M ] = 1.0 - U_g_mesh( I - 1, jy, 0 );
        F[ Theta * size + M ] = 1.0 - U_g_mesh( I - 1, jy, 0 );
      }
      // y = 0 boundary u_c = 1 - U_g (bottom)
      if ( j == 0 && i != 0 )
      {
        F[ Phi * size + M ] = 1.0 - U_g_mesh( ix, J - 1, 0 );
        F[ Psi * size + M ] = 1.0 - U_g_mesh( ix, J - 1, 0 );
        F[ U * size + M ] = 1.0 - U_g_mesh( ix, J - 1, 0 );
        F[ Theta * size + M ] = 1.0 - U_g_mesh( ix, J - 1, 0 );
      }
      // Internal nodes
      if ( i != 0 && j != 0 )
      {
        // - U_g_xx - U_g_yy - U_g_x * U_g_y
        // + eps^2 * U_g * ( x^2 + y^2 + x * y * U_g )
        F[ Phi * size + M ] = - U_g_mesh( ix, jy, 3 ) - U_g_mesh( ix, jy, 4 )
                              - U_g_mesh( ix, jy, 1 ) * U_g_mesh( ix, jy, 2 )
                              + eps2 * U_g_mesh( ix, jy, 0 ) * ( xi * xi + yj
                              * yj + xi * yj * U_g_mesh( ix, jy, 0 ) );
        F[ Psi * size + M ] = - U_g_mesh( ix, jy, 3 ) - U_g_mesh( ix, jy, 4 )
                              - U_g_mesh( ix, jy, 1 ) * U_g_mesh( ix, jy, 2 )
                              + eps2 * U_g_mesh( ix, jy, 0 ) * ( xi * xi + yj
                              * yj + xi * yj * U_g_mesh( ix, jy, 0 ) );
        F[ U * size + M ] = - U_g_mesh( ix, jy, 3 ) - U_g_mesh( ix, jy, 4 )
                              - U_g_mesh( ix, jy, 1 ) * U_g_mesh( ix, jy, 2 )
                              + eps2 * U_g_mesh( ix, jy, 0 ) * ( xi * xi + yj
                              * yj + xi * yj * U_g_mesh( ix, jy, 0 ) );
        F[ Theta * size + M ] = - U_g_mesh( ix, jy, 3 ) - U_g_mesh( ix, jy, 4 )
                              - U_g_mesh( ix, jy, 1 ) * U_g_mesh( ix, jy, 2 )
                              + eps2 * U_g_mesh( ix, jy, 0 ) * ( xi * xi + yj
                              * yj + xi * yj * U_g_mesh( ix, jy, 0 ) );
      }

      // BCs at infinity are natural boundary conditions

    }
    timer.print( "Matrix construction time" );
    // Solve the system for the spectral coefficients
    //cout << " mat = " << mat << endl;
    a_c = mat.solve_parallel( F );
    Vector<double> a_phi, a_psi, a_u, a_theta;
    for ( std::size_t i = 0; i < size; i++ )
    {
      a_phi.push_back( a_c[ Phi * size + i ] );
      a_psi.push_back( a_c[ Psi * size + i ] );
      a_u.push_back( a_c[ U * size + i ] );
      a_theta.push_back( a_c[ Theta * size + i ] );
    }
    U_g.update_coefficients( a_u );
    norm = a_c.norm_2();
    cout << "  * iter = " << iter << ", norm = " << std::scientific << norm
         << std::fixed << endl;
    ++iter;
    timer.print( "Total time" );
    timer.stop();

  }while( norm > tol && iter < max_iter );

  cout << "---------------------------------------------------------" << endl;

  // Output the base solution to a 1D mesh
  Mesh1D<double> base_solution( y_grid, 6 );
  for ( std::size_t j = 0; j < n; j++ )
  {
    double eta( y_grid[ j ] );
    // UB = f' = u' + 1
    base_solution( j, UB )      = u_base( eta, 1 ) + 1;
    // UB' = f'' = u''
    base_solution( j, UBd )     = u_base( eta, 2 );
    // PhiB = f = u + eta
    base_solution( j, PhiB )    = u_base( eta ) + eta;
    // ThetaB = ( 1 - beta ) * f'' = ( 1 - beta ) * u''
    base_solution( j, ThetaB )  = ( 1. - beta ) * u_base( eta, 2 );
    // ThetaBd = ( 1 - beta ) * f''' = ( 1 - beta ) * u'''
    base_solution( j, ThetaBd ) = ( 1. - beta ) * u_base( eta, 3 );
    // PsiB = ( 1 - beta ) * f' = ( 1 - beta ) * ( u' + 1 )
    base_solution( j, PsiB )    = ( 1. - beta ) * ( u_base( eta, 1 ) + 1. );
  }
  base_solution.output( "./DATA/Boundary_region_base.dat" );

  // Output the spectral solution to a 2D mesh
  U_g( solution, 0 );


  // Output the mesh to a file
  solution.output( "./DATA/Boundary_region.dat" );

  total_timer.print();
  total_timer.stop();

  cout << "--- FINISHED ---" << endl;
}
