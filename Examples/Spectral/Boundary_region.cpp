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
      //return - K * 0.5 * ( tanh( GAMMA * ( hzeta - 1. ) )
      //       - tanh( GAMMA * ( hzeta - 2. ) ) );
      // Gaussian
      return - K * exp( - hzeta * hzeta );
    }

    double Phi_w_hzeta_func( const double& hzeta ){
      // Gaussian
      return 2. * K * hzeta * exp( - hzeta * hzeta );
    }

    double initial_guess_x( const double& x )
    {
      //return exp( - EPS * x );
      return x * x * exp( - x );
    }

    double initial_guess_y( const double& y )
    {
      //return exp( - EPS * y );
      return y * exp( - y );
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
  int I( 25 );                    // Number of x collocation points
  int J( I );                     // Number of y collocation points ( J = I )
  int N_BASE( 60 );               // Number of coefficients in the base solution
  double L( 6.0 );                // Map parameter
  double L_BASE( 6.0 );           // Map parameter for Falkner-Skan
  double eps2( Example::EPS * Example::EPS );
  double tol( 1e-10 );            // Tolerance correction coefficients norm
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
  Vector<double> guess( size, 0.0 );
  Spectral2D<double> Phi_g( c, I, J, "RationalSemi", L );
  Spectral2D<double> Psi_g( c, I, J, "RationalSemi", L );
  Spectral2D<double> U_g( c, I, J, "RationalSemi", L );
  Spectral2D<double> Theta_g( c, I, J, "RationalSemi", L );

  //Theta_g( solution, 0 );

  // Solve the Falkner-Skan equation for the base flow
  cout << " * Solving the Falkner-Skan equation using " << N_BASE << endl
       << " * spectral coefficients." << endl;
  //TODO should this be timed too?
  Spectral<double> u_base;
  Falkner<double> falkner( beta );
  u_base = falkner.solve( N_BASE, L_BASE, tol, Example::base_guess, max_iter );
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

    // Make mesh for storing Var, Var_x, Var_y, Var_xx, Var_yy, Var_xy
    Mesh2D<double> Phi_g_mesh = Phi_g.mesh_derivatives( x, y );
    Mesh2D<double> Psi_g_mesh = Psi_g.mesh_derivatives( x, y );
    Mesh2D<double> U_g_mesh = U_g.mesh_derivatives( x, y );
    Mesh2D<double> Theta_g_mesh = Theta_g.mesh_derivatives( x, y );
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

        // hzeta = 0 boundary
        if ( i == 0 )
        {
          // Phi_hzeta = 0
          mat( Phi * size + M, Phi * size + N ) = TLx[ 1 ] * TLy[ 0 ];
          // Psi = 0
          mat( Psi * size + M, Psi * size + N ) = TLx[ 0 ] * TLy[ 0 ];
          // U_hzeta = 0
          //mat( U * size + M, U * size + N ) = rationalsemi( 0.0, f ) * TLy[ 0 ];
          mat( U * size + M, U * size + N ) = TLx[ 1 ] * TLy[ 0 ];
          // Theta = 0
          mat( Theta * size + M, Theta * size + N ) = TLx[ 0 ] * TLy[ 0 ];
        }
        // eta = 0 boundary
        if ( j == 0 && i != 0 )
        {
          // Phi = Phiw
          mat( Phi * size + M, Phi * size + N ) = TLx[ 0 ] * TLy[ 0 ];
          // Psi = 0
          mat( Psi * size + M, Psi * size + N ) = TLx[ 0 ] * TLy[ 0 ];
          // U = 0
          //mat( U * size + M, U * size + N ) = TLx[ 0 ] * rationalsemi( 0.0, g );
          mat( U * size + M, U * size + N ) = TLx[ 0 ] * TLy[ 0 ];
          // Theta = Psi_hzeta - ( 1 / zeta0^2 ) * Phiw_hzeta TODO
          /*mat( Theta * size + M, Theta * size + N ) = TLx[ 0 ] * TLy[ 0 ];
          mat( Theta * size + M, Psi * size + N ) = - TLx[ 0 ] * TLy[ 1 ];*/
          // Theta = 0 (temporary)
          mat( Theta * size + M, Theta * size + N ) = TLx[ 0 ] * TLy[ 0 ];
        }
        // Internal nodes
        if ( i != 0 && j != 0 )
        {
          // Phi equation
          // Phi_c_yy
          mat( Phi * size + M, Phi * size + N )  = TLx[ 0 ] * TLy[ 2 ];
          // + Phi_c_xx / zeta0^2
          mat( Phi * size + M, Phi * size + N ) += TLx[ 2 ] * TLy[ 0 ]
                                                   / ( zeta0 * zeta0 );
          // - ( 2 - beta ) * U_c_y
          mat( Phi * size + M, U * size + N ) += - ( 2. - beta )
                                                 * TLx[ 0 ] * TLy[ 1 ];
          // + Theta_c_x
          mat( Phi * size + M, Theta * size + N ) += TLx[ 1 ] * TLy[ 0 ];

          // Psi equation
          // Psi_c_yy
          mat( Psi * size + M, Psi * size + N )  = TLx[ 0 ] * TLy[ 2 ];
          // + Psi_c_xx / zeta0^2
          mat( Psi * size + M, Psi * size + N ) += TLx[ 2 ] * TLy[ 0 ]
                                                   / ( zeta0 * zeta0 );
          // - ( 2 - beta ) * U_c_x / zeta0^2
          mat( Psi * size + M, U * size + N ) += - ( 2. - beta )
                                                 * TLx[ 1 ] * TLy[ 0 ]
                                                 / ( zeta0 * zeta0 );
          // - Theta_c_y
          mat( Psi * size + M, Theta * size + N ) += TLx[ 0 ] * TLy[ 1 ];

          // U equation
          // U_c_yy
          mat( U * size + M, U * size + N )  = TLx[ 0 ] * TLy[ 2 ];
          // + U_c_xx / zeta0^2
          mat( U * size + M, U * size + N ) += TLx[ 2 ] * TLy[ 0 ]
                                               / ( zeta0 * zeta0 );
          // - 2 * beta * ( UB + U_g ) * U_c
          mat( U * size + M, U * size + N ) -= 2 * beta * ( base( jy, UB )
                                               + U_g_mesh( ix, jy, 0 ) )
                                               * TLx[ 0 ] * TLy[ 0 ];
          // + ( hzeta * PsiB + Psi_g ) * U_c_x
          //mat( U * size + M, U * size + N ) += ( xi * base( jy, PsiB )
          //                                    + Psi_g_mesh( ix, jy, 0 ) )
          //                                    * TLx[ 1 ] * TLy[ 0 ];
          // + hzeta * PsiB * U_c_x
          mat( U * size + M, U * size + N ) += xi * base( jy, PsiB )
                                              * TLx[ 1 ] * TLy[ 0 ];
          // + Psi_g * U_c_x TODO
          //mat( U * size + M, U * size + N ) += Psi_g_mesh( ix, jy, 0 )
          //                                    * TLx[ 1 ] * TLy[ 0 ];
          // + U_g_x * Psi_c TODO
          //mat( U * size + M, Psi * size + N ) += U_g_mesh( ix, jy, 1 )
          //                                     * TLx[ 0 ] * TLy[ 0 ];
          // + PhiB * U_c_y
          mat( U * size + M, U * size + N ) += base( jy, PhiB )
                                               * TLx[ 0 ] * TLy[ 1 ];
          // + Phi_g * U_c_y TODO
          //mat( U * size + M, U * size + N ) += Phi_g_mesh( ix, jy, 0 )
          //                                    * TLx[ 0 ] * TLy[ 1 ];
          // + ( PhiB + Phi_g ) * U_c_y
          /*mat( U * size + M, U * size + N ) += ( base( jy, PhiB )
                                               + Phi_g_mesh( ix, jy, 0 ) )
                                               * TLx[ 0 ] * TLy[ 1 ];*/
          // + U_g_x * Psi_c
          //mat( U * size + M, Psi * size + N ) += U_g_mesh( ix, jy, 1 )
          //                                     * TLx[ 0 ] * TLy[ 0 ];
          // + UBd  * Phi_c
          mat( U * size + M, Phi * size + N ) += base( jy, UBd )
                                               * TLx[ 0 ] * TLy[ 0 ];
          // + U_g_y * Phi_c TODO
          //mat( U * size + M, Phi * size + N ) += U_g_mesh( ix, jy, 2 )
          //                                     * TLx[ 0 ] * TLy[ 0 ];
          // + ( UBd + U_g_y ) * Phi_c
          /*mat( U * size + M, Phi * size + N ) += ( base( jy, UBd )
                                               + U_g_mesh( ix, jy, 2 ) )
                                               * TLx[ 0 ] * TLy[ 0 ];
          // + Phi_g * U_c_y
          mat( U * size + M, U * size + N ) += Phi_g_mesh( ix, jy, 0 )
                                               * TLx[ 0 ] * TLy[ 1 ];*/

          // Theta equation
          // Theta_c_yy
          mat( Theta * size + M, Theta * size + N )  = TLx[ 0 ] * TLy[ 2 ];
          // + Theta_c_xx / zeta0^2
          mat( Theta * size + M, Theta * size + N ) += TLx[ 2 ] * TLy[ 0 ]
                                                       / ( zeta0 * zeta0 );
          /*// - 2 * ( 1 - beta ) * hzeta * ( UB + U_g ) * U_c_y
          mat( Theta * size + M, U * size + N ) -= 2 * ( 1. - beta ) * xi
                                                  * ( base( jy, UB )
                                                    + U_g_mesh( ix, jy, 0 ) )
                                                  * TLx[ 0 ] * TLy[ 1 ];
          // + 2 * ( 1 - beta ) * eta * ( UB + U_g ) * U_c_x / zeta0^2
          mat( Theta * size + M, U * size + N ) += 2 * ( 1. - beta ) * yj
                                                  * ( base( jy, UB )
                                                    + U_g_mesh( ix, jy, 0 ) )
                                                  * TLx[ 1 ] * TLy[ 0 ]
                                                  / ( zeta0 * zeta0 );
          // - 2 * ( 1 - beta ) * hzeta * ( UBd + U_g_y ) * U_c
          mat( Theta * size + M, U * size + N ) -= 2 * ( 1. - beta ) * xi
                                                  * ( base( jy, UBd )
                                                    + U_g_mesh( ix, jy, 2 ) )
                                                  * TLx[ 0 ] * TLy[ 0 ];
          // + 2 * ( 1 - beta ) * eta * U_g_x * U_c / zeta0^2
          mat( Theta * size + M, U * size + N ) += 2 * ( 1. - beta ) * yj
                                                  * U_g_mesh( ix, jy, 1 )
                                                  * TLx[ 0 ] * TLy[ 0 ]
                                                  / ( zeta0 * zeta0 );
          // + ( PhiB + Phi_g ) * Theta_c_y
          mat( Theta * size + M, Theta * size + N ) += ( base( jy, PhiB )
                                                    + Phi_g_mesh( ix, jy, 0 ) )
                                                    * TLx[ 0 ] * TLy[ 1 ];
          // + ( hzeta * ThetaBd + Theta_g_y ) * Phi_c
          mat( Theta * size + M, Phi * size + N ) += ( xi * base( jy, ThetaBd )
                                                  + Theta_g_mesh( ix, jy, 2 ) )
                                                  * TLx[ 0 ] * TLy[ 0 ];
          // + ( hzeta * PsiB + Psi_g ) * Theta_c_x
          mat( Theta * size + M, Theta * size + N ) += ( xi * base( jy, PsiB )
                                                    + Psi_g_mesh( ix, jy, 0 ) )
                                                    * TLx[ 1 ] * TLy[ 0 ];
          // + ( ThetaB + Theta_g_x ) * Psi_c
          mat( Theta * size + M, Psi * size + N ) += ( base( jy, ThetaB )
                                                  + Theta_g_mesh( ix, jy, 1 ) )
                                                  * TLx[ 0 ] * TLy[ 0 ];
          // + ( 2 - beta ) * ( UB + U_g ) * Theta_c
          mat( Theta * size + M, Theta * size + N ) += ( 2. - beta )
                                                    * ( base( jy, UB )
                                                    + U_g_mesh( ix, jy, 0 ) )
                                                    * TLx[ 0 ] * TLy[ 0 ];
          // + ( 2 - beta ) * ( hzeta * ThetaB + Theta_g ) * U_c
          mat( Theta * size + M, U * size + N ) += ( 2. - beta )
                                                * ( xi * base( jy, ThetaB )
                                                + Theta_g_mesh( ix, jy, 0 ) )
                                                * TLx[ 0 ] * TLy[ 0 ];*/
        }
      }

      // Residuals
      // hzeta = 0 boundary
      if ( i == 0 )
      {
        // Phi_hzeta = 0
        F[ Phi * size + M ]   = - Phi_g_mesh( I - 1, jy, 1 );
        // Psi = 0
        F[ Psi * size + M ]   = - Psi_g_mesh( I - 1, jy, 0 );
        // U_hzeta = 0
        F[ U * size + M ]     = - U_g_mesh( I - 1, jy, 1 );
        // Theta = 0
        F[ Theta * size + M ] = - Theta_g_mesh( I - 1, jy, 0 );
      }
      // eta = 0 boundary
      if ( j == 0 && i != 0 )
      {
        // Phi = Phiw
        /*F[ Phi * size + M ]   =   Example::Phi_w_func( xi )
                                - Phi_g_mesh( ix, J - 1, 0 );*/
        F[ Phi * size + M ]   = - Phi_g_mesh( ix, J - 1, 0 ); // temporary TODO

        // Psi = 0
        F[ Psi * size + M ]   = - Psi_g_mesh( ix, J - 1, 0 );
        // U = 0
        F[ U * size + M ]     = - U_g_mesh( ix, J - 1, 0 );

        // Theta = Psi_hzeta - ( 1 / zeta0^2 ) * Phiw_hzeta
        /*F[ Theta * size + M ] =   Psi_g_mesh( ix, J - 1, 1 )
                                - Theta_g_mesh( ix, J - 1, 0 )
                                - ( Example::Phi_w_hzeta_func( xi )
                                / ( zeta0 * zeta0 ) );*/
        F[ Theta * size + M ] = - Theta_g_mesh( ix, J - 1, 0 ); // temporary TODO
      }
      // Internal nodes
      if ( i != 0 && j != 0 )
      {
        // Phi equation
        F[ Phi * size + M ] = - Phi_g_mesh( ix, jy, 4 )
                              - ( Phi_g_mesh( ix, jy, 3 ) / ( zeta0 * zeta0 ) )
                              + ( 2. - beta ) * U_g_mesh( ix, jy, 2 )
                              - Theta_g_mesh( ix, jy, 1 )
        + ( xi * yj * ( xi - 2. ) + 2 * yj + 2 * xi * xi * ( yj - 2. ) ) * exp( - xi - yj ) ;

        // Psi equation
        F[ Psi * size + M ] = - Psi_g_mesh( ix, jy, 4 )
                              - ( Psi_g_mesh( ix, jy, 3 ) / ( zeta0 * zeta0 ) )
                              + ( ( 2. - beta ) * U_g_mesh( ix, jy, 1 )
                              / ( zeta0 * zeta0 ) )
                              + Theta_g_mesh( ix, jy, 2 )
        + ( 2 * ( xi * yj * ( xi - 2. ) + yj - xi * xi + ( xi - 2. ) * xi * yj ) + xi * xi * ( yj - 1. ) ) * exp( - xi - yj );

        // U equation
        F[ U * size + M ] = - U_g_mesh( ix, jy, 4 )
                            - ( U_g_mesh( ix, jy, 3 ) / ( zeta0 * zeta0 ) )
                            + beta * ( 2 * base( jy, UB )
                            + U_g_mesh( ix, jy, 0 ) ) * U_g_mesh( ix, jy, 0 )
                            - xi * base( jy, PsiB ) * U_g_mesh( ix, jy, 1 )
                            //  - Psi_g_mesh( ix, jy, 0 ) * U_g_mesh( ix, jy, 1 )
                            - base( jy, PhiB ) * U_g_mesh( ix, jy, 2 )
                            - base( jy, UBd ) * Phi_g_mesh( ix, jy, 0 )
                            //  - U_g_mesh( ix, jy, 2 ) * Phi_g_mesh( ix, jy, 0 )

        + 2 * ( xi * yj * ( xi - 2. ) + yj - xi * xi ) * exp( - xi - yj )
        - beta * ( 2 * base( jy, UB ) + xi * xi * yj * exp( - xi - yj ) ) * xi * xi * yj * exp( - xi - yj )
        + base( jy, PsiB ) * ( 2. - xi ) * xi * xi * yj * exp( - xi - yj )
        // TODO psi * u_x
        //+ xi * xi * yj * exp( - xi - yj ) * ( 2. - xi ) * xi * yj * exp( - xi - yj )
        + base( jy, PhiB ) * xi * xi * ( 1. - yj ) * exp( - xi - yj )
        + base( jy, UBd ) * xi * xi * yj * exp( - xi - yj );
        // TODO u_y * phi
        //+ xi * xi * ( 1. - yj ) * exp( - xi - yj ) * xi * xi * yj * exp( - xi - yj );


        // Theta equation
        F[ Theta * size + M ] = - Theta_g_mesh( ix, jy, 4 )
                              - ( Theta_g_mesh( ix, jy, 3 ) / ( zeta0 * zeta0 ) )

                            /*+ 2 * ( 1. - beta ) * (
        xi * ( base( jy, UB ) + U_g_mesh( ix, jy, 0 ) ) * U_g_mesh( ix, jy, 2 )
        + xi * base( jy, UBd ) * U_g_mesh( ix, jy, 0 )
        - ( yj * ( base( jy, UB ) + U_g_mesh( ix, jy, 0 ) )
        * U_g_mesh( ix, jy, 1 ) / ( zeta0 * zeta0 ) ) )
        - ( base( jy, PhiB ) + Phi_g_mesh( ix, jy, 0 ) )
        * Theta_g_mesh( ix, jy, 2 )
        - xi * base( jy, ThetaBd ) * Phi_g_mesh( ix, jy, 0 )
        - xi * base( jy, PsiB ) * Theta_g_mesh( ix, jy, 1 )
        - Psi_g_mesh( ix, jy, 0 ) * ( base( jy, ThetaB )
        + Theta_g_mesh( ix, jy, 1 ) )
        - ( 2. - beta ) * ( ( base( jy, UB ) + U_g_mesh( ix, jy, 0 ) )
        * Theta_g_mesh( ix, jy, 0 )
        + xi * base( jy, ThetaB ) * U_g_mesh( ix, jy, 0 ) );*/

        + 2 * ( xi * yj * ( xi - 2. ) + yj - xi * xi ) * exp( - xi - yj ) ;
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
    Phi_g.update_coefficients( a_phi );
    Psi_g.update_coefficients( a_psi );
    U_g.update_coefficients( a_u );
    Theta_g.update_coefficients( a_theta );

    //Vector<double> a_diff;
    //a_diff = a_psi - a_theta;
    //cout << " a_diff.norm_2() = " << std::scientific << a_diff.norm_2() << endl;

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
  Theta_g( solution, 0 );
  /*for ( std::size_t i = 0; i < x_grid.size(); i++ )
  {
    double x( x_grid[ i ] );
    for ( std::size_t j = 0; j < y_grid.size(); j++ )
    {
      double y( y_grid[ j ] );
      for ( std::size_t n = 0; n < size; n++ )
      {
        f = n / J;
        g = n % J;
        //solution( i, j, 0 ) += U_g.get_coefficients()[ n ]
        //                      * rationalsemi( x, f ) * rationalsemi( y, g );
        solution( i, j, 0 ) = //base_solution( j, UB ) +
                              U_g.get_coefficients()[ n ]
                              * rationalsemi( x, f ) * rationalsemi( y, g );
      }
    }
  }*/

  // Output the mesh to a file
  solution.output( "./DATA/Boundary_region.dat" );

  total_timer.print();
  total_timer.stop();

  cout << "--- FINISHED ---" << endl;
}
