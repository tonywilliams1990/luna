/// \file Falkner_3D_alternative.cpp
/// \ingroup Examples
/// \ingroup Spectral
/// \todo TODO Edit description + figure out how to converge to 3D solution more
/// accurately.

#include "Luna/Core"
#include "Luna/ODE"
#include "Luna/Spectral"

namespace Luna
{
  namespace Example
  {
    double initial_guess_H( const double& y )
    {
      /// \todo need K parameter
      return - y * exp( - y );
    }

    double initial_guess_G( const double& y )
    {
      //return 2. + ( y * y + 2 * y + 2. ) * exp( - y );
      //return 0.0;
      return 0.35 * ( 1.0 - exp( - y ) );
      //return y + exp( - y ) - std::erf( 1. / y ) - y * exp( - 1. / ( y * y ) );
    }
  } // End of namespace Example
} // End of namespace Luna

using namespace std;
using namespace Luna;
using namespace Example;

int main()
{
  cout << "------------------ Falkner_3D_alternative -----------------" << endl;

  //cout << " * The Falkner-Skan equation, for a boundary layer on flat " << endl
  //     << " * plate, is solved using a rational Chebyshev spectral " << endl
  //     << " * method." << endl;

  double beta( 0.0 );               // Hartree parameter
  /// \todo K parameter
  int n( 40 );                      // Number of spectral coefficients
  double L( 5.0 );                  // Map parameter
  double tol( 1e-10 );              // Tolerance correction coefficients norm
  Vector<double> y;                 // Vector for collocation grid
  y.rational_semi_grid( n, L );
  y.push_back( 0.0 );

  int N( 1000 );                        // Number of output points
  Vector<double> grid;                  // Output grid
  grid.powspace( 0.0, y[ 0 ], N, 3 );   // Squash nodes near to the boundary
  Mesh1D<double> solution( grid, 6 );   // Output mesh


  RationalSemi<double> rationalsemi( L );
  Vector<double> H_g_coeffs( n, 0.0 );           // Vector of coefficients
  Vector<double> G_g_coeffs( n, 0.0 );

  // Approximate the guess as a semi-infinite rational Chebyshev polynomial
  H_g_coeffs = rationalsemi.approximate( Example::initial_guess_H, n );
  G_g_coeffs = rationalsemi.approximate( Example::initial_guess_G, n );
  // Extra coefficients for boundary conditions
  for ( std::size_t i = 0; i < 2; i++ )
  {
    H_g_coeffs.push_back( 0.0 );
    G_g_coeffs.push_back( 0.0 );
  }
  //TODO cout << " H_g_coeffs = " << H_g_coeffs << endl;
  //TODO cout << " G_g_coeffs = " << G_g_coeffs << endl;
  // Make this into a spectral solution
  Spectral<double> H_g( H_g_coeffs, "RationalSemi", L );
  Spectral<double> G_g( G_g_coeffs, "RationalSemi", L );

  double norm;
  int max_iter( 100 );
  int iter( 0 );
  int size( 2 * n + 4 );

  Timer timer;
  timer.start();

  do
  {
    // Setup the system of equations to find the correction
    Matrix<double> M( size, size, 0.0 );
    Vector<double> F( size ), a_c( size );
    Vector<double> H_a_c, G_a_c;
    int row( 0 );

    // Don't need the BC at infinity as this is a "natural" BC

    // H equation
    for ( std::size_t i = 0; i < n; i++ )
    {
      double yi = y[ i ];
      for ( std::size_t j = 0; j < n + 2; j++ )
      {
        // H_c'''
        M( row, j )  = rationalsemi( yi, j, 3 );
        // (y + H_g + ( 2 - beta ) * G_g ) * H_c''
        M( row, j ) += ( yi + H_g( yi ) + ( 2. - beta ) * G_g( yi ) )
                       * rationalsemi( yi, j, 2 );
        // - 2 * beta * (1 + H_g') * H_c'
        M( row, j ) -= 2 * beta * ( 1 + H_g( yi, 1 ) )
                       * rationalsemi( yi, j, 1 );
        // H_g'' * H_c
        M( row, j ) += H_g( yi, 2 ) * rationalsemi( yi, j );
        // ( 2 - beta ) * H_g'' * G_c
        M( row, n + 2 + j ) += ( 2. - beta ) * H_g( yi, 2 )
                             * rationalsemi( yi, j );
      }
      // - H_g''' - ( y + H_g + ( 2 - beta ) * G_g ) * H_g''
      // + beta * H_g' ( 2 + H_g' )
      F[ row ]  = - H_g( yi, 3 ) - ( yi + H_g( yi ) + ( 2. - beta )
                  * G_g( yi ) ) * H_g( yi, 2 );
      F[ row ] += beta * H_g( yi, 1 ) * ( 2 + H_g( yi, 1 ) );

      ++row;
    }

    // y = 0 boundary H_c = - K - H_g TODO K parameter
    for ( std::size_t j = 0; j < n + 2; j++ )
    {
      double yi = 0.0;
      M( row, j ) = rationalsemi( yi, j );
    }
    F[ row ] = - H_g( 0.0 );

    ++row;

    // y = 0 boundary H_c' = -1 - H_g
    for ( std::size_t j = 0; j < n + 2; j++ )
    {
      double yi = 0.0;
      M( row, j ) = rationalsemi( yi, j, 1 );
    }
    F[ row ] = - 1.0 - H_g( 0.0, 1 );

    ++row;

    // G equation
    for ( std::size_t i = 0; i < n; i++ )
    {
      double yi = y[ i ];
      for ( std::size_t j = 0; j < n + 2; j++ )
      {
        // G_c'''
        M( row, n + 2 + j )  = rationalsemi( yi, j, 3 );
        // ( y + H_g + ( 2 - beta ) * G_g ) * G_c''
        M( row, n + 2 + j ) += ( yi + H_g( yi ) + ( 2. - beta ) * G_g( yi ) )
                               * rationalsemi( yi, j, 2 );
        // 2 * ( 1 - beta ) * ( 1 + H_g' ) * G_c'
        M( row, n + 2 + j ) += 2 * ( 1. - beta ) * ( 1. + H_g( yi, 1 ) )
                               * rationalsemi( yi, j, 1 );
        // - 2 * ( 2 - beta ) * G_g' * G_c'
        M( row, n + 2 + j ) += - 2 * ( 2. - beta ) * G_g( yi, 1 )
                               * rationalsemi( yi, j, 1 );
        // ( 2 - beta ) * G_g'' * G_c
        M( row, n + 2 + j ) += ( 2. - beta ) * G_g( yi, 2 )
                               * rationalsemi( yi, j );
        // 2 * ( 1 - beta ) * G_g' * H_c'
        M( row, j ) += 2 * ( 1. - beta ) * G_g( yi, 1 )
                               * rationalsemi( yi, j, 1 );
        // G_g'' * H_c
        M( row, j ) += G_g( yi, 2 ) * rationalsemi( yi, j );

      }
      // - G_g''' - ( y + H_g + ( 2 - beta ) * G_g ) * G_g''
      // - 2 * ( 1 - beta ) * ( 1 + H_g' ) * G_g' + ( 2 - beta ) * G_g'^2
      F[ row ]  = - G_g( yi, 3 ) - ( yi + H_g( yi ) + ( 2. - beta )
                        * G_g( yi ) ) * G_g( yi, 2 );
      F[ row ] += - 2 * ( 1. - beta ) * ( 1. + H_g( yi, 1 ) ) * G_g( yi, 1 );
      F[ row ] += ( 2. - beta ) * G_g( yi, 1 ) * G_g( yi, 1 );

      ++row;

    }

    // y = 0 boundary G_c = - G_g
    for ( std::size_t j = 0; j < n + 2; j++ )
    {
      double yi = 0.0;
      M( row, n + 2 + j ) = rationalsemi( yi, j );
    }
    F[ row ] = - G_g( 0.0 );

    ++row;

    // y = 0 boundary G_c' = - G_g'
    for ( std::size_t j = 0; j < n + 2; j++ )
    {
      double yi = 0.0;
      M( row, n + 2 + j ) = rationalsemi( yi, j, 1 );
    }
    F[ row ] = - G_g( 0.0, 1 );

    ++row;


    // Solve the system for the correction spectral coefficients
    a_c = M.solve( F );
    for ( std::size_t j = 0; j < n + 2; j++ )
    {
      H_a_c.push_back( a_c[ j ] );
      G_a_c.push_back( a_c[ n + 2 + j ] );
    }

    H_g.update_coefficients( H_a_c );
    G_g.update_coefficients( G_a_c );
    norm = a_c.norm_2();
    ++iter;

  }while( norm > tol && iter < max_iter );

  cout << " * The final spectral coefficient for n = " << n << " is: "
       << std::scientific << H_g.get_coefficients()[ n-1 ] << endl;
  cout << " * F''(0) = " << std::fixed << H_g( 0.0, 2 ) << endl;

  // Put spectral solution into output mesh
  for ( std::size_t i = 0; i < N; i++ )
  {
    double y( grid[ i ] );
    solution( i, 0 ) = y + H_g( y );              // F   = eta + H
    solution( i, 1 ) = 1. + H_g( y, 1 );          // F'  = 1 + H'
    solution( i, 2 ) = H_g( y, 2 );               // F'' = H''
    solution( i, 3 ) = G_g( y );                  // G
    solution( i, 4 ) = G_g( y, 1 );               // G'
    solution( i, 5 ) = G_g( y, 2 );               // G''
  }

  solution.output( "./DATA/Falkner_3D_alternative.dat" );

  timer.print();
  timer.stop();

  cout << " * To see the spectral solution and its derivatives run: " << endl;
  cout << "python Plotting/Falkner_3D_alternative_plot.py" << endl;

  cout << "--- FINISHED ---" << endl;
}
