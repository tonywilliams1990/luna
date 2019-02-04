/// \file  Heat_equation.cpp
/// \ingroup Examples
/// TODO Description + tex equations

#include "Luna/Core"
#include "Luna/Sparse"

using namespace std;
using namespace Luna;

std::string stringify( const int &val )
{
  std::stringstream temp;
  temp << val;
  return temp.str();
}

std::string stringify( const double &val, int p )
{
  std::stringstream temp;
  temp.precision( p );
  temp << val;
  return temp.str();
}

int main()
{
  cout << "----- Heat equation -----" << endl;

  double T_max( 1.0 );                // Maxium simulation time
  int N( 500 );                        // Number of time steps
  int J( 50 );                        // Number of spacial steps
  double dt( T_max / ( 1. * N ) );    // Time step
  double dx( 1.0 / J );               // Spatial step
  double alpha( 1.0 );                // Thermal diffusivity

  double mu( alpha * dt / ( dx * dx ) );

  // Mesh for storing the data
  Vector<double> t_nodes, x_nodes;
  t_nodes.linspace( 0.0, T_max, N + 1 );
  x_nodes.linspace( 0.0, 1.0, J + 1 );
  Mesh2D<double> solution( x_nodes, t_nodes, 1 );

  // Matrices for Crank-Nicolson method including BCs u(0,t) = u(1,t) = 0
  Tridiagonal<double> Implicit( -mu, 2. + 2. * mu, -mu, J + 1 );
  Implicit( 0, 0 )     = 1.0;
  Implicit( 0, 1 )     = 0.0;
  Implicit( J, J - 1 ) = 0.0;
  Implicit( J, J )     = 1.0;

  Tridiagonal<double> Explicit( mu, 2. - 2. * mu, mu, J + 1 );
  Explicit( 0, 0 )     = 1.0;
  Explicit( 0, 1 )     = 0.0;
  Explicit( J, J - 1 ) = 0.0;
  Explicit( J, J )     = 1.0;

  // Initial condition u_0(x) = 100 * sin( pi * x )
  Vector<double> current( J + 1, 0.0 );
  for ( std::size_t j = 0; j < J + 1; j++ )
  {
    current[ j ] = 100 * sin( M_PI * j * dx );
    solution( j, 0, 0 ) = current[ j ];
  }

  cout << " * dx = " << dx << endl;
  cout << " * dt = " << dt << endl;
  //cout << " * u_0(x) = " << current << endl;

  Vector<double> next( J + 1 );

  for ( std::size_t i = 1; i < N + 1; i++ )
  {
    current = Explicit * current;
    next = Implicit.solve( current );
    //cout << " * u_"<< i <<"(x) = " << next << endl;
    current = next;
    for ( std::size_t j = 0; j < J + 1; j++ )
    {
      solution( j, i, 0 ) = current[ j ];
    }
  }

  string outstring( "Heat_equation" );
  outstring += "_J_" + stringify( J );
  outstring += "_N_" + stringify( N );
  outstring += ".dat";
  solution.output( outstring );


  //TODO Utility stringify functions

  cout << "--- FINISHED ---" << endl;

}
