/// \file  Heat_equation.cpp
/// \ingroup Examples
/// Solve the heat equation  \f[ u_t = \alpha u_{xx}, \f] where \f$ \alpha \f$
/// is the thermal diffusivity, subject to the boudary conditions
/// \f$ u(0,t)=u(1,t)=0 \f$ and the initial condition
/// \f[ u(x,0) = u_0(x) = 100 \sin( \pi x ). \f]
/// The numerical solution is found using the Crank-Nicolson method.

#include "Luna/Core"
#include "Luna/Sparse"

using namespace std;
using namespace Luna;


int main()
{
  cout << "---------------------- Heat equation ----------------------" << endl;

  double T_max( 1.0 );                // Maximum simulation time
  int N( 200 );                       // Number of time steps
  int J( 50 );                        // Number of spatial steps
  double dt( T_max / ( 1. * N ) );    // Time step
  double dx( 1.0 / J );               // Spatial step
  double alpha( 1.0 );                // Thermal diffusivity

  cout << " * Solving the heat equation u_t = a u_{xx} in the domain" << endl
       << " * t=[0," + Utility::stringify(T_max) + "] and x=[0,1] subject to "
       << " the boundary conditions " << endl << " * u(0,t)=u(1,t)=0 and the "
       << "initial condition u(x,0)=u_0(x)." << endl << " * The initial "
       << "condition u_0(x) is specified in the code."
       << endl ;

  cout << " * The spatial and time steps are:" << endl;
  cout << " * dx = " << dx << endl;
  cout << " * dt = " << dt << endl;
  cout << " * and the number of spatial and time steps are:" << endl;
  cout << " * J = " << J << endl;
  cout << " * N = " << N << endl;

  double mu( alpha * dt / ( dx * dx ) );

  // Mesh for storing the solution
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
  Explicit( 0, 0 )     =   1.0;
  Explicit( 0, 1 )     =   0.0;
  Explicit( J, J - 1 ) =   0.0;
  Explicit( J, J )     = - 1.0;

  // Initial condition u_0(x) = 100 * sin( pi * x )
  Vector<double> current( J + 1, 0.0 );
  current[ 0 ] = 0.0;
  solution( 0, 0, 0 ) = current[ 0 ];
  for ( std::size_t j = 1; j < J; j++ )
  {
    current[ j ] = 100 * sin( M_PI * j * dx );
    solution( j, 0, 0 ) = current[ j ];
  }
  current[ J ] = 0.0;
  solution( J, 0, 0 ) = current[ J ];

  Vector<double> next( J + 1 );

  cout << " * Calculating the solution..." << endl;

  for ( std::size_t i = 1; i < N + 1; i++ )
  {
    current = Explicit * current;
    next = Implicit.solve( current );
    current = next;
    for ( std::size_t j = 0; j < J + 1; j++ )
    {
      solution( j, i, 0 ) = current[ j ];
    }
  }
  cout << " * Solution complete." << endl;

  string outstring( "./DATA/Heat_equation" );
  outstring += "_J_" + Utility::stringify( J );
  outstring += "_N_" + Utility::stringify( N );
  outstring += ".dat";
  solution.output( outstring );

  cout << " * For an animation of the solution run:" << endl;
  cout << "python Plotting/Heat_plot.py " << J << " " << N << endl;
  cout << "--- FINISHED ---" << endl;
}
