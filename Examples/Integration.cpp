/// \file  Integration.cpp
/// \ingroup Examples
/// Some numerical integration examples

#include "Luna/Core"

using namespace std;
using namespace Luna;

int main()
{
  cout << "----------------- Integration -----------------" << endl;
  cout << "  * Create a 1D mesh with 2 variables and set "  << endl;
  cout << "  * one equal 2x and the other to x^2. "  << endl;

  Vector<double> nodes;
  nodes.linspace( 0, 1, 101 );
  std::size_t nvars( 2 );
  Mesh1D<double> mesh( nodes, nvars );

  for ( std::size_t i = 0; i< nodes.size(); ++i )
  {
    double x( mesh.coord( i ) );
    mesh( i, 0 ) = 2 * x;
    mesh( i, 1 ) = x * x;
  }

  cout << "  * # nodes = " << mesh.nnodes() << endl;
  cout << "  * # variables = " << mesh.nvars() << endl;
  Vector<double> vars;
  vars = mesh.get_interpolated_vars( 0.35 );
  cout << "  * Interpolate the values of each of the " << endl;
  cout << "  * variables at a specified point." << endl;

  cout << "  * vars( 0.35 ) = " << vars << endl;

  cout << "  * Numerically integrate the variables over " << endl;
  cout << "  * the domain (from 0 to 1). " << endl;

  cout << "  * mesh.integral( 0 ) = " << mesh.integral( 0 ) << endl;
  cout << "  * mesh.integral( 1 ) = " << mesh.integral( 1 ) << endl;

  cout << "-----------------------------------------------" << endl;
  cout << "  * Create a 2D mesh with 2 variables and set "  << endl;
  cout << "  * one equal 2x^2 and the other to 3y. "  << endl;
  cout << "  * The domain is split uniformly from 0 to 1." << endl;

  std::size_t N( 101 );
  Vector<double> x_nodes, y_nodes;
  x_nodes.linspace( 0, 1, N );
  y_nodes.linspace( 0, 1, N );
  Mesh2D<double> mesh2d( x_nodes, y_nodes, nvars );

  for ( std::size_t i = 0; i < x_nodes.size(); ++i )
  {
    for ( std::size_t j = 0; j < y_nodes.size(); ++j )
    {
      double x( mesh2d.coord( i, j ).first );
      double y( mesh2d.coord( i, j ).second );
      mesh2d( i, j, 0 ) = 2 * x * x;
      mesh2d( i, j, 1 ) = 3 * y;
    }
  }

  cout << "  * # x nodes = " << mesh2d.get_nnodes().first << endl;
  cout << "  * # y nodes = " << mesh2d.get_nnodes().second << endl;
  cout << "  * # variables = " << mesh2d.get_nvars() << endl;
  cout << "  * Integrate each variable over the domain." << endl;
  cout << "  * mesh2d.integral2D( 0 ) = " << mesh2d.integral2D( 0 ) << endl;
  cout << "  * mesh2d.integral2D( 1 ) = " << mesh2d.integral2D( 1 ) << endl;
  cout << "-----------------------------------------------" << endl;

  std::complex<double> x( 1.0 / 3.0, 1.0 );
  cout << " gamma( 1/3 + i ) = " << std::setprecision( 8 )
       << Luna::gamma( x ) << endl;

  cout << " lngamma( 0.5 ) = " << std::setprecision( 6 )
       << Luna::lngamma( 0.5 ) << endl;

  cout << "--- FINISHED ---" << endl;

}
