/// \file  Integration.cpp
/// \ingroup Examples
/// Some numerical integration examples

#include "Luna/Core"

using namespace std;
using namespace Luna;

int main()
{
  cout << "----------------- Integration -----------------" << endl;

  Vector<double> nodes;
  nodes.linspace( 0, 1, 101 );
  std::size_t nvars( 2 );
  Mesh1D<double> mesh( nodes, nvars );

  for ( std::size_t i = 0; i< nodes.size(); ++i )
  {
    // Set the first variable equal to 2 * x
    mesh( i, 0 ) = 2 * mesh.coord( i );
    // Set the second variable equal to x^2
    mesh( i, 1 ) = mesh.coord( i ) * mesh.coord( i );
  }

  cout << "  * mesh.nnodes() = " << mesh.nnodes() << endl;
  cout << "  * mesh.nvars() = " << mesh.nvars() << endl;

  Vector<double> vars;
  vars = mesh.get_interpolated_vars( 0.35 );
  cout << "  * vars( 0.35 ) = " << vars << endl;

  //mesh.output( "./test_output.txt" );

  cout << "  * mesh.integral( 0 ) = " << mesh.integral( 0 ) << endl;
  cout << "  * mesh.integral( 1 ) = " << mesh.integral( 1 ) << endl;

  cout << "-----------------------------------------------" << endl;

  Vector<double> x_nodes, y_nodes;
  x_nodes.linspace( 0, 1, 11 );
  y_nodes.linspace( 0, 1, 11 );
  Mesh2D<double> mesh2d( x_nodes, y_nodes, nvars );

  for ( std::size_t i = 0; i < x_nodes.size(); ++i )
  {
    for ( std::size_t j = 0; j < y_nodes.size(); ++j )
    {
      double x( mesh2d.coord( i, j ).first );
      double y( mesh2d.coord( i, j ).second );
      // Set the first variable equal to 2 * x^2
      mesh2d( i, j, 0 ) = 2 * x * x;
      // Set the second variable equal to 3 * y
      mesh2d( i, j, 1 ) = 3 * y;
    }
  }

  cout << " mesh2d( 1, 2, 0 ) = " << mesh2d( 1, 2, 0 ) << endl;
  cout << " mesh2d( 1, 2, 1 ) = " << mesh2d( 1, 2, 1 ) << endl;

  Vector<double> vec( 2, 1.1 );
  mesh2d.set_nodes_vars( 1, 2, vec );
  cout << " mesh2d( 1, 2, 0 ) = " << mesh2d( 1, 2, 0 ) << endl;
  cout << " mesh2d( 1, 2, 1 ) = " << mesh2d( 1, 2, 1 ) << endl;

  //mesh2d.assign( 3.14 );
  vec =  mesh2d.get_nodes_vars( 3, 3 );
  cout << " vec = " << vec << endl;

  cout << " # x nodes = " << mesh2d.get_nnodes().first << endl;
  cout << " # y nodes = " << mesh2d.get_nnodes().second << endl;
  cout << " # vars = " << mesh2d.get_nvars() << endl;
  cout << " mesh2d.xnodes() = " << mesh2d.xnodes() << endl;
  cout << " mesh2d.ynodes() = " << mesh2d.ynodes() << endl;

  cout << " mesh2d.get_var_as_matrix( 0 ) = "
       << mesh2d.get_var_as_matrix( 0 ) << endl;

  Vector<double> new_x_nodes, new_y_nodes;
  new_x_nodes.linspace( 0, 1, 6 );
  new_y_nodes.linspace( 0, 1, 6 );
  mesh2d.remesh( new_x_nodes, new_y_nodes );

  cout << " # x nodes = " << mesh2d.get_nnodes().first << endl;
  cout << " # y nodes = " << mesh2d.get_nnodes().second << endl;

  cout << " mesh2d.get_var_as_matrix( 0 ) = "
       << mesh2d.get_var_as_matrix( 0 ) << endl;

  //mesh2d.output_gnu( "./test.dat" );

  cout << " mesh2d.get_interpolated_vars( 0.5, 0.5 ) = "
       << mesh2d.get_interpolated_vars( 0.5, 0.5 ) << endl;

  cout << " mesh2d.integral2D( 0 ) = " << mesh2d.integral2D( 0 ) << endl;
  cout << " mesh2d.square_integral2D( 0 ) = "
       << mesh2d.square_integral2D( 0 ) << endl;


  cout << "--- FINISHED ---" << endl;

}
