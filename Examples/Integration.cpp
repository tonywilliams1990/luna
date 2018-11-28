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

  

  cout << "--- FINISHED ---" << endl;

}
