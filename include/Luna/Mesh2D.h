/// \file Mesh2D.h
/// A class specifying a two dimensional mesh object used for storing and
/// manipulating data

#ifndef MESH2D_H
#define MESH2D_H

#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cassert>

#include "Error.h"
#include "Vector.h"


namespace Luna
{
  template <class T>

  class Mesh2D
  {
  private:
    std::size_t NX, NY, NV;       // Number of nodes and variables
    Vector<double> X_NODES;       // Vector for storing nodal points (x)
    Vector<double> Y_NODES;       // Vector for storing nodal points (y)
    Vector<T> VARS;               // Vector for storing variables at each node

  public:

    /// Constructor for an empty mesh
    Mesh2D() : NX( 0 ), NY( 0 ), NV( 0 )
    {}

    /// Constructor (given nodal distribution)
    Mesh2D( const Vector<double>& x_nodes, const Vector<double>& y_nodes,
            const std::size_t& nvars ) : NX( x_nodes.size() ),
            NY( y_nodes.size() ), NV( nvars ), X_NODES( x_nodes ), Y_NODES( y_nodes )
      {
        // Create the storage Vector and set the contents to zero
        VARS = Vector<T>( NV * NX * NY, T(0.0) );
      }


  };	// End of class Mesh2D



}  // End of namespace Luna

#endif
