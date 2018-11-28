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
#include <utility>

#include "Error.h"
#include "Vector.h"
#include "Mesh1D.h"
#include "Matrix.h"


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
    /// \param x_nodes Vector of nodal positions along the x-axis
    /// \param y_nodes Vector of nodal positions along the y-axis
    /// \param nvars The number of variables
    Mesh2D( const Vector<double>& x_nodes, const Vector<double>& y_nodes,
            const std::size_t& nvars ) : NX( x_nodes.size() ),
            NY( y_nodes.size() ), NV( nvars ), X_NODES( x_nodes ),
            Y_NODES( y_nodes )
    {
      // Create the storage Vector and set the contents to zero
      VARS = Vector<T>( NV * NX * NY, T( 0.0 ) );
    }

    /// Destructor
    ~Mesh2D() {}

    /* ----- Operator overloading ----- */

    /// Access operator for a nodal point/variable in the mesh
    /// \param nodex Node index ( x-axis )
    /// \param nodey Node index ( y-axis )
    /// \param var Variable index
    /// \return A reference to the variable stored at the given node
    T& operator()( const std::size_t nodex, const std::size_t nodey,
                   const std::size_t var );

    /// Const access operator for a nodal point/variable in the mesh
    /// \param nodex Node index ( x-axis )
    /// \param nodey Node index ( y-axis )
    /// \param var Variable index
    /// \return A constant reference to the variable stored at the given node
    const T& operator()( const std::size_t nodex, const std::size_t nodey,
                         const std::size_t var ) const;

    /// Copy assignment
    /// \param original 2D mesh to be copied
    /// \return The new 2D mesh
    Mesh2D<T>& operator=( const Mesh2D<T>& original );

    /* ----- Methods ----- */

    /// Return the spatial position of a given given node as a pair
    /// \param nodex Node index ( x-axis )
    /// \param nodey Node index ( y-axis )
    /// \return The nodal value pair at the specified index
    std::pair<double, double> coord( const std::size_t nodex,
                                     const std::size_t nodey ) const;

    /// Set the variables stored at a specified node
    /// \param nodex Node index ( x-axis )
    /// \param nodey Node index ( y-axis )
    /// \param vec Vector of variables
    void set_nodes_vars( const std::size_t nodex, const std::size_t nodey,
                         const Vector<T>& vec );

    /// Get the variables stored at a specified node
    /// \param nodex Node index ( x-axis )
    /// \param nodey Node index ( y-axis )
    /// \return The Vector of variables
    Vector<T> get_nodes_vars( const std::size_t nodex,
                              const std::size_t nodey ) const;

    /// Get a cross section of the 2D mesh at a specified (constant) x node
    /// \param nodex Node index ( x-axis )
    /// \return A 1D mesh cross section of the 2D mesh
    Mesh1D<T> get_cross_section_xnode( const std::size_t nodex ) const;

    /// Get a cross section of the 2D mesh at a specified (constant) y node
    /// \param nodey Node index ( y-axis )
    /// \return A 1D mesh cross section of the 2D mesh
    Mesh1D<T> get_cross_section_ynode( const std::size_t nodey ) const;

    /// Assign an element to all entries in the mesh
    /// \param elt Element to assign to all entries in the mesh
    void assign( const T elt );

    //TODO the rest of the 2D mesh methods


  };	// End of class Mesh2D

  template <typename T>
  inline T& Mesh2D<T>::operator()( const std::size_t nodex,
                                   const std::size_t nodey,
                                   const std::size_t var )
  {
    if ( nodex > NX - 1 || nodey > NY - 1 )
    {
      std::string problem;
      problem = " The Mesh2D.operator() method is trying to \n";
      problem += " access a nodal point that is not in the mesh. \n";
      throw Error( problem );
    }
    if ( var > NV - 1 )
    {
      std::string problem;
      problem = " The Mesh2D.operator() method is trying to \n";
      problem += " access a variable index that is not in the mesh. \n";
      throw Error( problem );
    }

    return VARS[ ( nodex * NY + nodey ) * NV + var ];
  }

  template <typename T>
  inline const T& Mesh2D<T>::operator()( const std::size_t nodex,
                                         const std::size_t nodey,
                                         const std::size_t var ) const
  {
    if ( nodex > NX - 1 || nodey > NY - 1 )
    {
      std::string problem;
      problem = " The Mesh2D.operator() method is trying to \n";
      problem += " access a nodal point that is not in the mesh. \n";
      throw Error( problem );
    }
    if ( var > NV - 1 )
    {
      std::string problem;
      problem = " The Mesh2D.operator() method is trying to \n";
      problem += " access a variable index that is not in the mesh. \n";
      throw Error( problem );
    }
    return VARS[ ( nodex * NY + nodey ) * NV + var ];
  }

  template <typename T>
  inline Mesh2D<T>& Mesh2D<T>::operator=( const Mesh2D<T>& original )
  {
    NX = original.NX;
    NY = original.NY;
    NV = original.NV;
    X_NODES = original.X_NODES;
    Y_NODES =  original.Y_NODES;
    VARS = original.VARS;
    return *this;
  }

  template <typename T>
  inline std::pair<double, double> Mesh2D<T>::coord( const std::size_t nodex,
                                              const std::size_t nodey ) const
  {
    if ( nodex > NX - 1 || nodey > NY - 1 || nodex < 0 || nodey < 0 )
    {
      std::string problem;
      problem = " The Mesh2D.coord method is trying to \n";
      problem += " access a nodal point that is not in the mesh. \n";
      throw Error( problem );
    }
    std::pair< double, double > pos;
    pos.first = X_NODES[ nodex ];
    pos.second = Y_NODES[ nodey ];
    return pos;
  }

  template <typename T>
  void Mesh2D<T>::set_nodes_vars( const std::size_t nodex,
                                  const std::size_t nodey,
                                  const Vector<T>& vec )
  {
    if ( nodex > NX - 1 || nodey > NY - 1 || nodex < 0 || nodey < 0 )
    {
      std::string problem;
      problem = " The Mesh2D.set_nodes_vars method is trying to \n";
      problem += " access a nodal point that is not in the mesh. \n";
      throw Error( problem );
    }
    if ( vec.size() > NV )
    {
      std::string problem;
      problem = " The Mesh2D.set_nodes_vars method is trying to use a \n";
      problem += " vector that has more entries than variables stored in the mesh. \n";
      throw Error( problem );
    }
    // assign contents of vec to the member data
    std::size_t offset( ( nodex * NY + nodey ) * NV );
    for ( std::size_t var = 0; var < NV; ++var )
    {
      VARS[ offset++ ] = vec[ var ];
    }
  }

  template <typename T>
  Vector<T> Mesh2D<T>::get_nodes_vars( const std::size_t nodex,
                                       const std::size_t nodey ) const
  {
    if ( nodex > NX - 1 || nodey > NY - 1 )
    {
      std::string problem;
      problem = " The Mesh2D.get_nodes_vars method is trying to \n";
      problem += " access a nodal point that is not in the mesh. \n";
      throw Error( problem );
    }
    //Construct a vector with NV elements
    Vector<T> nodes_vars;
    for ( std::size_t var = 0; var < NV; ++var )
    {
      nodes_vars.push_back( VARS[ ( nodex * NY + nodey ) * NV + var ] );
    }
    return nodes_vars;
  }

  template<typename T>
  Mesh1D<T> Mesh2D<T>::get_cross_section_xnode( const std::size_t nodex ) const
  {
    Mesh1D<T> section( Y_NODES, NV );
    for ( std::size_t nodey = 0; nodey < NY; ++nodey )
    {
      section.set_nodes_vars( nodey, this -> get_nodes_vars( nodex, nodey ) );
    }
    return section;
  }

  template<typename T>
  Mesh1D<T> Mesh2D<T>::get_cross_section_ynode( const std::size_t nodey ) const
  {
    Mesh1D<T> section( X_NODES, NV );
    for ( std::size_t nodex = 0; nodex < NX; ++nodex )
    {
      section.set_nodes_vars( nodex, this -> get_nodes_vars( nodex, nodey ) );
    }
    return section;
  }

  template <typename T>
  void Mesh2D<T>::assign( const T elt )
  {
    VARS.assign( NX * NY * NV, elt );
  }

}  // End of namespace Luna

#endif
