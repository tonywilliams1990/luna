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
#include <complex>

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

    /// Get the number of nodes in the two directions of the 2D mesh
    /// \return A std::pair containing the number of nodes in each direction
    std::pair< std::size_t, std::size_t > get_nnodes() const;

    /// Get the number of variables that are stored at each node
    /// \return The number of variables stored in the mesh
    std::size_t get_nvars() const;

    /// Return A vector of the x-nodal positions for this mesh
    /// \return A the mesh x nodes Vector
    const Vector<double>& xnodes() const;

    /// Return A vector of the y-nodal positions for this mesh
    /// \return A the mesh y nodes Vector
    const Vector<double>& ynodes() const;

    /// Return a Matrix corresponding to each nodal point in the mesh
    /// \param var Variable index
    /// \return A Matrix for a specified variable at each nodal point
    Matrix<T> get_var_as_matrix( std::size_t var ) const;

    /// Interpolate the mesh data (bilinearly) into a new mesh
    /// \param newX New x nodal points Vector
    /// \param newY New y nodal points Vector
    void remesh( const Vector<double>& newX, const Vector<double>& newY );

    /// Output mesh data to std::cout
    void output() const;

    /// Output mesh data to a file
    /// \param filename Name of the file to output the mesh data to
    void output( std::string filename ) const;

    /// Output a single variable to a file with no nodal information
    /// \param filename Name of the file to output the data to
    /// \param var Variable index
    void output_var( std::string filename, const std::size_t var ) const;

    /// A simple method for reading data from a file
    /// \param filename Name of the file to read the data from
    /// \param reset Reset the mesh nodes using data from file
    void read( std::string filename, const bool reset = false );

    /// Output data to a file for gnuplot surface plotting
    /// \param filename Name of the file to output the mesh data to
    void output_gnu( std::string filename ) const;

    /// Get a bilinearly interpolated value at a specified point
    /// \param x The x-coordinate
    /// \param y The y-coordinate
    /// \return The interpolated Vector of variables at a specified point
    Vector<T> get_interpolated_vars( const double& x, const double& y );

    /// Numerically integrate a given variable over the domain
    /// \param var Variable index
    /// \return The integral over the domain
    T integral2D( std::size_t var = 0 ) const;

    /// Numerically integrate the square of the absolute value over the domain
    /// \param var Variable index
    /// \return The integral of the square of the absolute value over the domain
    T square_integral2D( std::size_t var = 0 ) const;

    /// Return the Vector of variables
    /// \return A copy of the Vector of variables
    Vector<T> get_vars();

    /// Fill the mesh points with values determined by a given function
    /// \param func The 2 variable function used to assign values
    /// \param var The variable to assign the values to
    void apply( T func ( const T&, const T& ), const int& var );

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
  inline void Mesh2D<T>::set_nodes_vars( const std::size_t nodex,
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
      problem = " The Mesh2D.set_nodes_vars method is using a vector that \n";
      problem += " has more entries than variables stored in the mesh. \n";
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

  template <typename T>
  std::pair< std::size_t, std::size_t > Mesh2D<T>::get_nnodes() const
  {
    std::pair< std::size_t, std::size_t > nodes;
    nodes.first = NX;
    nodes.second = NY;
    return nodes;
  }

  template <typename T>
  std::size_t Mesh2D<T>::get_nvars() const
  {
    return NV;
  }

  template <typename T>
  const Vector<double>& Mesh2D<T>::xnodes() const
  {
    return X_NODES;
  }

  template <typename T>
  const Vector<double>& Mesh2D<T>::ynodes() const
  {
    return Y_NODES;
  }

  template <typename T>
  Matrix<T> Mesh2D<T>::get_var_as_matrix( std::size_t var ) const
  {
    if ( var > NV - 1 )
    {
      std::string problem;
      problem = " The Mesh2D.get_var_as_matrix method is trying to use an \n";
      problem += " index bigger than the number of variables in the mesh. \n";
      throw Error( problem );
    }
    Matrix<T> temp( NX, NY, 0.0 );
    for ( std::size_t i = 0; i < NX; ++i )
    {
      for  ( std::size_t j = 0; j < NY; ++j )
      {
        temp( i, j ) = VARS[ ( i * NY + j ) * NV + var ];
      }
    }
    return temp;
  }

  template<class T>
  void Mesh2D<T>::remesh( const Vector<double>& newX,
                          const Vector<double>& newY )
  {
    // check start & end
    if ( std::abs( X_NODES[ 0 ] - newX[ 0 ] ) > 1.e-10 ||
         std::abs( X_NODES[ X_NODES.size() - 1 ] - newX[ newX.size() - 1 ] )
         > 1.e-10 )
    {
      std::string problem;
      problem = " The Mesh2D.remesh1 method has been called with a passed \n";
      problem += " X coordinate vector that has different start and/or \n";
      problem += " end points from the instantiated object. \n";
      throw Error( problem );
    }
    // check monotonic node positions
    for ( std::size_t i = 0; i < newX.size() - 1; ++i )
    {
      if ( newX[ i ] >= newX[ i + 1 ] )
      {
        std::string problem;
        problem = " The Mesh2D.remesh1 method has been passed \n";
        problem += " a non-monotonic X coordinate vector. \n";
        throw Error( problem );
      }
    }
    // check start and end
    if ( std::abs( Y_NODES[ 0 ] - newY[ 0 ] ) > 1.e-10 ||
         std::abs( Y_NODES[ Y_NODES.size() - 1 ] - newY[ newY.size() - 1 ] )
         > 1.e-10 )
    {
      std::string problem;
      problem = " The Mesh2D.remesh1 method has been called with a passed \n";
      problem += " Y coordinate vector that has different start and/or \n";
      problem += " end points from the instantiated object. \n";
      throw Error( problem );
    }
    // check monotonic node positions
    for ( std::size_t i = 0; i < newY.size() - 1; ++i )
    {
      if ( newY[ i ] >= newY[ i + 1 ] )
      {
        std::string problem;
        problem = " The Mesh2D.remesh1 method has been passed \n";
        problem += " a non-monotonic Y coordinate vector. \n";
        throw Error( problem );
      }
    }

    // new variables storage
    Vector<T> newvars( newX.size() * newY.size() * NV, 0.0 );

    // left boundary
    {
      std::size_t xnode( 0 );
      // bottom left corner copy
      for ( unsigned var = 0; var < NV; ++var )
      {
        newvars[ ( xnode * newY.size() + 0 ) * NV + var ] =
                                                  get_nodes_vars( 0, 0 )[ var ];
      }
      for ( std::size_t ynode = 1; ynode < newY.size() - 1; ++ynode )
      {
        std::size_t left_i( 0 );  // bracketing index
        std::size_t below_j( 0 ); // bracketing index
        double deltaY( 0.0 );
        // loop through the source mesh and find the bracket-nodes
        for ( std::size_t j = 0; j < Y_NODES.size() - 1; ++j )
        {
          if ( ( Y_NODES[ j ] <= newY[ ynode ] ) &&
               ( newY[ ynode ] < Y_NODES[ j + 1 ] ) )
          {
            below_j = j;
            deltaY = newY[ ynode ] - Y_NODES[ j ];
          }
        }
        Vector<T> dvarsdY = ( get_nodes_vars( left_i, below_j + 1 )
                              - get_nodes_vars( left_i, below_j ) )
                           / ( coord( left_i, below_j + 1 ).second
                              - coord( left_i, below_j ).second );
        Vector<T> interpolated_vars =   get_nodes_vars( left_i, below_j )
                                      + dvarsdY * deltaY;
        for ( std::size_t var = 0; var < NV; ++var )
        {
          newvars[ ( xnode * newY.size() + ynode ) * NV + var ] =
                                                       interpolated_vars[ var ];
        }
      }
      // top left corner copy
      for ( std::size_t var = 0; var < NV; ++var )
      {
        newvars[ ( xnode * newY.size() + newY.size() - 1 ) * NV + var ] =
        get_nodes_vars( 0, NY - 1 )[ var ];
      }
    }
    // right boundary
    {
      std::size_t xnode( newX.size() - 1 );
      // bottom right corner copy
      for ( std::size_t var = 0; var < NV; ++var )
      {
        newvars[ ( xnode * newY.size() + 0 ) * NV + var ] =
        get_nodes_vars( NX - 1, 0 )[ var ];
      }
      for ( std::size_t ynode = 1; ynode < newY.size() - 1; ++ynode )
      {
        std::size_t left_i( X_NODES.size() - 1 );  // bracketing index
        std::size_t below_j( 0 ); // bracketing index
        double deltaY( 0.0 );
        // loop through the source mesh and find the bracket-nodes
        for ( std::size_t j = 0; j < Y_NODES.size() - 1; ++j )
        {
          if ( ( Y_NODES[ j ] <= newY[ ynode ] ) &&
               ( newY[ ynode ] < Y_NODES[ j + 1 ] ) )
          {
            below_j = j;
            deltaY = newY[ ynode ] - Y_NODES[ j ];
          }
        }
        Vector<T> dvarsdY = ( get_nodes_vars( left_i, below_j + 1 )
                              - get_nodes_vars( left_i, below_j ) )
                           / ( coord( left_i, below_j + 1 ).second
                              - coord( left_i, below_j ).second );
        Vector<T> interpolated_vars =   get_nodes_vars( left_i, below_j )
                                      + dvarsdY * deltaY;
        for ( std::size_t var = 0; var < NV; ++var )
        {
          newvars[ ( xnode * newY.size() + ynode ) * NV + var ] =
                                                       interpolated_vars[ var ];
        }
      }
      // bottom right corner copy
      for ( std::size_t var = 0; var < NV; ++var )
      {
        newvars[ ( xnode * newY.size() + newY.size() - 1 ) * NV + var ] =
        get_nodes_vars( NX - 1, NY - 1 )[ var ];
      }
    }
    // bottom boundary
    {
      std::size_t ynode( 0 );
      for ( std::size_t xnode = 1; xnode < newX.size() - 1; ++xnode )
      {
        std::size_t left_i( 0 );  // bracketing index
        std::size_t below_j( 0 ); // bracketing index
        double deltaX( 0.0 );
        // loop through the source mesh and find the bracket-nodes
        for ( std::size_t i = 0; i < X_NODES.size() - 1; ++i )
        {
          if ( ( X_NODES[ i ] <= newX[ xnode ] ) &&
               ( newX[ xnode ] < X_NODES[ i + 1 ] ) )
          {
            left_i = i;
            deltaX = newX[ xnode ] - X_NODES[ i ];
          }
        }
        Vector<T> dvarsdX = ( get_nodes_vars( left_i + 1, below_j )
                              - get_nodes_vars( left_i, below_j ) )
                           / ( coord( left_i + 1, below_j ).first
                              - coord( left_i, below_j ).first );
        Vector<T> interpolated_vars =   get_nodes_vars( left_i, below_j )
                                      + dvarsdX * deltaX;
        for ( std::size_t var = 0; var < NV; ++var )
        {
          newvars[ ( xnode * newY.size() + ynode ) * NV + var ] =
                                                       interpolated_vars[ var ];
        }
      }
    }
    // top boundary
    {
      std::size_t ynode( newY.size() - 1 );
      for ( std::size_t xnode = 1; xnode < newX.size() - 1; ++xnode )
      {
        std::size_t left_i( 0 );  // bracketing index
        std::size_t below_j( Y_NODES.size() - 1 ); // bracketing index
        double deltaX( 0.0 );
        // loop through the source mesh and find the bracket-nodes
        for ( std::size_t i = 0; i < X_NODES.size() - 1; ++i )
        {
          if ( ( X_NODES[ i ] <= newX[ xnode ] ) &&
               ( newX[ xnode ] < X_NODES[ i + 1 ] ) )
          {
            left_i = i;
            deltaX = newX[ xnode ] - X_NODES[ i ];
          }
        }
        Vector<T> dvarsdX = ( get_nodes_vars( left_i + 1, below_j )
                              - get_nodes_vars( left_i, below_j ) )
                           / ( coord( left_i + 1, below_j ).first
                              - coord( left_i, below_j ).first );
        Vector<T> interpolated_vars =   get_nodes_vars( left_i, below_j )
                                      + dvarsdX * deltaX;
        for ( std::size_t var = 0; var < NV; ++var )
        {
          newvars[ ( xnode * newY.size() + ynode ) * NV + var ] =
                                                       interpolated_vars[ var ];
        }
      }
    }
    // loop thru interior nodes of the destination mesh one node at a time
    for ( std::size_t xnode = 1; xnode < newX.size() - 1; ++xnode )
    {
      for ( std::size_t ynode = 1; ynode < newY.size() - 1; ++ynode )
      {
        std::size_t left_i( 0 );  // bracketing index
        std::size_t below_j( 0 ); // bracketing index
        // loop through the source mesh and find the bracket-nodes
        for ( std::size_t i = 0; i < X_NODES.size() - 1; ++i )
        {
          if ( ( X_NODES[ i ] <= newX[ xnode ] ) &&
               ( newX[ xnode ] < X_NODES[ i + 1 ] ) )
          {
            left_i = i;
          }
        }
        // loop through the source mesh and find the bracket-nodes
        for ( std::size_t j = 0; j < Y_NODES.size() - 1; ++j )
        {
          if ( ( Y_NODES[ j ] <= newY[ ynode ] ) &&
               ( newY[ ynode ] < Y_NODES[ j + 1 ] ) )
          {
            below_j = j;
          }
        }
        Vector<T> dvarsdX = ( get_nodes_vars( left_i + 1, below_j )
                              - get_nodes_vars( left_i, below_j ) )
                           / ( coord( left_i + 1, below_j ).first
                              - coord( left_i, below_j ).first );
        Vector<T> dvarsdY = ( get_nodes_vars( left_i, below_j + 1 )
                              - get_nodes_vars( left_i, below_j ) )
                           / ( coord( left_i, below_j + 1 ).second
                              - coord( left_i, below_j ).second );

        Vector<T> interpolated_vars_bottom =
                    ( get_nodes_vars( left_i, below_j ) *
                      ( coord( left_i + 1, below_j ).first - newX[ xnode ] )
                      + get_nodes_vars( left_i + 1, below_j )
                      * ( newX[ xnode ] - coord( left_i, below_j ).first ) )
                    / ( coord( left_i + 1, below_j ).first
                      - coord( left_i, below_j ).first );

        Vector<T> interpolated_vars_top =
          ( get_nodes_vars( left_i, below_j + 1 )
          * ( coord( left_i + 1, below_j + 1 ).first - newX[ xnode ] )
          + get_nodes_vars( left_i + 1, below_j + 1 )
          * ( newX[ xnode ] - coord( left_i, below_j + 1 ).first ) ) /
          ( coord( left_i + 1, below_j + 1 ).first
          - coord( left_i, below_j + 1 ).first );

        Vector<T> interpolated_vars =
          (  interpolated_vars_bottom * ( coord( left_i, below_j + 1 ).second
          - newY[ ynode ] ) +  interpolated_vars_top
          * ( newY[ ynode ] - coord( left_i, below_j ).second ) ) /
          ( coord( left_i, below_j + 1 ).second
          - coord( left_i, below_j ).second );

        for ( std::size_t var = 0; var < NV; ++var )
        {
          newvars[ ( xnode * newY.size() + ynode ) * NV + var ] =
                                                       interpolated_vars[ var ];
        }
      }
    }
    // finally replace the old nodes with the new ones
    X_NODES = newX;
    Y_NODES = newY;
    NX = newX.size();
    NY = newY.size();
    VARS = newvars;
  }

  template <typename T>
  void Mesh2D<T>::output() const
  {
    for ( std::size_t var = 0; var < NV; ++var )
    {
      std::cout << "Variable : " << var << "\n";
      std::cout << " x = ";
      for ( std::size_t i = 0; i < NX; ++i )
      {
        std::cout << X_NODES[ i ] << ", ";
      }
      std::cout << "\n";
      for ( std::size_t j = 0; j < NY; ++j )
      {
        std::cout << " y = " << Y_NODES[ j ] << "\n";
        for ( std::size_t i = 0; i < NX; ++i )
        {
          std::cout << VARS[ ( i * NY + j ) * NV + var ] << ", ";
        }
        std::cout << "\n";
      }
    }
  }

  template <typename T>
  void Mesh2D<T>::output( std::string filename ) const
  {
    std::ofstream out;
    out.open( filename.c_str() );
    out.precision( 12 );
    out.setf( std::ios::showpoint );
    out.setf( std::ios::showpos );
    out.setf( std::ios::scientific );
    out.precision( 12 );
    for ( std::size_t j = 0; j < NY; ++j )
    {
      for ( std::size_t i = 0; i < NX; ++i )
      {
        out << X_NODES[ i ] << " " << Y_NODES[ j ] << " ";
        for ( std::size_t var = 0; var < NV; ++var )
        {
          out << VARS[ ( i * NY + j ) * NV + var ] << " ";
        }
        out << "\n";
      }
      out << "\n";
    }
  }

  template<typename T>
  void Mesh2D<T>::output_var( std::string filename,
                              const std::size_t var ) const
  {
    std::ofstream out;
    out.open( filename.c_str() );
    out.precision( 12 );
    out.setf( std::ios::showpoint );
    out.setf( std::ios::showpos );
    out.setf( std::ios::scientific );
    out.precision( 12 );
    for ( std::size_t j = 0; j < NY; ++j )
    {
      for ( std::size_t i = 0; i < NX; ++i )
      {
        out << VARS[ ( i * NY + j ) * NV + var ] << "\n";
      }
    }
  }

  template<>
  void Mesh2D<double>::read( std::string filename, const bool reset )
  {
    std::ifstream inFile;
    inFile.open( filename.c_str() );
    inFile.precision( 12 );
    inFile.setf( std::ios::showpoint );
    inFile.setf( std::ios::showpos );
    inFile.setf( std::ios::scientific );
    if (!inFile) {
      throw Error( "Mesh2D unable to open file" );
    }

    double x, y;
    Vector<double> vals( NV );

    for ( std::size_t i = 0; i < NX; ++i )
    {
      for ( std::size_t j = 0; j < NY; ++j )
      {
        inFile >> x;
        inFile >> y;
        std::size_t offset( ( i * NY + j ) * NV );
        for ( std::size_t v = 0; v < NV; ++v )
        {
          inFile >> vals[ v ];
          VARS[ offset++ ] = vals[ v ];
        }
        // Check the mesh
        if ( !reset )
  			{
  				if ( ( std::abs( x - X_NODES[ i ] ) > 1.e-6 ) ||
               ( std::abs( y - Y_NODES[ j ] ) > 1.e-6 ) )
  				{
  					std::cout << " Read x = " << x << " Expected x = "
                      << X_NODES[ i ] << std::endl;
            std::cout << " Read y = " << y << " Expected y = "
                      << Y_NODES[ i ] << std::endl;
            std::string problem;
            problem  = " Mesh2D trying to read from file with different nodal";
            problem += " points.";
  					throw Error( problem );
  				}
  			}
  			else
  			{
  				X_NODES[ i ] = x;
          Y_NODES[ j ] = y;
  			}
      }
    }
    inFile.close();
  }

  template<>
  void Mesh2D<double>::output_gnu( std::string filename ) const
  {
    std::ofstream out;
    out.open( filename.c_str() );
    out.precision( 12 );
    out.setf( std::ios::showpoint );
    out.setf( std::ios::showpos );
    out.setf( std::ios::scientific );

    for ( std::size_t i = 0; i < NX; ++i )
    {
      for ( std::size_t j = 0; j < NY; ++j )
      {
        out << X_NODES[ i ] << " " << Y_NODES[ j ] << " ";
        for ( std::size_t var = 0; var < NV; ++var )
        {
          out << VARS[ ( i * NY + j ) * NV + var ] << " ";
        }
        out << "\n";
      }
      out << "\n";
    }
    out.close();
  }

  template <>
  void Mesh2D< std::complex<double> >::output_gnu( std::string filename ) const
  {
    std::ofstream out;
    out.open( filename.c_str() );
    out.precision( 12 );
    out.setf( std::ios::showpoint );
    out.setf( std::ios::showpos );
    out.setf( std::ios::scientific );
    //
    for ( std::size_t i = 0; i < NX; ++i )
    {
      for ( std::size_t j = 0; j < NY; ++j )
      {
        out << X_NODES[ i ] << " " << Y_NODES[ j ] << " ";
        for ( std::size_t var = 0; var < NV; ++var )
        {
          out << real(VARS[ ( i * NY + j ) * NV + var ]) << " ";
          out << imag(VARS[ ( i * NY + j ) * NV + var ]) << " ";
        }
        out << "\n";
      }
      out << "\n";
    }
    out.close();
  }

  template <typename T>
  Vector<T> Mesh2D<T>::get_interpolated_vars( const double& x, const double& y)
  {
    if ( ( x < X_NODES[ 0 ] ) || ( x > X_NODES[ NX - 1 ] ) )
    {
      std::string problem;
      problem = "The Mesh2D.get_interpolated_vars method has been called \n";
      problem += "with an x coordinate that lies outside the mesh. \n";
      throw Error( problem );
    }
    if ( ( y < Y_NODES[ 0 ] ) || ( y > Y_NODES[ NY - 1 ] ) )
    {
      std::string problem;
      problem = "The Mesh2D.get_interpolated_vars method has been called \n";
      problem += "with a y coordinate that lies outside the mesh. \n";
      throw Error( problem );
    }
    int bottom_j( -2 );
    for ( unsigned j = 0; j < NY - 1; ++j )
    {
      if ( ( y > Y_NODES[ j ] ) && ( y < Y_NODES[ j + 1 ] ) )
      {
        bottom_j = j;
      }
      if ( ( abs( y - Y_NODES[ j ] ) < 1.e-10 ) ||
           ( abs( y - Y_NODES[ j + 1 ] ) < 1.e-10 ) )
      {
        bottom_j = j;
      }
    }
    if ( bottom_j == -1 )
    {
      std::string problem;
      problem = " The Mesh2D.get_interpolated_vars method is broken.\n";
      throw Error( problem );
    }

    //std::cout << y << " " << Y_NODES[bottom_j] << " "
    //          << Y_NODES[bottom_j+1] << "\n";

    Mesh1D<T> bottom_row = get_cross_section_ynode( bottom_j );
    Mesh1D<T> top_row = get_cross_section_ynode( bottom_j + 1 );
    const double y1 = Y_NODES[ bottom_j ];
    const double y2 = Y_NODES[ bottom_j + 1 ];
    Vector<T> result = top_row.get_interpolated_vars( x )*( y2-y )/( y2-y1 )
      + bottom_row.get_interpolated_vars(x)*( y-y1 )/( y2-y1 );
    //std::cout << "x,y,interp: " << x << " " << y << " " << result[0] << "\n";
    return result;
  }

  template<typename T>
  T Mesh2D<T>::integral2D( std::size_t var ) const
  {
    T sum( 0.0 );
    double dx( 0.0 );
    double dy( 0.0 );
    // Sum over the elements
    for ( std::size_t i=0; i < NX - 1; ++i )
    {
      dx = X_NODES[ i + 1 ] - X_NODES[ i ];
      for ( std::size_t j=0; j < NY - 1; ++j )
      {
        dy = Y_NODES[ j + 1 ] - Y_NODES[ j ];
        sum += 0.25 * dx * dy * ( VARS[ ( i * NY + j ) * NV + var ]
            + VARS[ ( ( i + 1 ) * NY + j ) * NV + var ]
            + VARS[ ( i * NY + j + 1 ) * NV + var ]
            + VARS[ ( ( i + 1 ) * NY + j + 1 ) * NV + var ] );
      }
    }
    return sum;
  }

  template <typename T>
  T Mesh2D<T>::square_integral2D( std::size_t var ) const
  {
    T sum( 0.0 );
    double dx( 0.0 );
    double dy( 0.0 );
    // Sum over the elements
    for ( std::size_t i=0; i < NX - 1; ++i )
    {
      dx = X_NODES[ i + 1 ] - X_NODES[ i ];
      for ( std::size_t j=0; j < NY - 1; ++j )
      {
        dy = Y_NODES[ j + 1 ] - Y_NODES[ j ];
        sum += 0.25 * dx * dy * (
            std::pow( std::abs( VARS[ ( i * NY + j ) * NV + var ] ), 2 )
          + std::pow( std::abs( VARS[ ( ( i + 1 ) * NY + j ) * NV + var ] ), 2 )
          + std::pow( std::abs( VARS[ ( i * NY + j + 1 ) * NV + var ] ), 2 )
          + std::pow( std::abs( VARS[ ( ( i + 1 ) * NY + j + 1 ) * NV + var ] )
                                                                        , 2 ) );
      }
    }
    return sum;
  }

  template <typename T>
  inline Vector<T> Mesh2D<T>::get_vars()
  {
    Vector<T> temp;
    temp = VARS;
    return temp;
  }

  template <typename T>
  inline void Mesh2D<T>::apply( T func ( const T&, const T& ), const int& var )
  {
    for ( std::size_t i = 0; i < NX; i++ )
    {
      double x( X_NODES[ i ] );
      for ( std::size_t j = 0; j < NY; j++ )
      {
        double y( Y_NODES[ j ] );
        VARS[ ( i * NY + j ) * NV + var ] = func( x, y );
      }
    }
  }

}  // End of namespace Luna

#endif
