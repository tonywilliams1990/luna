/// \file Mesh1D.h
/// A class specifying a one dimensional mesh object used for storing and
/// manipulating data

#ifndef MESH1D_H
#define MESH1D_H

#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <complex>

#include "Error.h"
#include "Vector.h"


namespace Luna
{

  /// A templated class
  template <class T, class X = double>

  class Mesh1D
  {
  private:
    std::size_t NV;                // Number of variables
    Vector<X> NODES;               // Vector for storing nodal points
    Vector<T> VARS;                // Vector for storing variables at each node

  public:

    /// Constructor for an empty mesh
    Mesh1D() : NV( 0 )
    {}

    /// Constructor (given nodal distribution)
    /// \param nodes The Vector of nodal points
    /// \param nvars The number of variables
    Mesh1D( const Vector<X>& nodes, const std::size_t& nvars ) :
        NV( nvars ), NODES( nodes )
    {
      // Create the storage Vector and set the contents to zero
      VARS = Vector<T>( NV * NODES.size(), T( 0.0 ) );
    }

    /// Destructor
    ~Mesh1D() {}

    /* ----- Operator overloading ----- */

    /// Read only access of a variable at a given node
    /// \param node Node index
    /// \param var Variable index
    /// \return A reference to the variable stored at the given node
    const T& operator()( const std::size_t node, const std::size_t var ) const
    {
      return VARS[ node * NV + var ];
    }

    /// Read and write access of a variable at a given node
    /// \param node Node index
    /// \param var Variable index
    /// \return A reference to the variable var stored at the index i
    T& operator()( const std::size_t node, const std::size_t var )
    {
      return VARS[ node * NV + var ];
    }

    /* ----- Methods ----- */

    /// Read only access to the nodal value
    /// \param node Node index
    /// \return The nodal value at the specified index
    const X& coord( const std::size_t& node ) const;

    /// Read and write access to the nodal value
    /// \param node Node index
    /// \return The nodal value at the specified index
    X& coord( const std::size_t& node );

    /// Set the variables stored at a specific node
    /// \param node Node index
    /// \param vec Vector of variables
    void set_nodes_vars( const std::size_t& node, const Vector<T>& vec );

    /// Get the variables stored at a specific node
    /// \param node Node index
    /// \return The Vector of variables
    Vector<T> get_nodes_vars( const std::size_t& node ) const;

    /// Return the number of nodal points in the mesh
    /// \return The number of nodal points
    std::size_t nnodes() const;

    /// Return the number of variables stored in the mesh
    /// \return The number of variables
    std::size_t nvars() const;

    /// Return the vector of nodal positions
    /// \return The Vector of nodal points
    const Vector<X>& nodes() const;

    /// Get the variables at an interpolated position ( first order scheme )
    /// \param pos Position within the nodal domain
    /// \return The Vector of interpolated variables
    Vector<T> get_interpolated_vars( const X& pos ) const;

    /// Output data to a file
    /// \param filename Name of the file to write to
    /// \param precision Decimal precision of the output
    void output( std::string filename, int precision = 15 ) const;

    /// Read data from a file
    /// \param filename Name of the file to read from
    /// \param precision Decimal precision of the input
    /// \param reset Reset the nodal positions after reading from file
    void read( std::string filename, int precision = 15, bool reset = false );

    /// Return a vector of all the variables stored in the mesh
    /// \return A reference to the Vector of the variables stored in the mesh
    const Vector<T>& vars_as_vector() const;

    /// Set the variables of this mesh from a vector
    /// \param A Vector of the variables to be stored in the mesh
    void set_vars_from_vector( const Vector<T>& vec );

    /// Integrate a given variable over the domain
    /// \param var Variable index
    T integral( std::size_t var = 0 ) const;


  };	// End of class Mesh1D

  template <typename T, typename X>
  inline const X& Mesh1D<T, X>::coord( const std::size_t& node ) const
  {
    return NODES[ node ];
  }

  template <typename T, typename X>
  inline X& Mesh1D<T, X>::coord( const std::size_t& node )
  {
    return NODES[ node ];
  }

  template <typename T, typename X>
  inline void Mesh1D<T, X>::set_nodes_vars( const std::size_t& node,
                                            const Vector<T>& vec )
  {
    if ( vec.size() != NV ) { throw Error( "Mesh error: set_nodes_vars " );}
    for ( std::size_t var = 0; var < vec.size(); ++var )
    {
      VARS[ node * NV + var ] = vec[ var ];
    }
  }

  template <typename T, typename X>
  inline Vector<T> Mesh1D<T, X>::get_nodes_vars( const std::size_t& node ) const
  {
    if ( ( node >= NODES.size() ) || ( node < 0 ) )
    {
      throw Error( "Mesh error: get_nodes_vars range error." );
    }
    Vector<T> nodes_vars;
    for ( std::size_t var = 0; var < NV; ++var )
    {
      nodes_vars.push_back( VARS[ node * NV + var ] );
    }
    return nodes_vars;
  }

  template <typename T, typename X>
  inline std::size_t Mesh1D<T, X>::nnodes() const
  {
    return NODES.size();
  }

  template <typename T, typename X>
  inline std::size_t Mesh1D<T, X>::nvars() const
  {
    return NV;
  }

  template <typename T, typename X>
  inline const Vector<X>& Mesh1D<T, X>::nodes() const
  {
    return NODES;
  }

  template <>
  inline Vector<double> Mesh1D<double, double>::get_interpolated_vars(
    const double& x_pos ) const
  {
    for ( unsigned node = 0; node < NODES.size() - 1; ++node )
    {
      // find bracketing nodes
      if ( ( ( NODES[ node ] < x_pos ) && ( NODES[ node + 1 ] > x_pos ) )
              ||  ( std::abs( NODES[ node ] - x_pos ) < 1.e-7  )
              || ( std::abs( NODES[ node + 1 ] - x_pos ) < 1.e-7 ) )
      {
        // distance from left node
        double delta_x( x_pos - NODES[ node ] );
        // empty data to return
        Vector<double> left;
        Vector<double> right;
        Vector<double> deriv;
        // interpolate data linearly
        left = get_nodes_vars( node );
        right = get_nodes_vars( node + 1 );
        deriv = ( right - left ) / ( NODES[ node + 1 ] - NODES[ node ] );
        // overwrite right
        right = left + deriv * delta_x;
        return right;
      }
    }
    // If the position is not in the range throw an error
    throw Error( "Mesh error: interpolation not in range " );
  }

  template <>
  inline Vector<std::complex<double> > Mesh1D<std::complex<double>,
                double>::get_interpolated_vars( const double& x_pos ) const
  {
    for ( unsigned node = 0; node < NODES.size() - 1; ++node )
    {
      // find bracketing nodes
      if ( ( NODES[ node ] < x_pos
        || std::abs( NODES[ node ] - x_pos ) < 1.e-7 ) &&
           ( NODES[ node + 1 ] > x_pos
        || std::abs( NODES[ node + 1 ] - x_pos ) < 1.e-7 ) )
      {
        // distance from left node
        double delta_x( x_pos - NODES[ node ] );
        // empty data to return
        Vector<std::complex<double> > left;
        Vector<std::complex<double> > right;
        Vector<std::complex<double> > deriv;
        // interpolate data linearly
        left = get_nodes_vars( node );
        right = get_nodes_vars( node + 1 );
        deriv = ( right - left ) / ( NODES[ node + 1 ] - NODES[ node ] );
        // overwrite right
        right = left + deriv * delta_x;
        return right;
      }
    }
    // If the position is not in the range throw an error
    throw Error( "Mesh error: interpolation not in range " );
  }

  template <>
  inline Vector<std::complex<double> > Mesh1D<std::complex<double>,
         std::complex<double> >::get_interpolated_vars(
         const std::complex<double>& pos ) const
  {
    double x_pos( pos.real() );
    for ( unsigned node = 0; node < NODES.size() - 1; ++node )
    {
      // find bracketing nodes
      if ( ( NODES[ node ].real() < x_pos
        || std::abs( NODES[ node ].real() - x_pos ) < 1.e-7 ) &&
           ( NODES[ node + 1 ].real() > x_pos
        || std::abs( NODES[ node + 1 ].real() - x_pos ) < 1.e-7 ) )
      {
        // distance from left node -- real coordinate is given. We also need to
        // interpolate between the two complex nodes
        std::complex<double> delta_z = ( NODES[ node + 1 ] - NODES[ node ] )
                                       * ( x_pos - NODES[ node ].real() )
                                       / ( NODES[ node + 1 ].real()
                                         - NODES[ node ].real() );
        // empty data to return
        Vector<std::complex<double> > left;
        Vector<std::complex<double> > right;
        Vector<std::complex<double> > deriv;
        // interpolate data linearly
        left = get_nodes_vars( node );
        right = get_nodes_vars( node + 1 );
        // derivative of the data
        deriv = ( right - left ) / ( NODES[ node + 1 ] - NODES[ node ] );
        // overwrite right
        right = left + deriv * delta_z;
        return right;
      }
    }
    // If the position is not in the range throw an error
    throw Error( "Mesh error: interpolation not in range " );
  }

  template <typename T, typename X>
  inline void Mesh1D<T, X>::output( std::string filename, int precision ) const
  {
    std::ofstream out;
    out.open( filename.c_str() );
    out.precision( precision );
    out.setf( std::ios::showpoint );
    out.setf( std::ios::showpos );
    out.setf( std::ios::scientific );
    for ( std::size_t i = 0; i < NODES.size(); ++i )
    {
      out << NODES[ i ] << " ";
      for ( std::size_t var = 0; var < NV; ++var )
      {
        out << VARS[ i * NV + var ] << " ";
      }
      out << "\n";
    }
  }

  template <>
	inline void Mesh1D<double, double>::read( std::string filename,
                                            int precision, bool reset)
	{
		std::ifstream inFile;
		inFile.open( filename.c_str() );
		inFile.precision( precision );
		inFile.setf( std::ios::showpoint );
		inFile.setf( std::ios::showpos );
		inFile.setf( std::ios::scientific );
		if (!inFile) {
			throw Error( "Mesh1D unable to open file" );
		}

		double x;
		Vector<double> vals( NV );

		for ( std::size_t i=0; i<NODES.size(); ++i )
		{
			inFile >> x;
			for ( std::size_t v=0; v<vals.size(); ++v )
			{
				inFile >> vals[ v ];
				VARS[ i * NV + v ] = vals[ v ];
			}
			// Check the mesh
			if ( !reset )
			{
				if ( ( std::abs( x - NODES[ i ] ) > 1.e-6 ) )
				{
					std::cout << " Read x = " << x
                    << " Expected x = " << NODES[ i ] << std::endl;
					throw Error( "Mesh1D: file has different nodal points." );
				}
			}
			else
			{
				NODES[ i ] = x;
			}
		}
		inFile.close();
	}

  //TODO complex version of reading from a file

  template <typename T, typename X>
  inline const Vector<T>& Mesh1D<T, X>::vars_as_vector() const
  {
    return VARS;
  }

  template <typename T, typename X>
  inline void Mesh1D<T, X>::set_vars_from_vector( const Vector<T>& vec )
  {
    if ( vec.size() != NV * NODES.size() )
    { throw Error( "Mesh1D: set_vars_from_vector sizes do not agree." ); }
    VARS = vec;
  }

  template <typename T, typename X>
  inline T Mesh1D<T, X>::integral( std::size_t var ) const
  {
    T sum( 0.0 );
    X dx( 0.0 );

    for ( std::size_t node = 0; node < NODES.size() - 1; ++node )
    {
      dx = ( NODES[ node + 1 ] - NODES[ node ] );
      sum += 0.5 * dx * ( VARS[ node * NV + var ]
                        + VARS[ ( node + 1 ) * NV + var ] );
    }
    return sum;
  }

}  // End of namespace Luna

#endif
