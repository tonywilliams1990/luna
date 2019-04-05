/// \file Residual_with_coords.h
/// A specification of a  residual class that defines a vector residual of a
/// vector of state variables. The residual may also depend upon additional
/// coordinate variables.

#ifndef RESIDUAL_WITH_COORDS_H
#define RESIDUAL_WITH_COORDS_H

#include "Residual.h"
#include "Vector.h"
#include "Matrix.h"

namespace Luna
{
  // A doubly templated class 
  template <class T, class X = double >

  // A base class to define residuals with extra coordinate variables
  class Residual_with_coords : public Residual<T>
  {
  public:
    /// Constructor for a square residual object, N residuals for N unknowns.
    /// \param order The order of the residual vector
    /// \param ncoords The number of coordinates to store
    Residual_with_coords( const unsigned& order, const unsigned& ncoords );

    /// Constructor for a non-square residual object, there are less residual
    ///constraints than unknowns.
    /// \param order The number of residuals
    /// \param nvars The number of unknowns/variables
    /// \param ncoords The number of coordinates to store
    Residual_with_coords( const unsigned& order, const unsigned& nvars,
                          const unsigned& ncoords );

    /// An empty destructor
    virtual ~Residual_with_coords();

    /// General handle access to the coordinates
    /// \return A handle to the i-th coordinate
    X& coord( const unsigned& i );

    /// General handle access to the coordinates
    /// \return A handle to the i-th coordinate
    const X& coord( const unsigned& i ) const;

  protected:

    /// The coordinates stored for this residual
    std::vector<X> coords;

  };  // End of class Residual_with_coords

  template <typename T, typename X>
  Residual_with_coords<T, X>::Residual_with_coords( const unsigned& order,
                              const unsigned& ncoords ) : Residual<T>( order )
  {
    coords = std::vector<X>( ncoords, 0.0 );
  }

  template <typename T, typename X>
  Residual_with_coords<T, X>::Residual_with_coords( const unsigned& order,
                                                    const unsigned& nvars,
                        const unsigned& ncoords ) : Residual<T>( order, nvars )
  {
    coords = std::vector<X>( ncoords, 0.0 );
  }

  template <typename T, typename X>
  Residual_with_coords<T, X>::~Residual_with_coords()
  {
  }

  template <typename T, typename X>
  inline X& Residual_with_coords<T, X>::coord( const unsigned& i )
  {
    return coords[ i ];
  }

  template <typename T, typename X>
  inline const X& Residual_with_coords<T, X>::coord( const unsigned& i ) const
  {
    return coords[ i ];
  }

}  // End of namespace Luna

#endif
