/// \file Vector.h
/// A templated numeric Vector class that encapsulates std::vector.

#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <functional>

#include "Error.h"

namespace Luna
{
  // Forward declare the Matrix class that will be friends with Vector
  template <class T>
  class Matrix;

  // Forward declare the SparseMatrix class that will be friends with Vector
  template <class T>
  class SparseMatrix;

	// Templated vector class
	template <class T>

  /// A Vector class for use with double and std::complex<double>
  class Vector
  {
    private:
    std::vector<T> VECTOR;            // Data storage encapsulates std::vector

    protected:
    friend class Matrix<T>;			      // Friend class can access VECTOR directly
    friend class SparseMatrix<T>;	    // Friend class can access VECTOR directly

    public:

    /// Constructor for an empty Vector of unspecified size
    Vector(){}

    /// Constructor for a Vector of specified size
    /// \param size Size of the Vector
    Vector( const std::size_t& size );

    /// Constructor of specified size and elements
    /// \param size Size of the Vector
    /// \param elem Element to fill each entry of the Vector with
    Vector( const std::size_t& size, const T& elem );

    /// Copy constructor
    /// \param source Vector to be copied when initialising
    Vector( const Vector<T>& source );

    /// Destructor
    ~Vector(){}

    /* ----- Operator overloading ----- */

    /// Output operator <<
    template <class Type>
    friend std::ostream& operator<<( std::ostream& os,
                                     const Vector<Type>& vec );

    /// Indexing operator ( read only )
    /// \param i Index
    const T& operator[] ( const std::size_t& i ) const;

    /// Indexing operator ( read / write )
    T& operator[] ( const std::size_t& i );

    /// Copy assignment
    /// \param original Vector to be copied
    /// \return The new Vector
    Vector<T>& operator=( const Vector<T>& original );

    /// Unary +
    /// \return The Vector
    Vector<T> operator+() const;

    /// Unary -
    /// \return The negation of the Vector
    Vector<T> operator-() const;

    /// Binary +
    /// \param v_plus The Vector to be added
    /// \return The sum of this Vector and v_plus
    Vector<T> operator+( const Vector<T>& v_plus ) const;

    /// Binary -
    /// \param v_minus The Vector to be subtracted
    /// \return The subtraction of v_minus from this Vector
    Vector<T> operator-( const Vector<T>& v_minus ) const;

    /// Scalar multiplication
    /// \param scalar The scalar to multiply the Vector by
    /// \return The Vector multiplied by the scalar
    Vector<T> operator*( const T& scalar ) const;
    friend Vector<T> operator*( const T& scalar, Vector<T>& vec )
    {
      return vec * scalar;
    }

    /// Scalar division
    /// \param scalar The scalar to divide the Vector by
    /// \return The Vector divided by the scalar
    Vector<T> operator/( const T& scalar ) const;

    /// Addition assignment
    /// \param v_plus The Vector to be added
    /// \return A reference to the Vector after addition by v_plus
    Vector<T>& operator+=( const Vector<T>& v_plus );

    /// Subtraction assignment
    /// \param v_minus The Vector to be subtracted
    /// \return A reference to the Vector after subtraction by v_minus
    Vector<T>& operator-=( const Vector<T>& v_minus );

    /// Scalar multiplication assignment
    /// \param scalar The scalar to multiply the Vector by
    /// \return A reference to the Vector after multiplication by a scalar
    Vector<T>& operator*=( const T& scalar );

    /// Scalar division assigment
    /// \param scalar The scalar to divide the Vector by
    /// \return A reference to the Vector after division by a scalar
    Vector<T>& operator/=( const T& scalar );

    /// Addition assignment ( add a constant to all elements )
    /// \param add The constant to be added to each element
    /// \return A reference to the Vector after addition by a constant
    Vector<T>& operator+=( const T& add );

    /* ----- Methods ----- */

    /// Size of the Vector
    /// \return The number of elements in the Vector
    std::size_t size() const;

    /// Push back a new element into the end of the Vector
    /// \param new_elem New element to be appended to the end of the Vector
    void push_back( const T& new_elem );

    /// TODO Remove the last element in the Vector 

    /// Resize the Vector
    /// \param size New size of the Vector
    void resize( const std::size_t& size );

    /// Clear all the elements from the Vector
    void clear();

    /// Absolute values
    /// \return Vector absolute of values of the elements
    Vector<double> abs() const;

    /// Swap elements i and j
    /// \param i The index of the first element to swap
    /// \param j The index of the second element to swap
    void swap( const std::size_t& i, const std::size_t& j );

    /// Reverse the order of the elements

    /// Assign new contents to the Vector

    /// Insert new elements into the Vector

    /// Real part of the elements of the Vector

    /// Conjugate of the elements of the Vector

    /// Create a linearly spaced Vector (of doubles) with n elements

    /// Create a nonuniform Vector using a power law with exponent p (p=1 ->linear)

    /// Product of the elements in the Vector (from index start to end)

    /// Sum of the elements in the Vector (from index start to end)

    /// Return the dot product of two Vectors v.dot(w)

    /// Output the Vector to a file

    /// Read the Vector from a file

    /// Maximum element

    /// Minimum element

    /// Maximum element index

    /// Minimum element index



    /* ----- Norms ----- */

    /// L1 norm: sum of absolute values

    /// L2 norm: square root of the sum of the squares

    /// Lp norm: p-th root of the sum of the absolute values to the power p

    /// Inf norm: largest absolute value element (p -> infinity)

    /* ----- Iterators ----- */




  }; // End of class Vector

  template <typename T>
  inline Vector<T>::Vector( const std::size_t& size )
  {
    VECTOR.resize( size );
  }

  template <typename T>
  inline Vector<T>::Vector( const std::size_t& size,
                            const T& elem ) : VECTOR( size, elem )
  {}

  template <typename T>
  inline Vector<T>::Vector( const Vector<T>& source )
  {
    VECTOR = source.VECTOR;
  }

  /* ----- Operator overloading ----- */

  template <class Type>
  inline std::ostream& operator<<( std::ostream& os, const Vector<Type>& vec )
  {
    for (std::size_t i = 0; i < vec.VECTOR.size(); ++i )
    {
      os << vec.VECTOR[ i ] << " ";
    }
    return os;
  }

  template <typename T>
  inline const T& Vector<T>::operator[]( const std::size_t& i ) const
  {
    if ( i < 0 || VECTOR.size() <= i )
    {
      throw Error( "Vector operator [] range error" );
    }
    return VECTOR[ i ];
  }

  template <typename T>
  inline T& Vector<T>::operator[] ( const std::size_t& i )
  {
    if ( i < 0 || VECTOR.size() <= i )
    {
      throw Error( "Vector operator [] range error" );
    }
    return VECTOR[ i ];
  }

  template <typename T>
  inline Vector<T>& Vector<T>::operator=( const Vector<T>& original )
  {
    if ( this == &original ){ return *this; }
    VECTOR = original.VECTOR;
    return *this;
  }

  template <typename T>
  inline Vector<T> Vector<T>::operator+() const
  {
    return *this;
  }

  template <typename T>
  inline Vector<T> Vector<T>::operator-() const
  {
    Vector<T> temp( *this );
    std::transform ( temp.VECTOR.cbegin(), temp.VECTOR.cend(),
                     temp.VECTOR.begin(), std::negate<T>() );
    return temp;
  }

  template <typename T>
  inline Vector<T> Vector<T>::operator+( const Vector<T>& v_plus ) const
  {
    Vector<T> temp( *this );
    if ( v_plus.VECTOR.size() != temp.VECTOR.size() )
    {
      throw Error( "Vector dimension error in + operator." );
    }
    std::transform ( temp.VECTOR.cbegin(), temp.VECTOR.cend(),
                     v_plus.VECTOR.begin(), temp.VECTOR.begin(),
                     std::plus<T>() );
    return temp;
  }

  template <typename T>
  inline Vector<T> Vector<T>::operator-( const Vector<T>& v_minus ) const
  {
    Vector<T> temp( *this );
    if ( v_minus.VECTOR.size() != temp.VECTOR.size() )
    {
      throw Error( "Vector dimension error in - operator." );
    }
    std::transform ( temp.VECTOR.cbegin(), temp.VECTOR.cend(),
                     v_minus.VECTOR.begin(), temp.VECTOR.begin(),
                     std::minus<T>() );
    return temp;
  }

  template <typename T>
  inline Vector<T> Vector<T>::operator*( const T& scalar ) const
  {
    Vector<T> temp( *this );
    std::transform ( temp.VECTOR.cbegin(), temp.VECTOR.cend(),
                     temp.VECTOR.begin(),
                     std::bind1st( std::multiplies<T>(), scalar ) );
    return temp;
  }

  template <typename T>
  inline Vector<T> Vector<T>::operator/( const T& scalar ) const
  {
    Vector<T> temp( *this );
    std::transform ( temp.VECTOR.cbegin(), temp.VECTOR.cend(),
                     temp.VECTOR.begin(),
                     std::bind2nd( std::divides<T>(), scalar ) );
    return temp;
  }

  template <typename T>
  inline Vector<T>& Vector<T>::operator+=( const Vector<T>& v_plus )
  {
    if ( v_plus.VECTOR.size() != VECTOR.size() )
    {
      throw Error( "Vector dimension error in += operator." );
    }
    std::transform ( VECTOR.cbegin(), VECTOR.cend(),
                     v_plus.VECTOR.begin(), VECTOR.begin(),
                     std::plus<T>() );
    return *this;
  }

  template <typename T>
  inline Vector<T>& Vector<T>::operator-=( const Vector<T>& v_minus )
  {
    if ( v_minus.VECTOR.size() != VECTOR.size() )
    {
      throw Error( "Vector dimension error in -= operator." );
    }
    std::transform ( VECTOR.cbegin(), VECTOR.cend(),
                     v_minus.VECTOR.begin(), VECTOR.begin(),
                     std::minus<T>() );
    return *this;
  }

  template <typename T>
  inline Vector<T>& Vector<T>::operator*=( const T& scalar )
  {
    std::transform ( VECTOR.cbegin(), VECTOR.cend(), VECTOR.begin(),
                     std::bind1st( std::multiplies<T>(), scalar ) );
    return *this;
  }

  template <typename T>
  inline Vector<T>& Vector<T>::operator/=( const T& scalar )
  {
    std::transform ( VECTOR.cbegin(), VECTOR.cend(), VECTOR.begin(),
                     std::bind2nd( std::divides<T>(), scalar ) );
    return *this;
  }

  template <typename T>
  inline Vector<T>& Vector<T>::operator+=( const T& add )
  {
    std::transform ( VECTOR.cbegin(), VECTOR.cend(), VECTOR.begin(),
                     std::bind1st( std::plus<T>(), add ) );
    return *this;
  }

  /* ----- Methods ----- */

  template <typename T>
  inline std::size_t Vector<T>::size() const
  {
    return VECTOR.size();
  }

  template <typename T>
  inline void Vector<T>::push_back( const T& new_elem )
  {
    VECTOR.push_back( new_elem );
  }

  template <typename T>
  inline void Vector<T>::resize( const std::size_t& size )
  {
    VECTOR.resize( size );
  }

  template <typename T>
  inline void Vector<T>::clear()
  {
    VECTOR.clear();
  }

  template <typename T>
  inline Vector<double> Vector<T>::abs() const
  {
    Vector<double> abs_vals( size() );
    for (size_t i=0; i < size(); ++i)
    {
      abs_vals.VECTOR[ i ] = std::abs( VECTOR[ i ] ) ;
    }
    return abs_vals;
  }

  template <typename T>
  inline void Vector<T>::swap( const std::size_t& i, const std::size_t& j )
  {
    if ( i<0 || size()<=i )	{ throw Error( "swap method Vector range error" );}
    if ( j<0 || size()<=j )	{ throw Error( "swap method Vector range error" );}
    if ( i == j ) { return; }
    std::swap<T>( VECTOR[ i ], VECTOR[ j ] );
  }

}  // End of namespace Luna

#endif
