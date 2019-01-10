/// \file Vector.h
/// A templated numeric Vector class that encapsulates std::vector.

#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <algorithm>
#include <functional>
#include <numeric>
#include <complex>

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
    /// \return A constant reference to the element located at the given index
    const T& operator[] ( const std::size_t& i ) const;

    /// Indexing operator ( read / write )
    /// \param i Index
    /// \return A reference to the element located at the given index
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
    /// \param scalar The divisor to divide the Vector by
    /// \return The Vector divided by the divisor
    Vector<T> operator/( const T& divisor ) const;

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
    /// \param scalar The divisor to divide the Vector by
    /// \return A reference to the Vector after division by a divisor
    Vector<T>& operator/=( const T& divisor );

    /// Constant addition assignment
    /// \param add The constant to be added to each element
    /// \return A reference to the Vector after addition by a constant
    Vector<T>& operator+=( const T& add );

    /// Constant subtraction assignment
    /// \param minus The constant to be subtracted from each element
    /// \return A reference to the Vector after subtraction by a constant
    Vector<T>& operator-=( const T& minus );

    /* ----- Methods ----- */

    /// Size of the Vector
    /// \return The number of elements in the Vector
    std::size_t size() const;

    /// Push back a new element into the end of the Vector
    /// \param new_elem New element to be appended to the end of the Vector
    void push_back( const T& new_elem );

    /// Remove the last element in the Vector
    void pop_back();

    /// Returns whether the Vector is empty
    /// \return true if the Vector is size 0 and false otherwise
    bool empty() const;

    /// Resize the Vector
    /// \param size New size of the Vector
    void resize( const std::size_t& size );

    /// Request a change in capacity
    /// \param size Minimum capacity for the Vector
    void reserve( const std::size_t& size );

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
    void reverse();

    /// Assign new contents to the Vector
    /// \param n The number of elements to assign
    /// \param elem The element to be assigned
    void assign( const std::size_t n, const T elem );

    /// Real part of the elements of the Vector
    /// \return Vector of real parts of the elements
    Vector<double> real() const;

    /// Conjugate of the elements of the Vector
    /// \return Vector of conjugates of the elements
    Vector<T> conjugate() const;

    /// Create a linearly spaced Vector (of doubles) with n elements
    /// \param a Start value
    /// \param b End value
    /// \param n Number of elements
    void linspace( const double& a, const double& b, const std::size_t& n );

    /// Create a nonuniform Vector using a power law with exponent p (1->linear)
    /// \param a Start value
    /// \param b End value
    /// \param n Number of elements
    /// \param p Exponent
    void powspace( const double& a, const double& b,
                   const std::size_t& n, const double& p);

    /// Product of the elements in the Vector (from index start to end)
    /// \param start Start index
    /// \param end End index
    /// \return Product of the elements in the Vector
    T product( const std::size_t& start, const std::size_t& end );

    // Product of the elements in the Vector
    /// \return Product of the elements in the Vector
    T product();

    /// Sum of the elements in the Vector (from index start to end)
    /// \param start Start index
    /// \param end End index
    /// \return Sum of the elements in the Vector
    T sum( const std::size_t& start, const std::size_t& end );

    /// Sum of the elements in the Vector
    /// \return Sum of the elements in the Vector
    T sum();

    /// Dot product of two Vectors v.dot(w)
    /// \param w The Vector to be dotted with
    /// \return The dot product of two vectors
    T dot( const Vector<T>& w );

    /// Output the Vector to a file
    /// \param filename The name of the file to output the Vector to
    /// \param precision The precision of the elements to be stored
    void output( std::string filename, int precision = 12 );

    /// Read the Vector from a file
    /// \param filename The name of the file to read the Vector from
    /// \param precision The precision of the elements to be read in
    void read( std::string filename, int precision = 12 );

    /// Maximum element
    /// \return The maximum element in the Vector
    T max();

    /// Minimum element
    /// \return The minimum element in the Vector
    T min();

    /// Maximum element index
    /// \return The index of the largest element (first index if more than one)
    std::size_t max_index();

    /// Minimum element index
    /// \return The index of the smallest element (first index if more than one)
    std::size_t min_index();

    /// Square all the elements in the Vector
    /// \return Vector of squares of the elements
    Vector<T> square() const;

    /// Raise all the elements in the Vector to the power p
    /// \param p The exponent to raise each of the elements to
    /// \return Vector of the elements raised to the power p
    Vector<T> power( const double& p ) const;

    /* ----- Norms ----- */

    /// L1 norm: sum of absolute values
    /// \return The L1 norm of the Vector
    double norm_1() const;

    /// L2 norm: square root of the sum of the squares
    /// \return The L2 norm of the Vector
    double norm_2() const;

    /// Lp norm: p-th root of the sum of the absolute values to the power p
    /// \param p The exponent for the p-norm
    /// \return The Lp norm of the Vector
    double norm_p( const double& p ) const;

    /// Inf norm: largest absolute value element (p -> infinity)
    /// \return The infinity norm of the Vector
    double norm_inf() const;

    /* ----- Iterators ----- */

    //template <typename Type>
    typedef typename std::vector<T>::iterator iter;
    typedef typename std::vector<T>::const_iterator citer;
    typedef typename std::vector<T>::reverse_iterator riter;
    typedef typename std::vector<T>::const_reverse_iterator criter;

    /// Pass through the std::vector begin iterator
    /// \return An iterator pointing to the first element
    iter begin()
    {
      return VECTOR.begin();
    }

    /// Pass through the std::vector begin reverse iterator
    /// \return A reverse iterator pointing to the last element
    riter rbegin()
    {
      return VECTOR.rbegin();
    }

    /// Pass through the std::vector begin constant iterator
    /// \return A constant iterator pointing to the first element
    citer begin() const
    {
      return VECTOR.begin();
    }

    /// Pass through the std::vector begin constant reverse iterator
    /// \return A constant reverse iterator pointing to the last element
    criter rbegin() const
    {
      return VECTOR.rbegin();
    }

    /// Pass through the std::vector end iterator
    /// \return An iterator referring to the past-the-end element
    iter end()
    {
      return VECTOR.end();
    }

    /// Pass through the std::vector end reverse iterator
    /// \return A reverse iterator pointing to the element before the first
    riter rend()
    {
      return VECTOR.rend();
    }

    /// Pass through the std::vector end constant iterator
    /// \return A constant iterator referring to the past-the-end element
    citer end() const
    {
      return VECTOR.end();
    }

    /// Pass through the std::vector end constant reverse iterator
    /// \return A constant reverse iterator pointing to the element before the first
    criter rend() const
    {
      return VECTOR.rend();
    }

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
      os << std::setw( 6 ) << std::setprecision( 2 ) << std::fixed
         << vec.VECTOR[ i ] << " ";
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
                     std::bind( std::multiplies<T>(),
                                std::placeholders::_1, scalar ) );
    return temp;
  }

  template <typename T>
  inline Vector<T> Vector<T>::operator/( const T& divisor ) const
  {
    Vector<T> temp( *this );
    std::transform ( temp.VECTOR.cbegin(), temp.VECTOR.cend(),
                     temp.VECTOR.begin(),
                     std::bind( std::divides<T>(),
                                std::placeholders::_1, divisor ) );
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
                     std::bind( std::multiplies<T>(),
                                std::placeholders::_1, scalar ) );
    return *this;
  }

  template <typename T>
  inline Vector<T>& Vector<T>::operator/=( const T& divisor )
  {
    std::transform ( VECTOR.cbegin(), VECTOR.cend(), VECTOR.begin(),
                     std::bind( std::divides<T>(),
                                std::placeholders::_1, divisor ) );
    return *this;
  }

  template <typename T>
  inline Vector<T>& Vector<T>::operator+=( const T& add )
  {
    std::transform ( VECTOR.cbegin(), VECTOR.cend(), VECTOR.begin(),
                     std::bind( std::plus<T>(), std::placeholders::_1, add ) );
    return *this;
  }

  template <typename T>
  inline Vector<T>& Vector<T>::operator-=( const T& minus )
  {
    std::transform ( VECTOR.cbegin(), VECTOR.cend(), VECTOR.begin(),
                     std::bind( std::minus<T>(), std::placeholders::_1, minus ));
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
  inline void Vector<T>::pop_back()
  {
    VECTOR.pop_back();
  }

  template <typename T>
  inline bool Vector<T>::empty() const
  {
    return VECTOR.empty();
  }

  template <typename T>
  inline void Vector<T>::resize( const std::size_t& size )
  {
    VECTOR.resize( size );
  }

  template <typename T>
  inline void Vector<T>::reserve( const std::size_t& size )
  {
    VECTOR.reserve( size );
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

  template <typename T>
  inline void Vector<T>::reverse()
  {
    std::reverse( VECTOR.begin(), VECTOR.end() );
  }

  template <typename T>
  inline void Vector<T>::assign( const std::size_t n, const T elem )
  {
    VECTOR.assign( n, elem );
  }

  template <typename T>
  inline Vector<double> Vector<T>::real() const
  {
    Vector<double> real_part( size() );
    for (size_t i=0; i < size(); ++i)
    {
      real_part[ i ] = std::real( VECTOR[ i ] ) ;
    }
    return real_part;
  }

  template <typename T>
  inline Vector<T> Vector<T>::conjugate() const
  {
     Vector<T> conjugate( *this );
     for (size_t i=0; i < size(); ++i)
     {
       conjugate.VECTOR[ i ] = std::conj( VECTOR[ i ] ) ;
     }
     return conjugate;
  }

  template <>
	inline void Vector<double>::linspace( const double& a, const double& b,
                                        const std::size_t& n )
	{
		VECTOR.resize( n );
		const double h = ( b - a ) / (n - 1)  ;
		for ( std::size_t i=0; i < n; ++i )
		{
			VECTOR[ i ] = a + h * i;
		}
	}

  template <>
  inline void Vector<double>::powspace( const double& a, const double& b,
                                     const std::size_t& n, const double& p)
  {
    VECTOR.resize( n );
    for ( std::size_t i = 0; i < n; ++i )
    {
      VECTOR[ i ] = a + ( b - a ) * std::pow( ( double )i / ( n - 1 ), p );
    }
  }

  template <typename T>
  inline T Vector<T>::product( const std::size_t& start,
                               const std::size_t& end )
  {
    if ( start > end )	{ throw Error( "Vector product: start > end" );}
    if ( start<0 || size()<=start )	{
      throw Error( "Vector product method range error" );}
    if ( end<0 || size()<=end )	{
      throw Error( "Vector product method range error" );}
    T prod( VECTOR[ start ] );
    for ( std::size_t i=start+1; i<=end; ++i )
    {
      prod *= VECTOR[ i ];
    }
    return prod;
  }

  template <typename T>
  inline T Vector<T>::product()
  {
    return this->product( 0, size() - 1 );
  }

  template <typename T>
  inline T Vector<T>::sum( const std::size_t& start, const std::size_t& end )
  {
    if ( start > end )	{ throw Error( "Vector sum: start > end" );}
    if ( start<0 || size()<=start )	{
      throw Error( "Vector sum method range error" );}
    if ( end<0 || size()<=end )	{
      throw Error( "Vector sum method range error" );}
    T Sum( VECTOR[ start ] );
    for ( std::size_t i=start+1; i<=end; ++i )
    {
      Sum += VECTOR[ i ];
    }
    return Sum;
  }

  template <typename T>
  inline T Vector<T>::sum()
  {
    return this->sum( 0, size() - 1 );
  }

  template <typename T>
  inline T Vector<T>::dot( const Vector<T>& w )
  {
    if ( size() != w.size() )	{ throw Error( "Vector dot product: size error" );}
    T init( 0.0 );
    return std::inner_product ( VECTOR.cbegin(), VECTOR.cend(),
                     w.VECTOR.begin(), init );
  }

  template <typename T>
  inline void Vector<T>::output( std::string filename, int precision )
  {
    std::ofstream out;
    out.open( filename.c_str() );
    out.precision( precision );
    out.setf( std::ios::showpoint );
    out.setf( std::ios::showpos );
    out.setf( std::ios::scientific );
    for ( std::size_t i = 0; i < size(); ++i )
    {
      out << VECTOR[ i ] << std::endl;
    }
  }

  template <typename T>
  inline void Vector<T>::read( std::string filename, int precision )
  {
    std::ifstream inFile;
		inFile.open( filename.c_str() );
		inFile.precision( precision );
		inFile.setf( std::ios::showpoint );
		inFile.setf( std::ios::showpos );
		inFile.setf( std::ios::scientific );
		if (!inFile) { throw Error( "Vector.read() unable to open file" ); }
    VECTOR.resize( 0 );
    T val;
    while ( inFile >> val )
    {
      VECTOR.push_back( val );
    }
    inFile.close();
  }

  template <typename T>
  inline T Vector<T>::max()
  {
    return VECTOR[ this-> max_index() ];
  }

  template <typename T>
  inline T Vector<T>::min()
  {
    return VECTOR[ this-> min_index() ];
  }

  template <typename T>
  inline std::size_t Vector<T>::max_index()
  {
    return std::distance( VECTOR.begin(),
                           std::max_element( VECTOR.begin(), VECTOR.end() ) );
  }

  template <typename T>
  inline std::size_t Vector<T>::min_index()
  {
    return std::distance( VECTOR.begin(),
                           std::min_element( VECTOR.begin(), VECTOR.end() ) );
  }

  template <typename T>
  inline Vector<T> Vector<T>::square() const
  {
    Vector<T> square;
    for (size_t i=0; i < size(); ++i)
    {
      square.VECTOR.push_back( VECTOR[ i ] * VECTOR[ i ] );
    }
    return square;
  }

  template <typename T>
  inline Vector<T> Vector<T>::power( const double& p ) const
  {
    Vector<T> power;
    for (size_t i=0; i < size(); ++i)
    {
      power.VECTOR.push_back( std::pow( VECTOR[ i ], p ) );
    }
    return power;
  }

  /* ----- Norms ----- */

  template <typename T>
  inline double Vector<T>::norm_1() const
  {
    return ( this->abs() ).sum();
  }

  template <typename T>
  inline double Vector<T>::norm_2() const
  {
    return std::sqrt( ( this->abs() ).square().sum() );
  }

  template <typename T>
  inline double Vector<T>::norm_p( const double& p ) const
  {
    return std::pow( ( this->abs() ).power( p ).sum(), 1.0 / p );
  }

  template <typename T>
  inline double Vector<T>::norm_inf() const
  {
    return ( this->abs() ).max();
  }

}  // End of namespace Luna

#endif
