/// \file Tridiagonal.h
/// A templated Tridiagonal matrix class which constructs a matrix as
/// three Vectors for the main, sub and super diagonals.

#ifndef TRIDIAGONAL_H
#define TRIDIAGONAL_H

#include <vector>
#include <random>
#include <chrono>
#include <complex>

#include "Error.h"
#include "Vector.h"
#include "Matrix.h"

namespace Luna
{

  // Templated Tridiagonal matrix class
  template <class T>

  /// A Tridiagonal matrix class for use with double and std::complex<double>
  class Tridiagonal
  {
  private:
    Vector<T> a;                // Subdiagonal of the matrix (size N-1)
    Vector<T> b;                // Main diagonal of the matrix (size N)
    Vector<T> c;                // Superdiagonal of the matrix (size N-1)
    std::size_t N;              // Size of the matrix (N x N) and main diagonal

  public:
    /// Constructor for an empty Tridiagonal matrix of unspecified size
    Tridiagonal() : N( 0 ) {}

    /// Constructor for a Tridiagonal matrix of specified size
    /// \param size Size of the main diagonal of the matrix
    Tridiagonal( const std::size_t& size );

    /// Constructor for a Tridiagonal matrix using specified Vectors
    /// \param sub The subdiagonal of the matrix as a Vector
    /// \param main The main diagonal of the matrix as a Vector
    /// \param super The superdiagonal of the matrix as a Vector
    Tridiagonal( const Vector<T>& sub, const Vector<T>& main,
                 const Vector<T>& super );

    /// Constructor for a Tridiagonal matrix using specified elements
    /// \param sub_elem The element to populate the subdiagonal
    /// \param main_elem The element to populate the main diagonal
    /// \param super_elem The element to populate the superdiagonal
    /// \param n The size of the main diagonal
    Tridiagonal( const T& sub_elem, const T& main_elem, const T& super_elem,
                 const std::size_t n );

    /// Destructor
    ~Tridiagonal(){}

    /* ----- Operator overloading ----- */

    /// Output operator <<
    template <class Type>
    friend std::ostream& operator<<( std::ostream& os,
                                     const Tridiagonal<Type>& m );

    /// Indexing operator ( read only )
    /// \param i Row index
    /// \param j Column index
    const T& operator() ( const std::size_t& i, const std::size_t& j ) const;

    /// Indexing operator ( read / write )
    /// \param i Row index
    /// \param j Column index
    T& operator() ( const std::size_t& i, const std::size_t& j );

    /// Copy assignment
    /// \param original Tridiagonal matrix to be copied
    /// \return The Tridiagonal matrix to which the original matrix is assigned
    Tridiagonal<T>& operator=( const Tridiagonal<T>& original );

    /// Unary +
    /// \return The Tridiagonal matrix
    Tridiagonal<T> operator+() const;

    /// Unary -
    /// \return The negation of the Tridiagonal matrix
    Tridiagonal<T> operator-() const;

    /// Binary +
    /// \param m_plus The Tridiagonal matrix to be added
    /// \return The sum of this Tridiagonal matrix and m_plus
    Tridiagonal<T> operator+( const Tridiagonal<T>& m_plus ) const;

    /// Binary -
    /// \param m_minus The Tridiagonal matrix to be subtracted
    /// \return The subtraction of m_minus from this Tridiagonal matrix
    Tridiagonal<T> operator-( const Tridiagonal<T>& m_minus ) const;

    /// Scalar multiplication
    /// \param scalar The scalar to multiply the Tridiagonal matrix by
    /// \return The Tridiagonal matrix multiplied by the scalar
    Tridiagonal<T> operator*( const T& scalar ) const;
    friend Tridiagonal<T> operator*( const T& scalar, Tridiagonal<T>& mat )
    {
      return mat * scalar;
    }

    /// Scalar division
    /// \param divisor The divisor to divide the Tridiagonal matrix by
    /// \return The Tridiagonal matrix divided by the divisor
    Tridiagonal<T> operator/( const T& divisor ) const;

    /// Addition assignment
    /// \param m_plus The Tridiagonal matrix to be added
    /// \return A reference to the Tridiagonal matrix after addition by m_plus
    Tridiagonal<T>& operator+=( const Tridiagonal<T>& m_plus );

    /// Subtraction assignment
    /// \param m_minus The Tridiagonal matrix to be subtracted
    /// \return A reference to the Tridiagonal matrix after subtraction by
    /// m_minus
    Tridiagonal<T>& operator-=( const Tridiagonal<T>& m_minus );

    /// Scalar multiplication assignment
    /// \param scalar The scalar to multiply the Tridiagonal matrix by
    /// \return A reference to the Tridiagonal matrix after multiplication by
    /// a scalar
    Tridiagonal<T>& operator*=( const T& scalar );

    /// Scalar division assigment
    /// \param divisor The divisor to divide the Tridiagonal matrix by
    /// \return A reference to the Tridiagonal matrix after division by a scalar
    Tridiagonal<T>& operator/=( const T& divisor );

    /// Constant addition assignment
    /// \param add The constant to be added to each element
    /// \return A reference to the Tridiagonal matrix after addition by a
    /// constant
    Tridiagonal<T>& operator+=( const T& add );

    /// Constant subtraction assignment
    /// \param minus The constant to be subtracted from each element
    /// \return A reference to the Tridiagonal matrix after subtraction by a
    /// constant
    Tridiagonal<T>& operator-=( const T& minus );

    /// Tridiagonal matrix Vector multiplication
    /// \param x The Vector which is to be multiplied
    /// \return The result vector T * x
    Vector<T> operator*( Vector<T>& x );

    /* ----- Methods ----- */

    /// Size of the main diagonal of the matrix
    /// \return The number of elements in the main diagonal
    std::size_t size() const;

    /// Set the subdiagonal of the matrix
    /// \param sub The subdiagonal of the matrix as a Vector
    void set_sub( const Vector<T>& sub );

    /// Set the main diagonal of the matrix
    /// \param main The main diagonal of the matrix as a Vector
    void set_main( const Vector<T>& main );

    /// Set the superdiagonal of the matrix
    /// \param super The superdiagonal of the matrix as a Vector
    void set_super( const Vector<T>& super );

    /// Set the subdiagonal of the matrix
    /// \param sub The subdiagonal element of the matrix
    void set_sub( const T& sub_elem );

    /// Set the main diagonal of the matrix
    /// \param main The main diagonal element of the matrix
    void set_main( const T& main_elem );

    /// Set the superdiagonal of the matrix
    /// \param super The superdiagonal element of the matrix
    void set_super( const T& super_elem );

    /// Get the subdiagonal of the matrix
    /// \return The subdiagonal of the matrix as a Vector
    Vector<T> get_sub();

    /// Get the main diagonal of the matrix
    /// \return The main diagonal of the matrix as a Vector
    Vector<T> get_main();

    /// Get the superdiagonal of the matrix
    /// \return The superdiagonal of the matrix as a Vector
    Vector<T> get_super();

    /// Get a row of the tridiagonal matrix as a Vector
    /// \param row The row index
    /// \return A Vector of the row of the Tridiagonal matrix
    Vector<T> get_row( const std::size_t& row );

    /// Replace the current Matrix with its transpose
    void transpose_in_place();

    /// Return the transpose of the Matrix
    /// \return The transpose of the Matrix
    Tridiagonal<T> transpose() const;

    /// Resize the Tridiagonal matrix
    /// \param size The new size of the main diagonal of the matrix
    void resize( const std::size_t& size );

    /// Take the complex conjugate of each element of the Tridiagonal matrix
    /// \return The complex conjugate Tridiagonal matrix
    Tridiagonal<T> conjugate() const;

    /// Convert the Tridiagonal matrix to a dense Matrix
    /// \return The dense Matrix form of the Tridiagonal matrix
    Matrix<T> convert();

    /* ----- Solve linear systems ----- */

    /// Solve the Tridiagonal system of equations Tu=r where u and r are Vectors
    /// \param r The right-side Vector of the system of equations
    /// \return The solution Vector u
    Vector<T> solve( const Vector<T>& r );

    /// Solve the Tridiagonal system TU=R where U and R are Matrices
    /// \param R The right-side Matrix of the system of equations
    /// \return The solution Matrix U
    Matrix<T> solve( const Matrix<T>& R );

    /* ----- Determinant ----- */

    /// Calculate the determinant of the matrix
    /// \return The determinant of the Tridiagonal matrix
    T det();

  };	// End of class Tridiagonal

  template <typename T>
  inline Tridiagonal<T>::Tridiagonal( const std::size_t& size )
  {
    a.resize( size - 1 );
    b.resize( size );
    c.resize( size - 1 );
    N = size;
  }

  template <typename T>
  inline Tridiagonal<T>::Tridiagonal( const Vector<T>& sub,
                                      const Vector<T>& main,
                                      const Vector<T>& super )
  {
    if ( main.size() != sub.size() + 1 || main.size() != super.size() + 1 )
    {
      std::string problem;
      problem = "Tridiagonal constructor error: ";
      problem += "size of Vectors are incompatible.";
      throw Error( problem );
    }
    a = sub;
    b = main;
    c = super;
    N = main.size();
  }

  template <typename T>
  inline Tridiagonal<T>::Tridiagonal( const T& sub_elem, const T& main_elem,
                                      const T& super_elem, const std::size_t n )
  {
    a.assign( n - 1, sub_elem );
    b.assign( n, main_elem );
    c.assign( n - 1, super_elem );
    N = n;
  }

  /* ----- Operator overloading ----- */

  template <class Type>
  inline std::ostream& operator<<( std::ostream& os,
                                   const Tridiagonal<Type>& m )
  {
    os << std::endl;
    if ( m.N == 0 )	{ return os; }
    if ( m.N > 0 )
    {
      for (std::size_t j = 0; j < m.N; ++j )
      {
        std::string tabs;
        for ( std::size_t i = 0; i < j; i++ ){ tabs += "\t"; }

        if ( j == 0 )
        {
          os << "\t" << m.b[ 0 ] << "\t" << m.c[ 0 ] << std::endl;
        }
        else if ( j == m.N - 1 )
        {
          os << tabs << m.a[ m.N - 2 ] << "\t" << m.b[ m.N - 1 ] << std::endl;
        }
        else
        {
          os << tabs << m.a[ j - 1 ] << "\t" << m.b[ j ]
             << "\t" << m.c[ j ] << std::endl;
        }
      }
    }
    return os;
  }

  template <typename T>
  inline const T& Tridiagonal<T>::operator() ( const std::size_t& i,
                                               const std::size_t& j ) const
  {
    if ( i<0 || N<=i )	{ throw Error( "Tridiagonal range error: row." ); }
    if ( j<0 || N<=j )	{ throw Error( "Tridiagonal range error: column." ); }
    // First row
    if ( i == 0 )
    {
      if ( j == 0 ) { return b[ 0 ]; }
      if ( j == 1 ) { return c[ 0 ]; }
      else{
        throw Error( "Tridiagonal error: off diagonal entry." );
      }
    }
    // Last row
    if ( i == N - 1 )
    {
      if ( j == N - 1 ) { return b[ N - 1 ]; }
      if ( j == N - 2 ) { return a[ N - 2 ]; }
      else {
        throw Error( "Tridiagonal error: off diagonal entry." );
      }
    }
    // Middle rows
    if ( i > 0 && i < N - 1 )
    {
      if ( i == j ) { return b[ j ]; }
      if ( i == j + 1 ) { return a[ j ]; }
      if ( i == j - 1 ) { return c[ i ]; }
      else {
        throw Error( "Tridiagonal error: off diagonal entry." );
      }
    }
  }

  template <typename T>
  inline T& Tridiagonal<T>::operator() ( const std::size_t& i,
                                         const std::size_t& j )
  {
    if ( i<0 || N<=i )	{ throw Error( "Tridiagonal range error: row." );}
    if ( j<0 || N<=j )	{ throw Error( "Tridiagonal range error: column." );}
    // First row
    if ( i == 0 )
    {
      if ( j == 0 ) { return b[ 0 ]; }
      if ( j == 1 ) { return c[ 0 ]; }
      else{
        throw Error( "Tridiagonal error: off diagonal entry." );
      }
    }
    // Last row
    if ( i == N - 1 )
    {
      if ( j == N - 1 ) { return b[ N - 1 ]; }
      if ( j == N - 2 ) { return a[ N - 2 ]; }
      else {
        throw Error( "Tridiagonal error: off diagonal entry." );
      }
    }
    // Middle rows
    if ( i > 0 && i < N - 1 )
    {
      if ( i == j ) { return b[ j ]; }
      if ( i == j + 1 ) { return a[ j ]; }
      if ( i == j - 1 ) { return c[ i ]; }
      else {
        throw Error( "Tridiagonal error: off diagonal entry." );
      }
    }
  }

  template <typename T>
  inline Tridiagonal<T>& Tridiagonal<T>::operator=(
                                              const Tridiagonal<T>& original )
  {
    if ( this == &original )
    {
      return *this;
    }
    a = original.a;
    b = original.b;
    c = original.c;
    N = original.N;
    return *this;
  }

  template <typename T>
  inline Tridiagonal<T> Tridiagonal<T>::operator+() const
  {
    return *this;
  }

  template <typename T>
  inline Tridiagonal<T> Tridiagonal<T>::operator-() const
  {
    Tridiagonal<T> temp( *this );
    temp.a = - a;
    temp.b = - b;
    temp.c = - c;
    return temp;
  }

  template <typename T>
  inline Tridiagonal<T> Tridiagonal<T>::operator+(
                                          const Tridiagonal<T>& m_plus ) const
  {
    Tridiagonal<T> temp( *this );
    if ( m_plus.N != N  )
    {
      throw Error( "Tridiagonal matrix dimension error in + operator." );
    }
    temp.a = a + m_plus.a;
    temp.b = b + m_plus.b;
    temp.c = c + m_plus.c;
    return temp;
  }

  template <typename T>
  inline Tridiagonal<T> Tridiagonal<T>::operator-(
                                          const Tridiagonal<T>& m_minus ) const
  {
    Tridiagonal<T> temp( *this );
    if ( m_minus.N != N )
    {
      throw Error( "Tridiagonal matrix dimension error in - operator." );
    }
    temp.a = a - m_minus.a;
    temp.b = b - m_minus.b;
    temp.c = c - m_minus.c;
    return temp;
  }

  template <typename T>
  inline Tridiagonal<T> Tridiagonal<T>::operator*( const T& scalar ) const
  {
    Tridiagonal<T> temp( *this );
    temp.a = a * scalar;
    temp.b = b * scalar;
    temp.c = c * scalar;
    return temp;
  }

  template <typename T>
  inline Tridiagonal<T> Tridiagonal<T>::operator/( const T& divisor ) const
  {
    Tridiagonal<T> temp( *this );
    temp.a = a / divisor;
    temp.b = b / divisor;
    temp.c = c / divisor;
    return temp;

  }

  template <typename T>
  inline Tridiagonal<T>& Tridiagonal<T>::operator+=(
                                                  const Tridiagonal<T>& m_plus )
  {
    if ( m_plus.N != N )
    {
      throw Error( "Tridiagonal matrix dimension error in += operator." );
    }
    a = a + m_plus.a;
    b = b + m_plus.b;
    c = c + m_plus.c;
    return *this;
  }

  template <typename T>
  inline Tridiagonal<T>& Tridiagonal<T>::operator-=(
                                                 const Tridiagonal<T>& m_minus )
  {
    if ( m_minus.N != N )
    {
      throw Error( "Tridiagonal matrix dimension error in -= operator." );
    }
    a = a - m_minus.a;
    b = b - m_minus.b;
    c = c - m_minus.c;
    return *this;
  }

  template <typename T>
  inline Tridiagonal<T>& Tridiagonal<T>::operator*=( const T& scalar )
  {
    a *= scalar;
    b *= scalar;
    c *= scalar;
    return *this;
  }

  template <typename T>
  inline Tridiagonal<T>& Tridiagonal<T>::operator/=( const T& divisor )
  {
    a /= divisor;
    b /= divisor;
    c /= divisor;
    return *this;
  }

  template <typename T>
  inline Tridiagonal<T>& Tridiagonal<T>::operator+=( const T& add )
  {
    a += add;
    b += add;
    c += add;
    return *this;
  }

  template <typename T>
  inline Tridiagonal<T>& Tridiagonal<T>::operator-=( const T& minus )
  {
    a -= minus;
    b -= minus;
    c -= minus;
    return *this;
  }

  template <typename T>
  inline Vector<T> Tridiagonal<T>::operator*( Vector<T>& x )
  {
    Vector<T> result( N );
    result[ 0 ] = b[ 0 ] * x[ 0 ] + c[ 0 ] * x[ 1 ];
    for ( std::size_t i = 1; i < N - 1; i++ )
    {
      result[ i ] = a[ i - 1 ] * x[ i - 1 ] + b[ i ] * x[ i ]
                  + c[ i ] * x[ i + 1 ];
    }
    result[ N - 1 ] = a[ N - 2 ] * x[ N - 2 ] + b[ N - 1 ] * x[ N - 1 ];
    return result;
  }


  /* ----- Methods ----- */

  template <typename T>
  inline std::size_t Tridiagonal<T>::size() const
  {
    return N;
  }

  template <typename T>
  inline void Tridiagonal<T>::set_sub( const Vector<T>& sub )
  {
    if ( sub.size() != N - 1 )
    {
      throw Error( "Tridiagonal error: subdiagonal wrong size." );
    }
    a = sub;
  }

  template <typename T>
  inline void Tridiagonal<T>::set_main( const Vector<T>& main )
  {
    if ( main.size() != N )
    {
      throw Error( "Tridiagonal error: main diagonal wrong size." );
    }
    b = main;
  }

  template <typename T>
  inline void Tridiagonal<T>::set_super( const Vector<T>& super )
  {
    if ( super.size() != N - 1 )
    {
      throw Error( "Tridiagonal error: superdiagonal wrong size." );
    }
    c = super;
  }

  template <typename T>
  inline void Tridiagonal<T>::set_sub( const T& sub_elem )
  {
    a.assign( N, sub_elem );
  }

  template <typename T>
  inline void Tridiagonal<T>::set_main( const T& main_elem )
  {
    b.assign( N, main_elem );
  }

  template <typename T>
  inline void Tridiagonal<T>::set_super( const T& super_elem )
  {
    c.assign( N, super_elem );
  }

  template <typename T>
  inline Vector<T> Tridiagonal<T>::get_sub()
  {
    Vector<T> temp( a );
    return temp;
  }

  template <typename T>
  inline Vector<T> Tridiagonal<T>::get_main()
  {
    Vector<T> temp( b );
    return temp;
  }

  template <typename T>
  inline Vector<T> Tridiagonal<T>::get_super()
  {
    Vector<T> temp( c );
    return temp;
  }

  template <typename T>
  inline Vector<T> Tridiagonal<T>::get_row( const std::size_t& row )
  {
    if ( row < 0 || N <= row ){
      throw Error( "Tridiagonal range error in get_row method." );
    }
    Vector<T> temp( N );
    if ( row == 0 )
    {
      temp[ 0 ] = b[ 0 ];
      temp[ 1 ] = c[ 0 ];
    }
    else if ( row == N - 1 )
    {
      temp[ N - 1 ] = b[ N - 1 ];
      temp[ N - 2 ] = a[ N - 2 ];
    }
    else
    {
      temp[ row - 1 ] = a[ row - 1 ];
      temp[ row ]     = b[ row ];
      temp[ row + 1 ] = c[ row ];
    }
    return temp;
  }

  template <typename T>
  inline void Tridiagonal<T>::transpose_in_place()
  {
    Vector<T> temp( a );
    a = c;
    c = temp;
  }

  template <typename T>
  inline Tridiagonal<T> Tridiagonal<T>::transpose() const
  {
    Tridiagonal<T> temp( *this );
    temp.transpose_in_place();
    return temp;
  }

  template <typename T>
  inline void Tridiagonal<T>::resize( const std::size_t& size )
  {
    N = size;
    a.resize( size - 1 );
    b.resize( size );
    c.resize( size - 1 );
  }

  template <typename T>
  inline Tridiagonal<T> Tridiagonal<T>::conjugate() const
  {
    Tridiagonal<T> temp( *this );
    temp.a = temp.a.conjugate();
    temp.b = temp.b.conjugate();
    temp.c = temp.c.conjugate();
    return temp;
  }

  template <typename T>
  inline Matrix<T> Tridiagonal<T>::convert()
  {
    Matrix<T> dense( N, N, 0 );
    if ( N == 0 ) { throw Error( "Tridiagonal error: zero size matrix." ); }
    if ( N == 1 ) {
      dense( 0, 0 ) = b[ 0 ];
      dense( 0, 1 ) = c[ 0 ];
    } else if ( N == 2 ) {
      dense( 0, 0 ) = b[ 0 ];
      dense( 0, 1 ) = c[ 0 ];
      dense( N - 1, N - 2 ) = a[ N - 2 ];
      dense( N - 1, N - 1 ) = b[ N - 1 ];
    } else {
      dense( 0, 0 ) = b[ 0 ];
      dense( 0, 1 ) = c[ 0 ];

      for ( std::size_t i = 1; i < N - 1; i++ )
      {
        dense( i, i - 1 ) = a[ i - 1 ];
        dense( i, i ) = b[ i ];
        dense( i, i + 1 ) = c[ i ];
      }

      dense( N - 1, N - 2 ) = a[ N - 2 ];
      dense( N - 1, N - 1 ) = b[ N - 1 ];
    }
    return dense;
  }


  /* ----- Solve linear systems ----- */

  template <typename T>
  inline Vector<T> Tridiagonal<T>::solve( const Vector<T>& r )
  {
    if ( r.size() != N ) {
      throw Error( "Tridiagonal solve error: Vector incorrect size." );
    }

    Vector<T> u( N );
    Vector<T> a_temp( a );
    Vector<T> c_temp( c );
    a_temp.push_front( 0.0 );
    c_temp.push_back( 0.0 );

    T beta( b[ 0 ] );
    Vector<T> gamma( N );

    if ( b[ 0 ] == 0.0 ) {
      throw Error( "Tridiagonal solve error: zero on leading diagonal." );
    }

    u[ 0 ] = r[ 0 ] / beta;
    for ( int j = 1; j < N; j++ )
    {
      gamma[ j ] = c_temp[ j - 1 ] / beta;
      beta = b[ j ] - a_temp[ j ] * gamma[ j ];
      if ( beta == 0.0 ) {
        throw Error( "Tridiagonal solve error: zero pivot." );
      }
      u[ j ] = ( r[ j ] - a_temp[ j ] * u[ j - 1 ] ) / beta;
    }

    for ( int j = (N - 2); j >= 0; j-- )
    {
      u[ j ] -= gamma[ j + 1 ] * u[ j + 1 ];
    }

    return u;
  }

  template <typename T>
  inline Matrix<T> Tridiagonal<T>::solve( const Matrix<T>& R )
  {
    if ( R.rows() != N ) {
      throw Error( "Tridiagonal solve error: Matrix rows incorrect size." );
    }
    Matrix<T> U( R );
    Vector<T> u( N );
    for ( std::size_t j = 0; j < R.cols(); ++j )
    {
      u = U.get_col( j );
      u = solve( u );
      U.set_col( j, u );
    }
    return U;
  }

  /* ----- Determinant ----- */

  template <typename T>
  inline T Tridiagonal<T>::det()
  {
    Vector<T> f( N + 1 );
    f[ 0 ] = 1;
    f[ 1 ] = b[ 0 ] * f[ 0 ];

    for ( std::size_t j = 2; j < N + 1; j++ )
    {
      f[ j ] = b[ j - 1 ] * f[ j - 1 ] - a[ j - 2 ] * c[ j - 2 ] * f[ j - 2];
    }
    return f[ N ];
  }

  /* ----- Inverse ----- */

}  // End of namespace Luna

#endif
