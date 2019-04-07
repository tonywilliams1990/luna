/// \file BandedMatrix.h
/// A templated banded matrix class.

#ifndef BANDEDMATRIX_H
#define BANDEDMATRIX_H

#include <vector>
#include <random>
#include <chrono>
#include <complex>
#include <algorithm>

#include "Error.h"
#include "Vector.h"
#include "Matrix.h"

namespace Luna
{

  // Templated BandedMatrix class
  template <class T>

  /// A banded matrix class for use with double and std::complex<double>
  class BandedMatrix
  {
  private:
    std::size_t N;        // Size of the matrix (N x N) and main diagonal
    std::size_t M1;       // Number of bands below the main diagonal
    std::size_t M2;       // Number of bands above the main diagonal

    Matrix<T> COMPACT;    // Banded matrix in compact form

    void decompose( Matrix<T>& au, Matrix<T>& al,
                    Vector<int>& index, double& d );

  public:

    /// Constructor for a BandedMatrix of specfied size and number of bands
    /// \param n The size of the matrix (n x n)
    /// \param m1 The number of bands below the main diagonal
    /// \param m2 The number of bands above the main diagonal
    BandedMatrix( const std::size_t& n, const std::size_t& m1,
                  const std::size_t& m2 );

    /// Constructor for a BandedMatrix from a dense Matrix
    /// \param mat The dense Matrix used to create the BandedMatrix
    /// \param m1 The number of bands below the main diagonal
    /// \param m2 The number of bands above the main diagonal
    BandedMatrix( const Matrix<T>& mat, const std::size_t& m1,
                  const std::size_t& m2 );

    /// Destructor
    ~BandedMatrix(){}

    /* ----- Operator overloading ----- */

    /// Output operator <<
    template <class Type>
    friend std::ostream& operator<<( std::ostream& os,
                                     const BandedMatrix<Type>& m );

    /// Indexing operator ( read only )
    /// \param i Row index
    /// \param j Column index
    const T& operator() ( const int& i, const int& j ) const;

    /// Indexing operator ( read / write )
    /// \param i Row index
    /// \param j Column index
    T& operator() ( const int& i, const int& j );

    /// Copy assignment
    /// \param original BandedMatrix to be copied
    /// \return The BandedMatrix to which the original matrix is assigned
    BandedMatrix<T>& operator=( const BandedMatrix<T>& original );

    /// Unary +
    /// \return The BandedMatrix
    BandedMatrix<T> operator+() const;

    /// Unary -
    /// \return The negation of the BandedMatrix
    BandedMatrix<T> operator-() const;

    /// Binary +
    /// \param m_plus The BandedMatrix to be added
    /// \return The sum of this BandedMatrix and m_plus
    BandedMatrix<T> operator+( const BandedMatrix<T>& m_plus ) const;

    /// Binary -
    /// \param m_minus The BandedMatrix to be subtracted
    /// \return The subtraction of m_minus from this BandedMatrix
    BandedMatrix<T> operator-( const BandedMatrix<T>& m_minus ) const;

    /// Scalar multiplication
    /// \param scalar The scalar to multiply the BandedMatrix by
    /// \return The BandedMatrix multiplied by the scalar
    BandedMatrix<T> operator*( const T& scalar ) const;
    friend BandedMatrix<T> operator*( const T& scalar, BandedMatrix<T>& mat )
    {
      return mat * scalar;
    }

    /// Scalar division
    /// \param divisor The divisor to divide the BandedMatrix by
    /// \return The BandedMatrix divided by the divisor
    BandedMatrix<T> operator/( const T& divisor ) const;

    /// Addition assignment
    /// \param m_plus The BandedMatrix to be added
    /// \return A reference to the BandedMatrix after addition by m_plus
    BandedMatrix<T>& operator+=( const BandedMatrix<T>& m_plus );

    /// Subtraction assignment
    /// \param m_minus The BandedMatrix to be subtracted
    /// \return A reference to the BandedMatrix after subtraction by m_minus
    BandedMatrix<T>& operator-=( const BandedMatrix<T>& m_minus );

    /// Scalar multiplication assignment
    /// \param scalar The scalar to multiply the BandedMatrix by
    /// \return A reference to the BandedMatrix after multiplication by a scalar
    BandedMatrix<T>& operator*=( const T& scalar );

    /// Scalar division assigment
    /// \param divisor The divisor to divide the BandedMatrix by
    /// \return A reference to the BandedMatrix after division by a scalar
    BandedMatrix<T>& operator/=( const T& divisor );

    /// Constant addition assignment
    /// \param add The constant to be added to each element
    /// \return A reference to the BandedMatrix after addition by a constant
    BandedMatrix<T>& operator+=( const T& add );

    /// Constant subtraction assignment
    /// \param minus The constant to be subtracted from each element
    /// \return A reference to the BandedMatrix after subtraction by a constant
    BandedMatrix<T>& operator-=( const T& minus );

    /// BandedMatrix Vector multiplication (B * x)
    /// \param x The Vector which is to be multiplied
    /// \return The result vector B * x
    Vector<T> operator*( Vector<T>& x );

    /* ----- Methods ----- */

    /// Return the BandedMatrix stored in compact form
    /// \return A Matrix containing the compact form of the BandedMatrix
    Matrix<T> compact();

    /// Return the size of the BandedMatrix
    /// \return The size of the BandedMatrix
    std::size_t size();

    /// Return the number of bands below the main diagonal
    /// \return The number of bands below the main diagonal
    std::size_t size_below();

    /// Return the number of bands above the main diagonal
    /// \return The number of bands above the main diagonal
    std::size_t size_above();

    /// Fill the Matrix with specified elements
    /// \param elem The element to fill the Matrix with
    void fill( const T& elem );

    /* ----- Solve linear systems ----- */

    /// Solve the banded system of equations Ax=b where x and b are Vectors
    /// \param b The right-side Vector of the system of equations
    /// \return The solution Vector x
    Vector<T> solve( const Vector<T>& b );

    /// Solve the banded system of equations AX=B where X and B are Matrices
    /// \param B The right-hand side Matrix of the system of equations
    /// \return The solution Matrix (each column is a solution Vector)
    Matrix<T> solve( const Matrix<T>& B );

    /* ----- Determinant ----- */

    /// Calculate the determinant of the matrix
    /// \return The determinant of the BandedMatrix
    T det();


  };	// End of class BandedMatrix

  /* ----- Constructors ----- */

  template <typename T>
  inline BandedMatrix<T>::BandedMatrix( const std::size_t& n,
    const std::size_t& m1, const std::size_t& m2 ) : N( n ), M1( m1 ), M2( m2 )
  {
    COMPACT.resize( n, m1 + m2 + 1 );
  }

  template <typename T>
  inline BandedMatrix<T>::BandedMatrix( const Matrix<T>& mat,
    const std::size_t& m1, const std::size_t& m2 ) : M1( m1 ), M2( m2 )
  {
    N = mat.rows();
    if ( N != mat.cols() )
    {
      throw Error( "BandedMatrix constructor error: Matrix is not square.");
    }
    COMPACT.resize( N, m1 + m2 + 1 );
    for ( std::size_t i = 0; i < N; i++ ) // Fill the BandedMatrix
    {
      for ( std::size_t j = 0; j < N; j++ )
      {
        if ( i <= j + M1 && j <= i + M2 )
        {
          COMPACT( i, j - i + M1 ) = mat( i, j );
        }
      }
    }
  }

  /* ----- Operator overloading ----- */

  template <class Type>
  inline std::ostream& operator<<( std::ostream& os,
                                   const BandedMatrix<Type>& m )
  {
    os << std::endl;
    if ( m.N == 0 )	{ return os; }
    if ( m.N > 0 )
    {
      for ( std::size_t i = 0; i < m.N; ++i )
      {
        for ( std::size_t j = 0; j < m.N; ++j )
        {
          if ( j > i + m.M2 || i > j + m.M1 )
          {
            os << "\t";
          } else {
            os << "\t" << m.COMPACT( i, j - i + m.M1 );
          }
        }
        os << std::endl;
      }
    }
    return os;
  }

  template <typename T>
  inline const T& BandedMatrix<T>::operator() ( const int& i,
                                                const int& j ) const
  {
    if ( i < 0 )	{ throw Error( "BandedMatrix range error: row < 0." ); }
    if ( j < 0 )	{ throw Error( "BandedMatrix range error: column < 0." ); }
    if ( j > i + M2 || i > j + M1 )
    {
      throw Error( "BandedMatrix range error: index not in a band." );
    }
    return COMPACT( i, j - i + M1 );
  }

  template <typename T>
  inline T& BandedMatrix<T>::operator() ( const int& i, const int& j )
  {
    if ( i < 0 )	{ throw Error( "BandedMatrix range error: row < 0." ); }
    if ( j < 0 )	{ throw Error( "BandedMatrix range error: column < 0." ); }
    if ( j > i + M2 || i > j + M1 )
    {
      throw Error( "BandedMatrix range error: index not in a band." );
    }
    return COMPACT( i, j - i + M1 );
  }

  template <typename T>
  inline BandedMatrix<T>& BandedMatrix<T>::operator=(
                                              const BandedMatrix<T>& original )
  {
    if ( this == &original )
    {
      return *this;
    }
    N = original.N;
    M1 = original.M1;
    M2 = original.M2;
    COMPACT = original.COMPACT;
    return *this;
  }

  template <typename T>
  inline BandedMatrix<T> BandedMatrix<T>::operator+() const
  {
    return *this;
  }

  template <typename T>
  inline BandedMatrix<T> BandedMatrix<T>::operator-() const
  {
    BandedMatrix<T> temp( *this );
    temp.COMPACT = - COMPACT;
    return temp;
  }

  template <typename T>
  inline BandedMatrix<T> BandedMatrix<T>::operator+(
                                          const BandedMatrix<T>& m_plus ) const
  {
    BandedMatrix<T> temp( *this );
    if ( m_plus.N != N || m_plus.M1 != M1 || m_plus.M2 != M2  )
    {
      throw Error( "BandedMatrix dimension error in + operator." );
    }
    temp.COMPACT = COMPACT + m_plus.COMPACT;
    return temp;
  }

  template <typename T>
  inline BandedMatrix<T> BandedMatrix<T>::operator-(
                                          const BandedMatrix<T>& m_minus ) const
  {
    BandedMatrix<T> temp( *this );
    if ( m_minus.N != N || m_minus.M1 != M1 || m_minus.M2 != M2 )
    {
      throw Error( "BandedMatrix dimension error in - operator." );
    }
    temp.COMPACT = COMPACT - m_minus.COMPACT;
    return temp;
  }

  template <typename T>
  inline BandedMatrix<T> BandedMatrix<T>::operator*( const T& scalar ) const
  {
    BandedMatrix<T> temp( *this );
    temp.COMPACT = COMPACT * scalar;
    return temp;
  }

  template <typename T>
  inline BandedMatrix<T> BandedMatrix<T>::operator/( const T& divisor ) const
  {
    BandedMatrix<T> temp( *this );
    temp.COMPACT = COMPACT / divisor;
    return temp;
  }

  template <typename T>
  inline BandedMatrix<T>& BandedMatrix<T>::operator+=(
                                                const BandedMatrix<T>& m_plus )
  {
    if ( m_plus.N != N || m_plus.M1 != M1 || m_plus.M2 != M2 )
    {
      throw Error( "BandedMatrix dimension error in += operator." );
    }
    COMPACT = COMPACT + m_plus.COMPACT;
    return *this;
  }

  template <typename T>
  inline BandedMatrix<T>& BandedMatrix<T>::operator-=(
                                               const BandedMatrix<T>& m_minus )
  {
    if ( m_minus.N != N || m_minus.M1 != M1 || m_minus.M2 != M2 )
    {
      throw Error( "BandedMatrix dimension error in -= operator." );
    }
    COMPACT = COMPACT - m_minus.COMPACT;
    return *this;
  }

  template <typename T>
  inline BandedMatrix<T>& BandedMatrix<T>::operator*=( const T& scalar )
  {
    COMPACT *= scalar;
    return *this;
  }

  template <typename T>
  inline BandedMatrix<T>& BandedMatrix<T>::operator/=( const T& divisor )
  {
    COMPACT /= divisor;
    return *this;
  }

  template <typename T>
  inline BandedMatrix<T>& BandedMatrix<T>::operator+=( const T& add )
  {
    COMPACT += add;
    return *this;
  }

  template <typename T>
  inline BandedMatrix<T>& BandedMatrix<T>::operator-=( const T& minus )
  {
    COMPACT -= minus;
    return *this;
  }

  template <typename T>
  inline Vector<T> BandedMatrix<T>::operator*( Vector<T>& x )
  {
    if ( N != x.size() )
    {
      throw Error( "BandedMatrix: dimensions do not agree in multiply method.");
    }
    Vector<T> result( N );
    int i, j, k, tmploop;
    for ( i = 0; i < N; i++ ) {
  		k = i - M1;
  		tmploop = std::min( M1 + M2 + 1, N - k );
  		result[ i ] = 0.0;
  		for ( j = std::max( 0, - k ); j < tmploop; j++ ){
        result[ i ] += COMPACT( i, j ) * x[ j + k ];
      }
  	}
    return result;
  }

  /* ----- Methods ----- */

  template <typename T>
  inline Matrix<T> BandedMatrix<T>::compact()
  {
    return COMPACT;
  }

  template <typename T>
  inline std::size_t BandedMatrix<T>::size()
  {
    return N;
  }

  template <typename T>
  inline std::size_t BandedMatrix<T>::size_below()
  {
    return M1;
  }

  template <typename T>
  inline std::size_t BandedMatrix<T>::size_above()
  {
    return M2;
  }

  template <typename T>
  inline void BandedMatrix<T>::fill( const T& elem )
  {
    COMPACT.fill( elem );
  }

  /* ----- Solve linear systems ----- */

  template <typename T>
  inline Vector<T> BandedMatrix<T>::solve( const Vector<T>& b )
  {
    if ( b.size() != N ) {
      throw Error( "BandedMatrix solve error: Vector incorrect size." );
    }

    // LU decomposition
    Matrix<T> au, al;
    Vector<int> index;
    double d;
    au = COMPACT;
    al.resize( N, M1 );
    index.resize( N );
    decompose( au, al, index, d );

    // Solve the system
    Vector<T> x( b );
    int i, j, k, l, mm;
    T dum;
    mm = M1 + M2 + 1;
    l = M1;
    for ( k = 0; k < N; k++ ) {
      j = index[ k ] - 1;
      if ( j != k ) std::swap( x[ k ], x[ j ] );
      if ( l < N ) l++;
      for ( j = k + 1; j < l; j++ ){
        x[ j ] -= al( k, j - k - 1 ) * x[ k ];
      }
    }
    l = 1;
    for ( i = N - 1; i >= 0; i-- ) {
      dum = x[ i ];
      for ( k = 1; k < l; k++ ) {
        dum -= au( i, k ) * x[ k + i ];
      }
      x[ i ] = dum / au( i, 0 );
      if ( l < mm ) l++;
    }
    return x;
  }

  template <typename T>
  inline Matrix<T> BandedMatrix<T>::solve( const Matrix<T>& B )
  {
    if ( B.rows() != N ) {
      throw Error( "BandedMatrix solve error: Matrix rows incorrect size." );
    }
    Matrix<T> X( B );

    // LU decomposition
    Matrix<T> au, al;
    Vector<int> index;
    double d;
    au = COMPACT;
    al.resize( N, M1 );
    index.resize( N );
    decompose( au, al, index, d );

    Vector<T> x( N );
    int i, j, k, l, mm;
    T dum;
    mm = M1 + M2 + 1;

    for ( std::size_t col = 0; col < B.cols(); ++col )
    {
      // Solve the system for each column
      x = X.get_col( col );
      l = M1;
      for ( k = 0; k < N; k++ ) {
        j = index[ k ] - 1;
        if ( j != k ) std::swap( x[ k ], x[ j ] );
        if ( l < N ) l++;
        for ( j = k + 1; j < l; j++ ){
          x[ j ] -= al( k, j - k - 1 ) * x[ k ];
        }
      }
      l = 1;
      for ( i = N - 1; i >= 0; i-- ) {
        dum = x[ i ];
        for ( k = 1; k < l; k++ ) {
          dum -= au( i, k ) * x[ k + i ];
        }
        x[ i ] = dum / au( i, 0 );
        if ( l < mm ) l++;
      }
      X.set_col( col, x );
    }
    return X;
  }

  /* ----- Determinant ----- */

  template <typename T>
  inline T BandedMatrix<T>::det()
  {
    Matrix<T> au, al;
    Vector<int> index;
    double d;
    au = COMPACT;
    al.resize( N, M1 );
    index.resize( N );
    decompose( au, al, index, d );
    T dd = d;
    for ( int i = 0; i < N; i++ ) {
      dd *= au( i, 0 );
    }
    return dd;
  }

  /* ----- Private ----- */

  template <typename T>
  inline void BandedMatrix<T>::decompose( Matrix<T>& au, Matrix<T>& al,
                                          Vector<int>& index, double& d )
  {
    const double TINY( 1.0e-40 );
    std::size_t i, j, k, l, mm;
    T dum;
    mm = M1 + M2 + 1;
    l = M1;
    for ( i = 0; i < M1; i++ ) {               // Rearrange the storage a bit
      for ( j = M1 - i; j< mm ; j++ ){
        au( i, j - l ) = au( i, j );
      }
      l--;
      for ( j = mm - l - 1; j < mm; j++ ){
        au( i, j ) = 0.0;
      }
    }
    d = 1.0;
    l = M1;
    for ( k = 0; k < N; k++ ) {
      dum = au( k, 0 );
      i = k;
      if ( l < N ){ l++; }
      for ( j = k + 1; j < l; j++ ) {          // Find the pivot element
        if ( std::abs( au( j, 0 ) ) > std::abs( dum ) ) {
          dum = au( j, 0 );
          i = j;
        }
      }
      index[ k ] = i + 1;
      if ( std::abs( dum ) == 0.0 ) au( k, 0 ) = TINY;
      // Singular matrix but proceed anyway with TINY pivot
      if (i != k) {
        d = -d;
        for ( j = 0; j < mm; j++ ){
          std::swap( au( k, j ), au( i, j ) );
        }
      }
      for ( i = k + 1; i < l; i++ ) {
        dum = au( i, 0 ) / au( k, 0 );
        al( k, i - k - 1 ) = dum;
        for ( j = 1; j < mm; j++ ){
          au( i, j - 1 ) = au( i, j ) - dum * au( k, j );
        }
        au( i, mm - 1 ) = 0.0;
      }
    }
  }

}  // End of namespace Luna

#endif
