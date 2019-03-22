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

    /// Constructor
    //TODO specified size but empty (need access operators)

    /// Constructor for a BandedMatrix from a dense Matrix
    /// \param mat The dense Matrix used to create the BandedMatrix
    /// \param m1 The number of bands below the main diagonal
    /// \param m2 The number of bands above the main diagonal
    BandedMatrix( const Matrix<T>& mat, const std::size_t& m1,
                  const std::size_t& m2 );

    /// Destructor
    ~BandedMatrix(){}

    /* ----- Operator overloading ----- */

    //TODO print the BandedMatrix to the screen -> convert COMPACT to nice output
    // leaving zero entries outside bands blank

    //TODO access operator, =, +, - (unary and binary), * (scalar and Vector),
    // /, +=, -=, *= (scalar), /=,

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

    /* ----- Solve linear systems ----- */

    /// Solve the banded system of equations Ax=b where x and b are Vectors
    /// \param b The right-side Vector of the system of equations
    /// \return The solution Vector x
    Vector<T> solve( const Vector<T>& b );

    //TODO solve with Matrix input -> just decompose once

    /* ----- Determinant ----- */

    /// Calculate the determinant of the matrix
    /// \return The determinant of the BandedMatrix
    T det();


  };	// End of class BandedMatrix

  /* ----- Constructors ----- */

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
    /*
    // Main diagonal
    for ( std::size_t i = 0; i < N; i++ )
    {
      COMPACT( i, m1 ) = mat( i, i );
    }
    // Below
    for ( std::size_t l = 1; l <= m1; l++ )
    {
      for ( std::size_t i = 0; i < N - l; i++ )
      {
        COMPACT( i + l, m1 - l ) = mat( i + l, i );
      }
    }
    // Above
    for ( std::size_t l = 1; l <= m2; l++ )
    {
      for ( std::size_t i = 0; i < N - l; i++ )
      {
        COMPACT( i, m1 + l ) = mat( i, i + l );
      }
    }
    */

    for ( std::size_t i = 0; i < N; i++ )
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
