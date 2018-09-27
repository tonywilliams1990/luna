/// \file Matrix.h
/// A templated matrix class which constructs a matrix as an std::vector of
/// Vectors.

#ifndef MATRIX_H
#define MATRIX_H

#include <vector>

#include "Error.h"
#include "Vector.h"

namespace Luna
{

  // Templated matrix class
  template <class T>

  /// A Matrix class for use with double and std::complex<double>
  class Matrix
  {
  private:
    std::vector< Vector<T> > MATRIX;    // Matrix as std::vector of Vectors
    std::size_t ROWS;					          // Number of rows
    std::size_t COLS;					          // Number of columns

  public:

    /// Constructor for an empty Matrix of unspecified size
    Matrix() : ROWS( 0 ), COLS( 0 ) {}

    /// Constructor for a matrix with specfied initial elements
    /// \param rows The number of rows in the matrix.
    /// \param cols The number of columns in the matrix.
    /// \param elem The entry to be placed in all elements.
    Matrix( const std::size_t& rows, const std::size_t& cols, const T& elem );

    /// Copy constructor
    /// \param source The source Matrix to be copied
    Matrix( const Matrix<T>& source );

    /// Destructor
    ~Matrix(){}

    /* ----- Operator overloading ----- */

    /// Output operator <<
    template <class Type>
    friend std::ostream& operator<<( std::ostream& os, const Matrix<Type>& m );

    /// Indexing operator ( read only )
    /// \param i Row index
    /// \param j Column index
    const T& operator() ( const std::size_t& i, const std::size_t& j ) const;

    /// Indexing operator ( read / write )
    /// \param i Row index
    /// \param j Column index
    T& operator() ( const std::size_t& i, const std::size_t& j );

    /// Copy assignment
    /// \param original Matrix to be copied
    /// \return The new Matrix to which the original matrix is assigned
    Matrix<T>& operator=( const Matrix<T>& original );

    /// Unary +
    /// \return The Matrix
    Matrix<T> operator+() const;

    /// Unary -
    /// \return The negation of the Matrix
    Matrix<T> operator-() const;

    /// Binary +
    /// \param m_plus The Matrix to be added
    /// \return The sum of this Matrix and m_plus
    Matrix<T> operator+( const Matrix<T>& m_plus ) const;

    /// Binary -
    /// \param m_minus The Matrix to be subtracted
    /// \return The subtraction of m_minus from this Matrix
    Matrix<T> operator-( const Matrix<T>& m_minus ) const;

    /// Scalar multiplication
    Matrix<T> operator*( const T& scalar ) const;

    /// Scalar division

    /// Addition assignment

    /// Subtraction assignment

    /// Scalar multiplication assignment

    /// Scalar division assigment

    /// Addition assignment ( add a constant to all elements )

    /// Matrix multiplication

    /// Matrix Vector multiplication

    /// Operator overloading for row access



    /* ----- Methods ----- */

    /// Return the number of rows in the Matrix
    /// \return The number of rows in the Matrix
    std::size_t rows() const;

    /// Return the number of columns in the Matrix
    /// \return The number of columns in the Matrix
    std::size_t cols() const;

    /// Return the number of elements in the Matrix
    /// \return The number of elements in the Matrix
    std::size_t numel() const;


    //TODO add method - add in place without copy (same for subtract)


  };	// End of class Matrix

  template <typename T>
  inline Matrix<T>::Matrix( const std::size_t& rows, const std::size_t& cols,
                            const T& elem ) : ROWS( rows ), COLS( cols )
  {
    const Vector<T> row( cols, elem );
    MATRIX.reserve( rows );
    for ( std::size_t i=0; i < rows; ++i )
    {
      MATRIX.push_back( row );
    }
  }

  template <typename T>
  inline Matrix<T>::Matrix( const Matrix<T>& source )
  {
    *this = source;
  }

  /* ----- Operator overloading ----- */

  template <class Type>
  inline std::ostream& operator<<( std::ostream& os, const Matrix<Type>& m )
  {
    os << std::endl;
    for (std::size_t j = 0; j < m.ROWS; ++j )
    {
      os << "\t" << m.MATRIX[ j ] << std::endl;
    }
    return os;
  }

  template <typename T>
  inline const T& Matrix<T>::operator()( const std::size_t& i,
                                         const std::size_t& j ) const
  {
    if ( i<0 || ROWS<=i )	{ throw Error( "Matrix range error: dimension 1" );}
    if ( j<0 || COLS<=j )	{ throw Error( "Matrix range error: dimension 2" );}
    return MATRIX[ i ][ j ];
  }

  template <typename T>
  inline T& Matrix<T>::operator()( const std::size_t& i, const std::size_t& j )
  {
    if ( i<0 || ROWS<=i )	{ throw Error( "Matrix range error: dimension 1" );}
    if ( j<0 || COLS<=j )	{ throw Error( "Matrix range error: dimension 2" );}
    return MATRIX[ i ][ j ];
  }

  template <typename T>
  inline Matrix<T>& Matrix<T>::operator=( const Matrix<T>& original )
  {
    if ( this == &original )
    {
      return * this;
    }
    MATRIX = original.MATRIX;
    ROWS = original.ROWS;
    COLS = original.COLS;
    return *this;
  }

  template <typename T>
  inline Matrix<T> Matrix<T>::operator+() const
  {
    return *this;
  }

  template <typename T>
  inline Matrix<T> Matrix<T>::operator-() const
  {
    Matrix<T> temp( *this );
    std::transform ( temp.MATRIX.cbegin(), temp.MATRIX.cend(),
                     temp.MATRIX.begin(), std::negate< Vector<T> >() );
    return temp;
  }

  template <typename T>
  inline Matrix<T> Matrix<T>::operator+( const Matrix<T>& m_plus ) const
  {
    Matrix<T> temp( *this );
    if ( m_plus.rows() != ROWS || m_plus.cols() != COLS )
    {
      throw Error( "Matrix dimension error in + operator." );
    }
    std::transform ( temp.MATRIX.cbegin(), temp.MATRIX.cend(),
                     m_plus.MATRIX.begin(), temp.MATRIX.begin(),
                     std::plus< Vector<T> >() );
    return temp;
  }

  template <typename T>
  inline Matrix<T> Matrix<T>::operator-( const Matrix<T>& m_minus ) const
  {
    Matrix<T> temp( *this );
    if ( m_minus.rows() != ROWS || m_minus.cols() != COLS )
    {
      throw Error( "Matrix dimension error in + operator." );
    }
    std::transform ( temp.MATRIX.cbegin(), temp.MATRIX.cend(),
                     m_minus.MATRIX.begin(), temp.MATRIX.begin(),
                     std::minus< Vector<T> >() );
    return temp;
  }

  template <typename T>
  struct scale_vector
  {
    scale_vector( T value )
    {
      SCALAR = value;
    }

    Vector<T> operator() ( Vector<T> vec )
    {
      return vec * SCALAR;
    }

    T SCALAR;
  };

  template <typename T>
  inline Matrix<T> Matrix<T>::operator*( const T& scalar ) const
  {
    Matrix<T> temp( *this );
    std::transform ( temp.MATRIX.begin(), temp.MATRIX.end(),
                     temp.MATRIX.begin(), scale_vector<T>( scalar ) );
    return temp;
  }


  /* ----- Methods ----- */

  template <typename T>
  inline std::size_t Matrix<T>::rows() const
  {
    return ROWS;
  }

  template <typename T>
  inline std::size_t Matrix<T>::cols() const
  {
    return COLS;
  }

  template <typename T>
  inline std::size_t Matrix<T>::numel() const
  {
    return ROWS * COLS;
  }



}  // End of namespace Luna

#endif
