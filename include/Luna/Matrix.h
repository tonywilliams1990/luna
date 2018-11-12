/// \file Matrix.h
/// A templated matrix class which constructs a matrix as an std::vector of
/// Vectors.

#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <random>
#include <chrono>

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
    /// \param scalar The scalar to multiply the Matrix by
    /// \return The Matrix multiplied by the scalar
    Matrix<T> operator*( const T& scalar ) const;
    friend Matrix<T> operator*( const T& scalar, Matrix<T>& mat )
    {
      return mat * scalar;
    }

    /// Scalar division
    /// \param divisor The divisor to divide the Matrix by
    /// \return The Matrix divided by the divisor
    Matrix<T> operator/( const T& divisor ) const;

    /// Addition assignment
    /// \param m_plus The Matrix to be added
    /// \return A reference to the Matrix after addition by m_plus
    Matrix<T>& operator+=( const Matrix<T>& m_plus );

    /// Subtraction assignment
    /// \param m_minus The Matrix to be subtracted
    /// \return A reference to the Matrix after subtraction by m_minus
    Matrix<T>& operator-=( const Matrix<T>& m_minus );

    /// Scalar multiplication assignment
    /// \param scalar The scalar to multiply the Matrix by
    /// \return A reference to the Matrix after multiplication by a scalar
    Matrix<T>& operator*=( const T& scalar );

    /// Scalar division assigment
    /// \param divisor The divisor to divide the Matrix by
    /// \return A reference to the Matrix after division by a scalar
    Matrix<T>& operator/=( const T& divisor );

    /// Constant addition assignment
    /// \param add The constant to be added to each element
    /// \return A reference to the Matrix after addition by a constant
    Matrix<T>& operator+=( const T& add );

    /// Constant subtraction assignment
    /// \param minus The constant to be subtracted from each element
    /// \return A reference to the Matrix after subtraction by a constant
    Matrix<T>& operator-=( const T& minus );

    /// Operator overloading for row access ( read only )
    /// \param i The row index
    /// \return A constant reference to the Vector located at the given row
    const Vector<T>& operator[] ( const std::size_t& i ) const;

    /// Operator overloading for row access ( read / write )
    /// \param i The row index
    /// \return A reference to the Vector located at the given row
    Vector<T>& operator[] ( const std::size_t& i );

    /// Matrix multiplication
    /// \param B The matrix which is to be multiplied
    /// \return The result matrix A * B
    Matrix<T> operator*( Matrix<T>& B );

    /// Matrix Vector multiplication
    /// \param x The vector which is to be multiplied
    /// \return The result vector A * x
    Vector<T> operator*( Vector<T>& x );


    /* ----- Methods ----- */

    /// Matrix multiplication
    /// \param B The matrix which is to be multiplied
    /// \return The result matrix A * B
    Matrix<T> multiply( Matrix<T>& B );

    /// Matrix Vector multiplication
    /// \param x The vector which is to be multiplied
    /// \return The result vector A * x
    Vector<T> multiply( const Vector<T>& x );

    /// TODO parallel multiplication in separate function?

    /// Return the number of rows in the Matrix
    /// \return The number of rows in the Matrix
    std::size_t rows() const;

    /// Return the number of columns in the Matrix
    /// \return The number of columns in the Matrix
    std::size_t cols() const;

    /// Return the number of elements in the Matrix
    /// \return The number of elements in the Matrix
    std::size_t numel() const;

    /// Set a row of the matrix using a Vector
    /// \param i The row index
    /// \param vec The Vector used to replace the specified row
    void set_row( const std::size_t& row, Vector<T>& vec );

    /// Set a column of the matrix using a Vector
    /// \param col The column index
    /// \param vec The Vector used to replace the specified column
    void set_col( const std::size_t& col, const Vector<T>& vec );

    //TODO get_row ( wrap [] operator )

    /// Get a column of the matrix as a Vector
    /// \param col The column index
    /// \return A Vector of the column of the Matrix
    Vector<T> get_col( const std::size_t& col );

    /// Return the transpose of the Matrix
    /// \return The transpose of the Matrix
    Matrix<T> transpose();

    /// Replace the current matrix with its transpose
    void transpose_in_place();

    //TODO adjoint ???

    //TODO conjugate

    //TODO conjugate transpose

    //TODO resize

    //TODO swap_rows

    //TODO fill, fill_diag, fill_band, fill_tridiag

    /// Fill the matrix with random elements (between 0 and 1)
    void random();

    /* ----- Norms ----- */

    /* ----- Determinant ----- */

    /* ---- Solve linear systems ----- */


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
      if ( m.MATRIX[ j ].size() != m.COLS ){
        throw Error( "Matrix dimension error: wrong sized row in Matrix." );
      }
      os << "\t" << m.MATRIX[ j ] << std::endl;
    }
    return os;
  }

  template <typename T>
  inline const T& Matrix<T>::operator()( const std::size_t& i,
                                         const std::size_t& j ) const
  {
    if ( i<0 || ROWS<=i )	{ throw Error( "Matrix range error: dimension 1." );}
    if ( j<0 || COLS<=j )	{ throw Error( "Matrix range error: dimension 2." );}
    return MATRIX[ i ][ j ];
  }

  template <typename T>
  inline T& Matrix<T>::operator()( const std::size_t& i, const std::size_t& j )
  {
    if ( i<0 || ROWS<=i )	{ throw Error( "Matrix range error: dimension 1." );}
    if ( j<0 || COLS<=j )	{ throw Error( "Matrix range error: dimension 2." );}
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
      throw Error( "Matrix dimension error in - operator." );
    }
    std::transform ( temp.MATRIX.cbegin(), temp.MATRIX.cend(),
                     m_minus.MATRIX.begin(), temp.MATRIX.begin(),
                     std::minus< Vector<T> >() );
    return temp;
  }

  template <typename T>
  struct scale_vector
  {
    scale_vector( const T& scalar )
    {
      SCALAR = scalar;
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

  template <typename T>
  struct div_vector
  {
    div_vector( const T& divisor )
    {
      DIV = divisor;
    }

    Vector<T> operator() ( Vector<T> vec )
    {
      return vec / DIV;
    }

    T DIV;
  };

  template <typename T>
  inline Matrix<T> Matrix<T>::operator/( const T& divisor ) const
  {
    Matrix<T> temp( *this );
    std::transform ( temp.MATRIX.begin(), temp.MATRIX.end(),
                     temp.MATRIX.begin(), div_vector<T>( divisor ) );
    return temp;

  }

  template <typename T>
  inline Matrix<T>& Matrix<T>::operator+=( const Matrix<T>& m_plus )
  {
    if ( m_plus.rows() != ROWS || m_plus.cols() != COLS )
    {
      throw Error( "Matrix dimension error in += operator." );
    }
    std::transform ( MATRIX.begin(), MATRIX.end(),
                     m_plus.MATRIX.begin(), MATRIX.begin(),
                     std::plus< Vector<T> >() );
    return *this;
  }

  template <typename T>
  inline Matrix<T>& Matrix<T>::operator-=( const Matrix<T>& m_minus )
  {
    if ( m_minus.rows() != ROWS || m_minus.cols() != COLS )
    {
      throw Error( "Matrix dimension error in -= operator." );
    }
    std::transform ( MATRIX.begin(), MATRIX.end(),
                     m_minus.MATRIX.begin(), MATRIX.begin(),
                     std::minus< Vector<T> >() );
    return *this;
  }

  template <typename T>
  inline Matrix<T>& Matrix<T>::operator*=( const T& scalar )
  {
    std::transform ( MATRIX.begin(), MATRIX.end(), MATRIX.begin(),
                     scale_vector<T>( scalar ) );
    return *this;
  }

  template <typename T>
  inline Matrix<T>& Matrix<T>::operator/=( const T& divisor )
  {
    std::transform ( MATRIX.begin(), MATRIX.end(), MATRIX.begin(),
                     div_vector<T>( divisor ) );
    return *this;
  }

  template <typename T>
  struct add_vector
  {
    add_vector( const T& add )
    {
      ADD = add;
    }

    Vector<T> operator() ( Vector<T> vec )
    {
      return vec += ADD;
    }

    T ADD;
  };

  template <typename T>
  inline Matrix<T>& Matrix<T>::operator+=( const T& add )
  {
    std::transform ( MATRIX.begin(), MATRIX.end(),  MATRIX.begin(),
                     add_vector<T>( add ) );
    return *this;
  }

  template <typename T>
  struct minus_vector
  {
    minus_vector( const T& minus )
    {
      MINUS = minus;
    }

    Vector<T> operator() ( Vector<T> vec )
    {
      return vec -= MINUS;
    }

    T MINUS;
  };

  template <typename T>
  inline Matrix<T>& Matrix<T>::operator-=( const T& minus )
  {
    std::transform ( MATRIX.begin(), MATRIX.end(),  MATRIX.begin(),
                     minus_vector<T>( minus ) );
    return *this;
  }

  template <typename T>
  inline const Vector<T>& Matrix<T>::operator[]( const std::size_t& i ) const
  {
    if ( i < 0 || ROWS <= i )
    {
      throw Error( "Matrix operator [] range error" );
    }
    return MATRIX[ i ];
  }

  template <typename T>
  inline Vector<T>& Matrix<T>::operator[]( const std::size_t& i )
  {
    if ( i < 0 || ROWS <= i )
    {
      throw Error( "Matrix operator [] range error" );
    }
    return MATRIX[ i ];
  }

  template <typename T>
  inline Matrix<T> Matrix<T>::operator*( Matrix<T>& B )
  {
    return multiply( B );
  }

  template <typename T>
  inline Vector<T> Matrix<T>::operator*( Vector<T>& x )
  {
    return multiply( x );
  }

  /* ----- Methods ----- */

  template <typename T>
  inline Matrix<T> Matrix<T>::multiply( Matrix<T>& B )
  {
    if ( B.rows() != COLS )
    {
      throw Error( "Matrix dimensions do not agree in multiply method." );
    }
    Matrix<T> temp( ROWS, B.cols(), 0.0 );
    for ( std::size_t col = 0; col < B.cols(); ++col )
    {
      temp.set_col( col, multiply( B.get_col( col ) ) );
    }
    return temp;
  }

  template <typename T>
  inline Vector<T> Matrix<T>::multiply( const Vector<T>& x )
  {
    if ( x.size() != COLS )
    {
      throw Error( "Dimensions do not agree in Matrix multiply method." );
    }
    Vector<T> temp;
    temp.reserve( ROWS );
    for ( std::size_t row = 0; row < ROWS; ++row )
    {
      temp.push_back( MATRIX[ row ].dot( x ) );
    }
    return temp;
  }

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

  template <typename T>
  inline void Matrix<T>::set_row( const std::size_t& row, Vector<T>& vec )
  {
    if ( vec.size() != COLS ){
      throw Error( "Size error in set_row method." );
    }
    if ( row < 0 || ROWS <= row ){
      throw Error( "Range error in set_row method." );
    }
    MATRIX[ row ] = vec;
  }

  template <typename T>
  inline void Matrix<T>::set_col( const std::size_t& col, const Vector<T>& vec )
  {
    for ( std::size_t row = 0; row < ROWS; ++row )
    {
      MATRIX[ row ][ col ] = vec[ row ];
    }
  }

  template <typename T>
  inline Vector<T> Matrix<T>::get_col( const std::size_t& col )
  {
    Vector<T> temp( ROWS, 0.0 );
    for ( std::size_t row = 0; row < ROWS; ++row )
    {
      temp[ row ] = MATRIX[ row ][ col ];
    }
    return temp;
  }


  template <typename T>
  inline Matrix<T> Matrix<T>::transpose()
  {
    Matrix<T> temp( *this );
    temp.transpose_in_place();
    return temp;
  }

  template <typename T>
  inline void Matrix<T>::transpose_in_place()
  {
    if( ROWS == COLS ) // Square matrix
    {
      for ( std::size_t i = 0; i < ROWS; ++i )
      {
        for ( std::size_t j = i + 1; j < COLS; ++j )
        {
          std::swap( MATRIX[ i ][ j ], MATRIX[ j ][ i ] );
        }
      }
    } else // Non-square matrix
    {
      std::vector< Vector<T> > temp;
      temp.resize( COLS );
      for ( std::size_t row = 0; row < COLS; ++row )
       {
         temp[ row ].resize( ROWS );
       }
       for ( std::size_t i = 0; i < ROWS; ++i )
       {
         for ( std::size_t j = 0; j < COLS; ++j )
         {
           temp[ j ][ i ] = MATRIX[ i ][ j ];
         }
       }
       MATRIX = temp;
       std::swap( ROWS, COLS );
    }
  }


  template <>
  inline void Matrix<double>::random()
  {
    std::mt19937_64 rng;
    // initialize the random number generator with time-dependent seed
    auto start = std::chrono::high_resolution_clock::now();
    uint64_t timeSeed = start.time_since_epoch().count();
    std::seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed>>32)};
    rng.seed(ss);
    // initialize a uniform distribution between 0 and 1
    std::uniform_real_distribution<double> unif(0, 1);
    for ( std::size_t i = 0; i < ROWS; ++i )
    {
      for ( std::size_t j = 0; j < COLS; ++j )
      {
        double num( unif(rng) );
        MATRIX[ i ][ j ] = num;
      }
    }
  }

  template <>
  inline void Matrix< std::complex<double> >::random()
  {
    std::mt19937_64 rng;
    // initialize the random number generator with time-dependent seed
    auto start = std::chrono::high_resolution_clock::now();
    uint64_t timeSeed = start.time_since_epoch().count();
    std::seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed>>32)};
    rng.seed(ss);
    // initialize a uniform distribution between 0 and 1
    std::uniform_real_distribution<double> unif(0, 1);
    for ( std::size_t i = 0; i < ROWS; ++i )
    {
      for ( std::size_t j = 0; j < COLS; ++j )
      {
        double real( unif(rng) );
        double imag( unif(rng) );
        MATRIX[ i ][ j ] = std::complex<double> ( real, imag );
      }
    }
  }

}  // End of namespace Luna

#endif
