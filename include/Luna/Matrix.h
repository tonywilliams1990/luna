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

    template <typename Type>
    struct scale_vector
    {
      scale_vector( const Type& scalar )
      {
        SCALAR = scalar;
      }

      Vector<Type> operator() ( Vector<Type> vec )
      {
        return vec * SCALAR;
      }

      Type SCALAR;
    };

    template <typename Type>
    struct div_vector
    {
      div_vector( const Type& divisor )
      {
        DIV = divisor;
      }

      Vector<Type> operator() ( Vector<Type> vec )
      {
        return vec / DIV;
      }

      Type DIV;
    };

    template <typename Type>
    struct add_vector
    {
      add_vector( const Type& add )
      {
        ADD = add;
      }

      Vector<Type> operator() ( Vector<Type> vec )
      {
        return vec += ADD;
      }

      Type ADD;
    };

    template <typename Type>
    struct minus_vector
    {
      minus_vector( const Type& minus )
      {
        MINUS = minus;
      }

      Vector<Type> operator() ( Vector<Type> vec )
      {
        return vec -= MINUS;
      }

      Type MINUS;
    };

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
    /// \param B The Matrix which is to be multiplied
    /// \return The result Matrix A * B
    Matrix<T> operator*( Matrix<T>& B );

    /// Matrix Vector multiplication
    /// \param x The vector which is to be multiplied
    /// \return The result vector A * x
    Vector<T> operator*( Vector<T>& x );


    /* ----- Methods ----- */

    /// Matrix multiplication
    /// \param B The Matrix which is to be multiplied
    /// \return The result Matrix A * B
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

    /// Set a row of the Matrix using a Vector
    /// \param row The row index
    /// \param vec The Vector used to replace the specified row
    void set_row( const std::size_t& row, Vector<T>& vec );

    /// Set a column of the Matrix using a Vector
    /// \param col The column index
    /// \param vec The Vector used to replace the specified column
    void set_col( const std::size_t& col, const Vector<T>& vec );

    /// Get a row of the Matrix as a Vector ( wrap [] operator )
    /// \param row The row index
    /// \return A Vector of the row of the Matrix
    Vector<T> get_row( const std::size_t& row );

    /// Get a column of the Matrix as a Vector
    /// \param col The column index
    /// \return A Vector of the column of the Matrix
    Vector<T> get_col( const std::size_t& col );

    /// Return the transpose of the Matrix
    /// \return The transpose of the Matrix
    Matrix<T> transpose() const;

    /// Replace the current Matrix with its transpose
    void transpose_in_place();

    /// Take the complex conjugate of each element of the Matrix
    /// \return The complex conjugate Matrix
    Matrix<T> conjugate() const;

    /// Transpose of the Matrix along with the complex conjugate of each entry
    /// \return The conjugate transpose Matrix
    Matrix<T> conjugate_transpose() const;

    /// Resize the Matrix while (empty entries are appended if necessary)
    /// \param rows The number of rows in the Matrix
    /// \param cols The number of columns in the Matrix
    void resize( const std::size_t& rows, const std::size_t& cols );

    /// Swap two rows of the Matrix
    /// \param i The index of the first row
    /// \param j The index of the second row
    void swap_rows( const std::size_t& i, const std::size_t& j );

    /// Fill the Matrix with specified elements
    /// \param elem The element to fill the Matrix with
    void fill( const T& elem );

    /// Fill the leading diagonal of the Matrix with specified elements
    /// \param elem The element to fill the leading diagonal with
    void fill_diag( const T& elem );

    /// Fill a diagonal band of the Matrix
    /// \param offset The offset from the leading diagonal (+ above, - below)
    /// \param elem The element to fill the diagonal band with
    void fill_band( const std::size_t& offset, const T& elem );

    /// Fill the main three diagonals of the Matrix
    /// \param lower The element to fill the lower band
    /// \param diag The element to fill the main diagonal
    /// \param upper The element to fill the upper band
    void fill_tridiag( const T& lower, const T& diag, const T& upper );

    /// Fill the Matrix with random elements (between -1 and 1)
    void random();

    /* ----- Norms ----- */

    /// The Matrix one-norm
    /// \return The maximum absolute column sum of the Matrix
    double norm_1() const;

    /// The Matrix infinity-norm
    /// \return The maximum absolute row sum of the Matrix
    double norm_inf() const;

    /// The Matrix p-norm (p=2 is Frobenius, p=inf is max norm)
    /// \param p The exponent for the p-norm
    /// \return The entrywise p-norm of the Matrix
    double norm_p( const double& p ) const;

    /// The Matrix Frobenius norm
    /// \return The Frobenius norm of the Matrix
    double norm_frob() const;

    /// The Matrix max-norm
    /// \return The entrywise max-norm of the Matrix
    double norm_max() const;

    /* ----- Determinant ----- */

    //TODO

    /* ---- Solve linear systems ----- */

    //TODO solve_basic, solve_parallel, solve -> string to choose method
    //TODO Vector and Matrix versions of each 


  };	// End of class Matrix

  template <typename T>
  inline Matrix<T>::Matrix( const std::size_t& rows, const std::size_t& cols,
                            const T& elem ) : ROWS( rows ), COLS( cols )
  {
    const Vector<T> row( cols, elem );
    MATRIX.reserve( rows );
    for ( std::size_t i = 0; i < rows; ++i )
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
  inline Matrix<T> Matrix<T>::operator*( const T& scalar ) const
  {
    Matrix<T> temp( *this );
    std::transform ( temp.MATRIX.begin(), temp.MATRIX.end(),
                     temp.MATRIX.begin(), scale_vector<T>( scalar ) );
    return temp;
  }

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
  inline Matrix<T>& Matrix<T>::operator+=( const T& add )
  {
    std::transform ( MATRIX.begin(), MATRIX.end(),  MATRIX.begin(),
                     add_vector<T>( add ) );
    return *this;
  }

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
  inline Vector<T> Matrix<T>::get_row( const std::size_t& row )
  {
    Vector<T> temp;
    temp = MATRIX[ row ];
    return temp;
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
  inline Matrix<T> Matrix<T>::transpose() const
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

  template <typename T>
  inline Matrix<T> Matrix<T>::conjugate() const
  {
     Matrix<T> temp( *this );
     for ( std::size_t i = 0; i < ROWS; ++i )
     {
       temp.MATRIX[ i ] = temp.MATRIX[ i ].conjugate();
     }
     return temp;
  }

  template <typename T>
  inline Matrix<T> Matrix<T>::conjugate_transpose() const
  {
    Matrix<T> temp( *this );
    temp.transpose_in_place();
    return temp.conjugate();
  }

  template <typename T>
  inline void Matrix<T>::resize( const std::size_t& rows,
                                 const std::size_t& cols )
  {
    Matrix<T> temp( *this );
    MATRIX.clear();
    // Create an empty Matrix
    const Vector<T> row( cols, 0 );
    MATRIX.reserve( rows );
    for ( std::size_t i = 0; i < rows; ++i )
    {
      MATRIX.push_back( row );
    }
    // Refill the Matrix using old values
    for ( std::size_t i = 0; i < rows; ++i )
    {
      for ( std::size_t j = 0; j < cols; ++j )
      {
        if ( i<ROWS && j<COLS )
        {
          MATRIX[ i ][ j ] = temp.MATRIX[ i ][ j ];
        }
      }
    }
    ROWS = rows;
    COLS = cols;
  }

  template <typename T>
  inline void Matrix<T>::swap_rows( const std::size_t& i, const std::size_t& j )
  {
    if ( i < 0 || ROWS <= i )
    {
      throw Error( "Matrix swap row range error in first index" );
    }
    if ( j < 0 || ROWS <= j )
    {
      throw Error( "Matrix swap row range error in second index" );
    }
    std::swap< Vector<T> > ( MATRIX[ i ], MATRIX[ j ] );
  }

  template <typename T>
  inline void Matrix<T>::fill( const T& elem )
  {
    for ( std::size_t i = 0; i < ROWS; ++i )
    {
      for ( std::size_t j = 0; j < COLS; ++j )
      {
        MATRIX[ i ][ j ] = elem;
      }
    }
  }

  template <typename T>
  inline void Matrix<T>::fill_diag( const T& elem )
  {
    std::size_t N( ROWS );
    if ( COLS < ROWS ) { N = COLS; }
    for ( std::size_t i = 0; i < N ; ++i )
    {
      MATRIX[ i ][ i ] = elem;
    }
  }

  template <typename T>
  inline void Matrix<T>::fill_band( const std::size_t& offset, const T& elem )
  {
    for ( std::size_t row = 0; row < ROWS; ++row )
    {
        if ( ( row + offset < COLS ) && ( row + offset >= 0 ) )
        {
            MATRIX[ row ][ row + offset ] = elem;
        }
    }
  }

  template <typename T>
  inline void Matrix<T>::fill_tridiag( const T& lower, const T& diag,
                                       const T& upper )
  {
    fill_band( -1, lower );
    fill_diag( diag );
    fill_band( 1, upper );
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
    // initialize a uniform distribution between -1 and 1
    std::uniform_real_distribution<double> unif(-1, 1);
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
    // initialize a uniform distribution between -1 and 1
    std::uniform_real_distribution<double> unif(-1, 1);
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

  /* ----- Norms ----- */

  template <typename T>
  inline double Matrix<T>::norm_1() const
  {
    double max( 0.0 );
    for ( std::size_t j = 0; j < COLS; ++j )
    {
        double sum( 0.0 );
        for ( std::size_t i = 0; i < ROWS; ++i )
        {
            sum += std::abs( MATRIX[ i ][ j ] );
        }
        max = std::max( max, sum );
    }
    return max;
  }

  template <typename T>
  inline double Matrix<T>::norm_inf() const
  {
    double max( 0.0 );
    for ( std::size_t i = 0; i < ROWS; ++i )
    {
        double sum( 0.0 );
        for ( std::size_t j = 0; j < COLS; ++j )
        {
            sum += std::abs( MATRIX[ i ][ j ] );
        }
        max = std::max( max, sum );
    }
    return max;
  }

  template <typename T>
  inline double Matrix<T>::norm_p( const double& p ) const
  {
    double sum( 0.0 );
    for ( std::size_t i = 0; i < ROWS; ++i )
    {
        for ( std::size_t j = 0; j < COLS; ++j )
        {
            sum += std::pow( std::abs( MATRIX[ i ][ j ] ), p );
        }
    }
    return std::pow( sum , 1.0/p );
  }

  template <typename T>
  inline double Matrix<T>::norm_frob() const
  {
    return this->norm_p( 2.0 );
  }

  template <typename T>
  inline double Matrix<T>::norm_max() const
  {
    double max( 0.0 );
    for ( std::size_t i = 0; i < ROWS; ++i )
    {
        for ( std::size_t j = 0; j < COLS; ++j )
        {
            max = std::max( max, std::abs( MATRIX[ i ][ j ] ) );
        }
    }
    return max;
  }

}  // End of namespace Luna

#endif
