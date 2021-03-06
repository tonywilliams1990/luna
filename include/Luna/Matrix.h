/// \file Matrix.h
/// A templated matrix class which constructs a Matrix as an std::vector of
/// Vectors.

#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <algorithm>
#include <random>
#include <chrono>
#include <complex>
#include <omp.h>

#include "Error.h"
#include "Vector.h"
#include "Eigen/Dense"

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

    std::size_t max_abs_in_column( std::size_t col, std::size_t start_row );

    void backsolve( Vector<T>& x );

    void Gauss_with_pivot( Vector<T>& x );

    void Gauss_with_pivot( Matrix<T>& X );

    void partial_pivot( Vector<T>& x, std::size_t k );

    void partial_pivot( Matrix<T>& X, std::size_t k );

    std::size_t LU_decomp_in_place( Matrix<T>& P );

    void rsolve( Vector<T>& b, Vector<T>& x, Matrix<T>& R );

    void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove)
    {
      unsigned int numRows = matrix.rows()-1;
      unsigned int numCols = matrix.cols();

      if( rowToRemove < numRows )
          matrix.block(rowToRemove,0,numRows-rowToRemove,numCols)
          = matrix.block(rowToRemove+1,0,numRows-rowToRemove,numCols);

      matrix.conservativeResize(numRows,numCols);
    }

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
    /// \param x The Vector which is to be multiplied
    /// \return The result vector A * x
    Vector<T> operator*( Vector<T>& x );

    /// Matrix Vector multiplication
    /// \param x The Vector which is to be multiplied
    /// \return The result vector A * x
    //Vector<std::complex<double>> operator*( Vector<std::complex<double>>& x );

    /* ----- Methods ----- */

    /// Matrix multiplication
    /// \param B The Matrix which is to be multiplied
    /// \return The result Matrix A * B
    Matrix<T> multiply( Matrix<T>& B );

    /// Matrix Vector multiplication
    /// \param x The vector which is to be multiplied
    /// \return The result vector A * x
    Vector<T> multiply( const Vector<T>& x );

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

    /// Delete a row from the matrix
    /// \param row The row index
    void delete_row( const std::size_t& row );

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

    /// Swap two elements of the Matrix
    /// \param i_1 The row index of the first element
    /// \param j_1 The column index of the first element
    /// \param i_2 The row index of the second element
    /// \param j_2 The column index of the second element
    void swap_elem( const std::size_t& i_1, const std::size_t& j_1,
                    const std::size_t& i_2, const std::size_t& j_2 );

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

    /// Set the Matrix to be the identity Matrix
    void eye();

    /// Real part of the elements of the Matrix
    /// \return Matrix of real parts of the elements
    Matrix<double> real() const;

    /// Imaginary part of the elements of the Matrix
    /// \return Matrix of imaginary parts of the elements
    Matrix<double> imag() const;

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

    /* ----- Solve linear systems ----- */

    /// Solve the system of equations Ax=b where x and b are Vectors
    /// \param b The right-hand side Vector of the system of equations
    /// \return The solution Vector
    Vector<T> solve_basic( const Vector<T>& b );

    /// Solve the system of equations AX=B where X and B are Matrices
    /// \param B The right-hand side Matrix of the system of equations
    /// \return The solution Matrix (each column is a solution Vector)
    Matrix<T> solve_basic( const Matrix<T>& B );

    /// LU decomposition of the Matrix ( with partial pivotting )
    /// \param L The lower triangular Matrix
    /// \param U The upper triangular Matrix
    /// \param P The permutation Matrix
    void LU_decomp( Matrix<T>& L, Matrix<T>& U, Matrix<T>& P );

    /// Solve the system of equations Ax=b where x and b are Vectors using LU
    /// decomposition
    /// \param b The right-hand side Vector of the system of equations
    /// \return The solution Vector
    Vector<T> solve_LU( const Vector<T>& b );

    /// Solve the system of equations AX=B where X and B are Matrices using LU
    /// decomposition
    /// \param B The right-hand side Matrix of the system of equations
    /// \return The solution Matrix (each column is a solution Vector)
    Matrix<T> solve_LU( const Matrix<T>& B );

    /// Solve the system of equations Ax=b where x and b are Vectors
    /// \param b The right-hand side Vector of the system of equations
    /// \return The solution Vector
    Vector<T> solve( const Vector<T>& b ){
      if ( ROWS < 256 ){
        return solve_LU( b );
      } else {
        return solve_parallel( b );
      }
    }

    /// Solve the system of equations AX=B where X and B are Matrices
    /// \param B The right-hand side Matrix of the system of equations
    /// \return The solution Matrix (each column is a solution Vector)
    Matrix<T> solve( const Matrix<T>& B ){ return solve_LU( B ); }

    /// QR decomposition of the Matrix A = QR
    /// \param QT The transpose of the orthogonal Matrix Q (Q^T Q = I)
    /// \param R The upper triangular Matrix R
    void QR_decomp( Matrix<T>& QT, Matrix<T>& R );

    /// Solve the system of equations Ax=b where x and b are Vectors using QR
    /// decomposition
    /// \param b The right-hand side Vector of the system of equations
    /// \return The solution Vector
    Vector<T> solve_QR( const Vector<T>& b );

    //TODO Matrix version of solve_QR

    /// \todo TODO solve_parallel use Eigen partialPivLu as in Eigen_linsys_test
    /// to solve linear systems on multiple cores -> is a bit inefficient in
    /// terms of matrix storage as a copy is made for Eigen - need to rework this.

    /// Solve the system of equations Ax=b where x and b are Vectors using Eigen
    /// parallel algorithm.
    /// \param b The right-hand side Vector of the system of equations
    /// \return The solution Vector
    Vector<T> solve_parallel( const Vector<T>& b )
    {
      if ( ROWS != b.size() )
      {
        throw Error( "solve_parallel error: ROWS != b.size() " );
      }
      if ( ROWS != COLS )
      {
        throw Error( "solve_parallel error: Matrix is not square.");
      }

      Eigen::Matrix<T, -1, 1> X, B;
      X.resize( b.size(), 1 );
      B.resize( b.size(), 1 );
      std::size_t rows( ROWS );
      std::size_t cols( COLS );

      // Create the Eigen matrix (whilst deleting the current Matrix)
      Eigen::Matrix<T, -1, -1> eigen_matrix;
      eigen_matrix.resize( rows, cols );

      for ( std::size_t i = 0; i < rows; ++i )
      {
        B( i, 0 ) = b[ i ];
        for ( std::size_t j = 0; j < cols; ++j )
        {
          //eigen_matrix( i, j ) = MATRIX[ 0 ][ j ];
          eigen_matrix( i, j ) = MATRIX[ i ][ j ];
        }
        //this -> delete_row( 0 );
      }

      // Solve the system of equations
      X = eigen_matrix.partialPivLu().solve( B );

      // Put solution into output Vector
      Vector<T> x( b.size() );
      for ( std::size_t i = 0; i < rows; ++i )
      {
        x[ i ] = X( i, 0 );
      }

      // Reconstruct the Matrix (whilst deleting the Eigen matrix)
      /*COLS = cols;
      ROWS = rows;
      MATRIX.reserve( ROWS );

      Vector<T> row( COLS );

      for ( std::size_t i = 0; i < ROWS; ++i )
      {
        for ( std::size_t j = 0; j < COLS; ++j )
        {
          row[ j ] = eigen_matrix( 0, j );
        }
        MATRIX.push_back( row );
        removeRow( eigen_matrix, 0 );
      }*/

      return x;
    }

    /* ----- Determinant ----- */

    /// Calculate the determinant of the Matrix ( via LU decomposition )
    /// \return The determinant of the Matrix
    T det();

    /* ----- Inverse ----- */

    /// Calculate the inverse of the Matrix ( via LU decomposition )
    /// \return The inverse Matrix
    Matrix<T> inverse();

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
      throw Error( "Matrix: Size error in set_row method." );
    }
    if ( row < 0 || ROWS <= row ){
      throw Error( "Matrix: Range error in set_row method." );
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
    if ( row < 0 || ROWS <= row ){
      throw Error( "Matrix: Range error in get_row method." );
    }
    Vector<T> temp;
    temp = MATRIX[ row ];
    return temp;
  }

  template <typename T>
  inline Vector<T> Matrix<T>::get_col( const std::size_t& col )
  {
    if ( col < 0 || COLS <= col ){
      throw Error( "Range error in get_col method." );
    }
    Vector<T> temp( ROWS, 0.0 );
    for ( std::size_t row = 0; row < ROWS; ++row )
    {
      temp[ row ] = MATRIX[ row ][ col ];
    }
    return temp;
  }

  template <typename T>
  inline void Matrix<T>::delete_row( const std::size_t& row )
  {
    if ( row < 0 || ROWS <= row ){
      throw Error( "Matrix: Range error in delete_row method." );
    }
    MATRIX.erase( MATRIX.begin() + row );
    --ROWS;
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
  inline void Matrix<T>::swap_elem( const std::size_t& i_1,
                                    const std::size_t& j_1,
                                    const std::size_t& i_2,
                                    const std::size_t& j_2 )
  {
    std::swap( MATRIX[ i_1 ][ j_1 ], MATRIX[ i_2 ][ j_2 ] );
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

  template <typename T>
  inline void Matrix<T>::eye()
  {
    this->fill( 0.0 );
    this->fill_diag( 1.0 );
  }

  template <typename T>
  inline Matrix<double> Matrix<T>::real() const
  {
    Matrix<double> real_part( ROWS, COLS, 0.0 );
    for (size_t i = 0; i < ROWS; ++i )
    {
      for (size_t j = 0; j < COLS; ++j )
      real_part( i, j ) = std::real( MATRIX[ i ][ j ] ) ;
    }
    return real_part;
  }

  template <typename T>
  inline Matrix<double> Matrix<T>::imag() const
  {
    Matrix<double> imag_part( ROWS, COLS, 0.0 );
    for (size_t i = 0; i < ROWS; ++i )
    {
      for (size_t j = 0; j < COLS; ++j )
      imag_part( i, j ) = std::imag( MATRIX[ i ][ j ] ) ;
    }
    return imag_part;
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

  /* ----- Solve linear systems ----- */

  template <typename T>
  inline Vector<T> Matrix<T>::solve_basic( const Vector<T>& b )
  {
    if ( ROWS != b.size() )
    {
      throw Error( "solve_basic error: ROWS != b.size() " );
    }
    if ( ROWS != COLS )
    {
      throw Error( "solve_basic error: Matrix is not square.");
    }
    Matrix<T> A( *this );
    Vector<T> x( b );
    A.Gauss_with_pivot( x );
    A.backsolve( x );
    return x;
  }

  template <typename T>
  inline Matrix<T> Matrix<T>::solve_basic( const Matrix<T>& B )
  {
    if ( ROWS != B.rows() )
    {
      throw Error( "solve_basic error: ROWS != B.rows() " );
    }
    if ( ROWS != COLS )
    {
      throw Error( "solve_basic error: Matrix is not square.");
    }
    Matrix<T> A( *this );
    Matrix<T> X( B );
    A.Gauss_with_pivot( X );

    Vector<T> x( ROWS );

    for ( std::size_t j = 0; j < B.COLS; ++j )
    {
      x = X.get_col( j );
      A.backsolve( x );
      X.set_col( j, x );
    }
    return X;
  }

  template <typename T>
  inline void Matrix<T>::LU_decomp( Matrix<T>& L, Matrix<T>& U, Matrix<T>& P )
  {
    if ( ROWS != COLS )
    {
      throw Error( "LU_decomp error: Matrix is not square.");
    }

    Matrix<T> temp( *this );
    std::size_t pivots;
    pivots = temp.LU_decomp_in_place( P );
    L.resize( ROWS, ROWS );
    L.fill( 0 );
    U.resize( ROWS, ROWS );
    U.fill( 0 );

    for ( std::size_t i = 0; i < ROWS; i++ )
    {
      for ( std::size_t j = 0; j < ROWS; j++ )
      {
        if ( i < j )
        {
          U( i, j ) = temp( i, j );
        }
        else if ( i == j )
        {
          L( i, j ) = 1;
          U( i, j ) = temp( i, j );
        }
        else
        {
          L( i, j ) = temp( i, j );
        }
      }
    }
  }

  template <typename T>
  inline Vector<T> Matrix<T>::solve_LU( const Vector<T>& b )
  {
    if ( ROWS != b.size() )
    {
      throw Error( "solve_LU error: ROWS != b.size() " );
    }
    if ( ROWS != COLS )
    {
      throw Error( "solve_LU error: Matrix is not square.");
    }
    Matrix<T> LU( *this );
    Matrix<T> P;
    std::size_t pivots;
    pivots = LU.LU_decomp_in_place( P );
    Vector<T> x( b );
    x = P * x;

    for ( std::size_t i = 0; i < ROWS; i++ ) {
        for ( std::size_t k = 0; k < i; k++ )
            x[ i ] -= LU.MATRIX[ i ][ k ] * x[ k ];
    }
    LU.backsolve( x );
    return x;
  }

  template <typename T>
  inline Matrix<T> Matrix<T>::solve_LU( const Matrix<T>& B )
  {
    if ( ROWS != B.rows() )
    {
      throw Error( "solve_basic error: ROWS != B.rows() " );
    }
    if ( ROWS != COLS )
    {
      throw Error( "solve_basic error: Matrix is not square.");
    }
    Matrix<T> LU( *this );
    Matrix<T> P;
    std::size_t pivots;
    pivots = LU.LU_decomp_in_place( P );
    Matrix<T> X( B );
    X = P * X;
    Vector<T> x( ROWS );
    for ( std::size_t j = 0; j < B.COLS; ++j )
    {
      x = X.get_col( j );
      for ( std::size_t i = 0; i < ROWS; i++ ) {
          for ( std::size_t k = 0; k < i; k++ )
              x[ i ] -= LU.MATRIX[ i ][ k ] * x[ k ];
      }
      LU.backsolve( x );
      X.set_col( j, x );
    }
    return X;
  }

  template <typename T>
  inline void Matrix<T>::QR_decomp( Matrix<T>& QT, Matrix<T>& R )
  {
    if ( ROWS != COLS )
    {
      throw Error( "QR_decomp error: Matrix is not square.");
    }
    R.MATRIX = MATRIX;
    R.ROWS = ROWS;
    R.COLS = COLS;

    QT.resize( ROWS, ROWS );

    bool singular( false );

    int i, j, k;
    Vector<T> c( ROWS ), d( ROWS );
    T sigma, sum, tau;
    double scale;

    for ( k = 0; k < ROWS - 1; k++ ) {
      scale = 0.0;
      for ( i = k; i < ROWS; i++ ) {
        scale = std::max( scale, std::abs( R.MATRIX[ i ][ k ] ) );
      }
      if ( scale == 0.0 ) {
        singular = true;
        c[k] = d[k] = 0.0;
      } else {
        for ( i = k; i < ROWS; i++ ){
          R.MATRIX[ i ][ k ] /= scale;
        }
        for ( sum = 0.0, i = k; i < ROWS; i++ ){
          sum += R.MATRIX[ i ][ k ] * R.MATRIX[ i ][ k ];
        }
        T a = std::sqrt( sum ); //TODO check this is okay
        T b = R.MATRIX[ k ][ k ];
        sigma = std::real(b) >= 0 ? ( std::real(a) >= 0 ? a : -a) :
                                    ( std::real(a) >= 0 ? -a : a); // sign(a,b)
        R.MATRIX[ k ][ k ] += sigma;
        c[k] = sigma * R.MATRIX[ k ][ k ];
        d[k] = - scale * sigma;
        for ( j = k + 1; j < ROWS; j++ ) {
          for ( sum = 0.0 , i = k; i < ROWS; i++ ){
            sum += R.MATRIX[ i ][ k ] * R.MATRIX[ i ][ j ];
          }
          tau = sum / c[ k ];
          for (i = k; i < ROWS; i++ ){
            R.MATRIX[ i ][ j ] -= tau * R.MATRIX[ i ][ k ];
          }
        }
      }
    }

    d[ ROWS - 1 ] = R.MATRIX[ ROWS - 1 ][ ROWS - 1 ];
    if ( d[ ROWS - 1 ] == 0.0) singular=true;
    QT.eye();

    for ( k = 0; k < ROWS - 1; k++ ) {
      if ( c[ k ] != 0.0 ) {
        for ( j = 0; j < ROWS; j++ ) {
          sum = 0.0;
          for ( i = k; i < ROWS; i++ ) {
            sum += R.MATRIX[ i ][ k ] * QT.MATRIX[ i ][ j ];
          }
          sum /= c[ k ];
          for ( i = k; i < ROWS; i++ ) {
            QT.MATRIX[ i ][ j ] -= sum * R.MATRIX[ i ][ k ];
          }
        }
      }
    }

    for ( i = 0; i < ROWS; i++ ) {
      R.MATRIX[ i ][ i ] = d[ i ];
      for ( j = 0; j < i; j++ ){
        R.MATRIX[ i ][ j ] = 0.0;
      }
    }

    if ( singular ) {
      std::cout << " *** WARNING! Matrix is singular in QR decomposition."
                << std::endl;
    }
  }

  template <typename T>
  inline Vector<T> Matrix<T>::solve_QR( const Vector<T>& b )
  {
    if ( ROWS != b.size() )
    {
      throw Error( "solve_QR error: ROWS != b.size() " );
    }
    if ( ROWS != COLS )
    {
      throw Error( "solve_QR error: Matrix is not square.");
    }
    Matrix<T> QT, R;
    this->QR_decomp( QT, R );
    Vector<T> x( ROWS, 0.0 );
    x = QT.multiply( b );
    this->rsolve( x, x, R );            // Solve R * x = Q^T * b
    return x;
  }

  /* ----- Determinant ----- */

  template <typename T>
  inline T Matrix<T>::det()
  {
    T det( 1 );
    Matrix<T> temp( *this );
    Matrix<T> P;
    std::size_t pivots;
    pivots = temp.LU_decomp_in_place( P );
    for ( std::size_t i = 0; i < ROWS; i++ )
    {
      det *= temp.MATRIX[ i ][ i ];
    }
    if ( pivots % 2 == 0 )
    {
      return det;
    }
    else
    {
      return - det;
    }
  }

  /* ----- Inverse ----- */

  template <typename T>
  inline Matrix<T> Matrix<T>::inverse()
  {
    if ( ROWS != COLS )
    {
      throw Error( "Matrix inverse error: Matrix is not square.");
    }
    Matrix<T> inv;
    Matrix<T> LU( *this );
    std::size_t pivots;
    pivots = LU.LU_decomp_in_place( inv );
    for ( int j = 0; j < ROWS; j++ )
    {
        for ( int i = 0; i < ROWS; i++ )
        {
            for ( int k = 0; k < i; k++ )
            {
                inv.MATRIX[ i ][ j ] -= LU.MATRIX[ i ][ k ]
                                      * inv.MATRIX[ k ][ j ];
            }
        }

        for (int i = ROWS - 1; i >= 0; i--)
        {
            for (int k = i + 1; k < ROWS; k++)
            {
                inv.MATRIX[ i ][ j ] -= LU.MATRIX[ i ][ k ]
                                      * inv.MATRIX[ k ][ j ];
            }
            inv.MATRIX[ i ][ j ] = inv.MATRIX[ i ][ j ] / LU.MATRIX[ i ][ i ];
        }
    }
    return inv;
  }

  /* ----- Private ----- */

  template <typename T>
  inline std::size_t Matrix<T>::max_abs_in_column( std::size_t col,
                                                   std::size_t start_row )
  {
    std::size_t max_index = 0;
    double max = -1.0e+10;

    for ( std::size_t i = start_row; i < ROWS; i++ ){
      if ( max < std::abs( MATRIX[ i ][ col ] ) ){
        max = std::abs( MATRIX[ i ][ col ] );
        max_index = i;
      }
    }
    return max_index;
  }

  template <typename T>
  inline void Matrix<T>::backsolve( Vector<T>& x )
  {
    x[ ROWS - 1 ] = x[ ROWS - 1 ] / MATRIX[ ROWS - 1 ][ ROWS - 1 ];

    for ( std::size_t n = 2; n < ROWS + 1; ++n )
    {
      for ( std::size_t j = ROWS - n + 1; j < ROWS; ++j )
      {
        x[ ROWS - n ] -= MATRIX[ ROWS - n ][ j ] * x[ j ];
      }
      x[ ROWS - n ] = x[ ROWS - n ] / MATRIX[ ROWS - n ][ ROWS - n ];
    }
  }

  template <typename T>
  inline void Matrix<T>::Gauss_with_pivot( Vector<T>& x )
  {
    std::size_t k;
    for ( k = 0; k < ROWS - 1; ++k )
    {
      this->partial_pivot( x, k );
      for ( std::size_t i = k + 1; i < ROWS; ++i )
      {
        T elem = MATRIX[ i ][ k ] / MATRIX[ k ][ k ];
        for ( std::size_t j = k; j < ROWS; ++j )
        {
          MATRIX[ i ][ j ] -= elem * MATRIX[ k ][ j ];
        }
        x[ i ] -= elem * x[ k ];
      }
    }
  }

  template <typename T>
  inline void Matrix<T>::Gauss_with_pivot( Matrix<T>& X )
  {
    std::size_t k;
    for ( k = 0; k < ROWS - 1; ++k )
    {
      this->partial_pivot( X, k );
      for ( std::size_t i = k + 1; i < ROWS; ++i )
      {
        T elem = MATRIX[ i ][ k ] / MATRIX[ k ][ k ];
        for ( std::size_t j = k; j < ROWS; ++j )
        {
          MATRIX[ i ][ j ] -= elem * MATRIX[ k ][ j ];
        }
        X.MATRIX[ i ] -= elem * X.MATRIX[ k ];
      }
    }
  }

  template <typename T>
  inline void Matrix<T>::partial_pivot( Vector<T>& x, std::size_t k )
  {
    std::size_t pivot;
    pivot = this->max_abs_in_column( k, k );
    this->swap_rows( pivot, k );
    x.swap( pivot, k );
  }

  template <typename T>
  inline void Matrix<T>::partial_pivot( Matrix<T>& X, std::size_t k )
  {
    std::size_t pivot;
    pivot = this->max_abs_in_column( k, k );
    this->swap_rows( pivot, k );
    X.swap_rows( pivot, k );
  }

  template <typename T>
  inline std::size_t Matrix<T>::LU_decomp_in_place( Matrix<T>& P )
  {

    if ( ROWS != COLS )
    {
      throw Error( "LU_decomp_in_place error: Matrix is not square.");
    }

    std::size_t i, j, k, imax, pivots( 0 );
    double maxA, absA;

    Vector<T> row( ROWS, 0.0 );

    P.resize( ROWS, ROWS );
    P.fill( 0 );
    P.fill_diag( 1 );

    for (i = 0; i < ROWS; i++) {
        maxA = 0.0;
        imax = i;

        for (k = i; k < ROWS; k++) {
          if ( ( absA = std::abs( MATRIX[ k ][ i ] ) ) > maxA ) {
              maxA = absA;
              imax = k;
          }
        }

        if ( maxA < 1e-40 )
        {
          throw Error( "LU_decomp_in_place error: Matrix is singular." );
        }

        if (imax != i) {
            P.swap_rows( i, imax );
            this->swap_rows( i, imax );
            pivots++;
        }

        for (j = i + 1; j < ROWS; j++) {
            MATRIX[ j ][ i ] /= MATRIX[ i ][ i ];

            for (k = i + 1; k < ROWS; k++)
                MATRIX[ j ][ k ] -= MATRIX[ j ][ i ] * MATRIX[ i ][ k ];
        }
    }
    return pivots;
  }

  template <typename T>
  inline void Matrix<T>::rsolve( Vector<T>& b, Vector<T>& x, Matrix<T>& R )
  {
    int i, j;
    T sum;
    for ( i = R.ROWS - 1; i >= 0; i-- ) {
      sum = b[ i ];
      for ( j = i + 1; j < R.ROWS; j++ ){
        sum -= R.MATRIX[ i ][ j ] * x[ j ];
      }
      x[ i ]=sum / R.MATRIX[ i ][ i ];
    }
  }


}  // End of namespace Luna

#endif
