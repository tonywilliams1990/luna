/// \file SparseMatrix.h
/// A templated sparse matrix class
/// TODO

#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H

#include <vector>
#include <algorithm>
#include <tuple>
#include <random>
#include <chrono>
#include <complex>
#include <omp.h>

#include "Error.h"
#include "Vector.h"
#include "Triplet.h"

namespace Luna
{

  // Templated sparse matrix class
  template <class T>

  /// A SparseMatrix class for use with double and std::complex<double>
  /// TODO column compressed storage?
  class SparseMatrix
  {
  private:
    std::size_t ROWS;					          // Number of rows
    std::size_t COLS;					          // Number of columns
    std::size_t N;                      // Number of nonzero elements

    Vector<T> VAL;                      // Vector of nonzero values
    Vector<std::size_t> ROW_INDEX;      // Row indices of nonzero values
    Vector<std::size_t> COL_START;      // Pointers to start of columns

  public:

    SparseMatrix() : ROWS( 0 ), COLS( 0 ), N( 0 ),
                  COL_START(), ROW_INDEX(), VAL()
    {}

    SparseMatrix( const std::size_t& rows, const std::size_t& cols,
                  const std::size_t& n ) : ROWS( rows ), COLS( cols ), N( n ),
                  COL_START( cols + 1, 0 ), ROW_INDEX( n, 0 ), VAL( n, 0.0 )
    {}

    /// Constructor from Vector of Triplets
    /// \param rows The number of rows in the matrix.
    /// \param cols The number of columns in the matrix.
    /// \param triplets Vector of Triplets
    SparseMatrix( const std::size_t& rows, const std::size_t& cols,
                  Vector< Triplet<T> >& triplets );

    /// Constructor from Vectors
    /// \param rows The number of rows in the matrix.
    /// \param cols The number of columns in the matrix.
    /// \param val Vector of nonzero values
    /// \param row_index Vector of row indices
    /// \param col_start Vector of column start pointers
    SparseMatrix( const std::size_t& rows, const std::size_t& cols,
                  const Vector<T>& val, const Vector<std::size_t>& row_index,
                  const Vector<std::size_t>& col_start );

    /// Copy constructor
    /// \param source The source SparseMatrix to be copied
    SparseMatrix( const SparseMatrix<T>& source );

    /// Destructor
    ~SparseMatrix(){}

    /* ----- Methods ----- */

    /// Return the number of rows in the SparseMatrix
    /// \return The number of rows in the SparseMatrix
    std::size_t rows() const;

    /// Return the number of columns in the SparseMatrix
    /// \return The number of columns in the SparseMatrix
    std::size_t cols() const;

    /// Return the Vector of nonzero values
    /// \return The Vector of nonzero values
    Vector<T> val() const;

    /// Return the Vector of row indices
    /// \return The Vector of row indices
    Vector<std::size_t> row_index() const;

    /// Return the Vector of column start pointers
    /// \return The column start Vector
    Vector<std::size_t> col_start() const;

    /// Return the number of nonzero values
    /// \return The number of nonzero values
    std::size_t nonzero();

    /// The Vector of column indices ( triplet form )
    /// \return A Vector of column indices relating to each value
    Vector<std::size_t> col_index() const;

    /// Calculate column start Vector for a given column index Vector
    /// \param col_index A column index Vector
    /// \return col_start A column start Vector
    Vector<std::size_t> col_start_from_index(
                                        const Vector<std::size_t>& col_index );

    /// Insert a nonzero element into the SparseMatrix
    /// \param row The row index of the element
    /// \param col The column index of the element
    /// \param val The value of the element
    void insert( const std::size_t& row, const std::size_t& col, const T& val );

    /// Scale the SparseMatrix by a constant factor
    /// \param scalar The scalar to multiply the SparseMatrix by
    void scale( const T& scalar );

    /// Get a specific value if it exists
    /// \param row The row index
    /// \param col The column index
    /// \return The value stored at that location
    const T& get( const std::size_t& row, const std::size_t& col ) const;


    /// Multiply the SparseMatrix A by a Vector to the right
    /// \param x The vector which is to be multiplied
    /// \return The result vector A * x
    Vector<T> multiply( const Vector<T>& x );

    /// Multiply the transpose of the SparseMatrix A^T by a Vector to the right
    /// \param x The vector which is to be multiplied
    /// \return The result vector A^T * x
    Vector<T> transpose_multiply( const Vector<T>& x );

    /// Return the transpose of the SparseMatrix
    /// \return The transpose of the SparseMatrix
    SparseMatrix<T> transpose() const;

  };	// End of class SparseMatrix

  template <typename T>
  inline SparseMatrix<T>::SparseMatrix( const std::size_t& rows,
                                        const std::size_t& cols,
                                        Vector< Triplet<T> >& triplets )
  {
    ROWS = rows;
    COLS = cols;
    // Sort the triplets for column major ordering
    std::sort( triplets.begin(), triplets.end() );

    Vector<std::size_t> col_index;
    std::size_t row, col;
    T val;
    N = 0;

    for ( std::size_t i = 0; i < triplets.size(); i++ )
    {
      row = triplets[ i ].get_row();
      col = triplets[ i ].get_col();
      if ( row<0 || ROWS<=row )	{
        throw Error( "SparseMatrix range error: dimension 1 in triplets." );
      }
      if ( col<0 || COLS<=col )	{
        throw Error( "SparseMatrix range error: dimension 2 in triplets." );
      }
      val = triplets[ i ].get_val();

      // Check there are no duplicates
      for ( std::size_t k = 0; k < N; k++ )
      {
        if ( ROW_INDEX[ k ] == row && col_index[ k ] == col )
        {
          throw Error( "SparseMatrix: duplicate entry in triplet list." );
        }
      }
      // Push into Vectors
      ROW_INDEX.push_back( row );
      col_index.push_back( col );
      VAL.push_back( val );
      N++;
    }
    // Convert col_index to COL_START
    COL_START = this->col_start_from_index( col_index );
  }

  template <typename T>
  inline SparseMatrix<T>::SparseMatrix( const std::size_t& rows,
               const std::size_t& cols, const Vector<T>& val,
               const Vector<std::size_t>& row_index,
               const Vector<std::size_t>& col_start )
  {
    ROWS = rows;
    COLS = cols;
    VAL = val;
    ROW_INDEX = row_index;
    COL_START = col_start;
    N = col_start[ col_start.size() - 1 ];
  }

  template <typename T>
  inline SparseMatrix<T>::SparseMatrix( const SparseMatrix<T>& source )
  {
    *this = source;
  }

  /* ----- Methods ----- */

  template <typename T>
  inline std::size_t SparseMatrix<T>::rows() const
  {
    return ROWS;
  }

  template <typename T>
  inline std::size_t SparseMatrix<T>::cols() const
  {
    return COLS;
  }

  template <typename T>
  inline Vector<T> SparseMatrix<T>::val() const
  {
    return VAL;
  }

  template <typename T>
  inline Vector<std::size_t> SparseMatrix<T>::row_index() const
  {
    return ROW_INDEX;
  }

  template <typename T>
  inline Vector<std::size_t> SparseMatrix<T>::col_start() const
  {
    return COL_START;
  }

  template <typename T>
  inline std::size_t SparseMatrix<T>::nonzero()
  {
    return N;
  }

  template <typename T>
  inline Vector<std::size_t> SparseMatrix<T>::col_index() const
  {
    Vector<std::size_t> temp;
    if ( N == 0 ){ return temp; }
    if ( COL_START.size() < COLS + 1 )
    {
      throw Error( "Some columns have no entries." );
    }
    Vector<std::size_t> gaps( COL_START.size() - 1 );
    for ( std::size_t k = 0; k < gaps.size(); k++ )
    {
      gaps[ k ] = COL_START[ k + 1 ] - COL_START[ k ];
      for ( std::size_t j = 0; j < gaps[ k ]; j ++ )
      {
        temp.push_back( k );
      }
    }
    return temp;
  }

  template <typename T>
  inline Vector<std::size_t> SparseMatrix<T>::col_start_from_index(
                                      const Vector<std::size_t>& col_index )
  {
    Vector<std::size_t> temp;
    std::size_t max_col_index;
    max_col_index = col_index[ col_index.size() - 1 ];
    std::size_t count( 0 );
    Vector<std::size_t> gaps;
    for ( std::size_t j = 0; j < max_col_index + 1; j++ )
    {
      count = std::count( col_index.begin(),
                          col_index.begin() + col_index.size(), j );
      gaps.push_back( count );
    }
    temp.push_back( 0 );
    for ( std::size_t k = 1; k < gaps.size() + 1; k++ )
    {
      temp.push_back( temp[ k - 1 ] + gaps[ k - 1 ] );
    }
    return temp;
  }

  template <typename T>
  inline void SparseMatrix<T>::insert( const std::size_t& row,
                                       const std::size_t& col, const T& val )
  {
    if ( row<0 || ROWS<=row )	{
      throw Error( "SparseMatrix range error: dimension 1." );
    }
    if ( col<0 || COLS<=col )	{
      throw Error( "SparseMatrix range error: dimension 2." );
    }
    if ( COL_START.size() < COLS + 1 )
    {
      std::string problem;
      problem  = "SparseMatrix insert error: Some columns have no entries.\n";
      problem += "Use Vector<Triplet> instead.";
      throw Error( problem );
    }
    // Convert to triplet format
    Vector<std::size_t> col_index = this->col_index();
    // Check if element already exists
    for ( std::size_t k = 0; k < N; k++)
    {
      if ( ROW_INDEX[ k ] == row && col_index[ k ] == col )
      {
        VAL[ k ] = val;
        return;
      }
    }
    if ( N == 0 )
    {
      ROW_INDEX.push_back( row );
      VAL.push_back( val );
      COL_START.push_back( 1 );
      N++;
      return;
    }
    // If not insert new element
    std::size_t find_col = col_index.find( col );
    col_index.insert( find_col, col );
    std::size_t find_row( find_col );
    while ( ROW_INDEX[ find_row ] < row )
    {
      find_row++;
    }
    ROW_INDEX.insert( find_row, row );
    VAL.insert( find_row, val );
    // Convert back to compressed column format
    COL_START = this->col_start_from_index( col_index );
    N++;
    return;
  }

  template <typename T>
  inline void SparseMatrix<T>::scale( const T& scalar )
  {
    VAL *= scalar;
  }

  template <typename T>
  inline const T& SparseMatrix<T>::get( const std::size_t& row,
                                        const std::size_t& col ) const
  {
    if ( row<0 || ROWS<=row )	{
      throw Error( "SparseMatrix range error in get method: dimension 1." );
    }
    if ( col<0 || COLS<=col )	{
      throw Error( "SparseMatrix range error in get method: dimension 2." );
    }
    if ( COL_START.size() < COLS + 1 )
    {
      std::string problem;
      problem  = "SparseMatrix get error: Some columns have no entries.";
      throw Error( problem );
    }

    Vector<std::size_t> col_index = this->col_index();

    for ( std::size_t k = 0; k < N; k++)
    {
      if ( ROW_INDEX[ k ] == row && col_index[ k ] == col )
      {
        return VAL[ k ];
      }
      else
      {
        throw Error( "SparseMatrix get error: entry does not exist." );
      }
    }
  }

  template <typename T>
  inline Vector<T> SparseMatrix<T>::multiply( const Vector<T>& x )
  {
    if ( x.size() != COLS )
    {
      throw Error( "SparseMatrix multiply: dimensions do not agree." );
    }
    Vector<T> y( ROWS, 0 );
    for ( std::size_t j = 0; j < COLS; j++ )
    {
      for ( std::size_t i = COL_START[ j ]; i < COL_START[ j + 1 ]; i++ )
      {
        y[ ROW_INDEX[ i ] ] += VAL[ i ] * x[ j ];
      }
    }
    return y;
  }

  template <typename T>
  inline Vector<T> SparseMatrix<T>::transpose_multiply( const Vector<T>& x )
  {
    if ( x.size() != ROWS )
    {
      std::string problem;
      problem  = "SparseMatrix transpose_multiply: dimensions do not agree.";
      throw Error( problem );
    }
    Vector<T> y( COLS, 0 );
    for ( std::size_t i = 0; i < COLS; i++ )
    {
      for ( std::size_t j = COL_START[ i ]; j < COL_START[ i + 1 ]; j++ )
      {
        y[ i ] += VAL[ j ] * x[ ROW_INDEX[ j ] ];
      }
    }
    return y;
  }

  template <typename T>
  inline SparseMatrix<T> SparseMatrix<T>::transpose() const
  {
    std::size_t i, j, k, index, m = ROWS, n = COLS;
    SparseMatrix<T> at( n, m, N );

    Vector<int> count( m, 0 );
    for ( i = 0; i < n; i++ )
    {
      for ( j = COL_START[ i ]; j < COL_START[ i + 1 ]; j++ )
      {
        k = ROW_INDEX[ j ];
        count[ k ]++;
      }
    }
    for ( j = 0; j < m; j++ )
    {
      at.COL_START[ j + 1 ] = at.COL_START[ j ] + count[ j ];
    }
    for ( j = 0; j < m; j++ )
    {
      count[ j ] = 0;
    }
    for ( i = 0; i < n; i++ )
    {
      for ( j = COL_START[ i ]; j < COL_START[ i + 1 ]; j++ )
      {
        k = ROW_INDEX[ j ];
        index = at.COL_START[ k ] + count[ k ];
        at.ROW_INDEX[ index ] = i;
        at.VAL[ index ] = VAL[ j ];
        count[ k ]++;
      }
    }
    return at;
  }


}  // End of namespace Luna

#endif
