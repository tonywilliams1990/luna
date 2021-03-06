/// \file SparseMatrix.h
/// A SparseMatrix class for use with double and std::complex<double>

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
#include "Timer.h"

namespace Luna
{

  // Templated sparse matrix class
  template <class T>

  /// A templated sparse matrix class for solving sparse linear systems.
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

    /// Solve the system of equations Ax=b using the biconjugate gradient method
    /// \param b The right-hand side Vector of the system of equations
    /// \param x Initial guess x_0 the solution is also stored here
    /// \param max_iter The maximum number of iterations to be performed
    /// ( the number of iterations actually taken is returned here )
    /// \param tol Tolerance for convergence ( the estimated error is returned
    /// here )
    /// \param itol Convergence test specifier
    /// \return Success or failure code ( 0 = success, 1 = exceeded max_iter )
    int solve_BiCG( const Vector<T>& b, Vector<T>& x, int& max_iter,
                    double& tol, int itol = 1 );

    /// Diagonal preconditioner for the solve_BiCG method
    /// \param b Right hand side Vector
    /// \param x Solution Vector
    void diagonal_preconditioner( const Vector<T>& b, Vector<T>& x );

    /// Identity preconditioner for the solve_BiCG method
    /// \param b Right hand side Vector
    /// \param x Solution Vector
    void identity_preconditioner( const Vector<T>& b, Vector<T>& x );

    /// Solve the system of equations Ax=b using the stabilised biconjugate
    /// gradient method BiCGSTAB
    /// \param b The right-hand side Vector of the system of equations
    /// \param x Initial guess x_0 the solution is also stored here
    /// \param max_iter The maximun number of iterations ( the number of
    /// iterations performed will also be returned here )
    /// \param tol Tolerance for convergence ( the error residual is also
    /// returned here )
    /// \return Success or failure code ( 0 = success, 1 = exceeded max_iter )
    int solve_BiCGSTAB( const Vector<T>& b, Vector<T>& x, int& max_iter,
                        double& tol );

    /// Solve the system of equations Ax=b using the conjugate gradient method
    /// \param b The right-hand side Vector of the system of equations
    /// \param x Initial guess x_0 the solution is also stored here
    /// \param max_iter The maximun number of iterations ( the number of
    /// iterations performed will also be returned here )
    /// \param tol Tolerance for convergence ( the error residual is also
    /// returned here )
    /// \return Success or failure code ( 0 = success, 1 = exceeded max_iter )
    int solve_CG( const Vector<T>& b, Vector<T>& x, int& max_iter,
                  double& tol );

    /// Solve the system of equations Ax=b using the Quasi-minimal residual
    /// method
    /// \param b The right-hand side Vector of the system of equations
    /// \param x Initial guess x_0 the solution is also stored here
    /// \param max_iter The maximun number of iterations ( the number of
    /// iterations performed will also be returned here )
    /// \param tol Tolerance for convergence ( the error residual is also
    /// returned here )
    /// \return Success or failure code ( 0 = success, 1 = exceeded max_iter )
    int solve_QMR( const Vector<T>& b, Vector<T>& x, int& max_iter,
                  double& tol );

    //TODO GMRES, SVD to find condition number ?
    //TODO General solve method - choose method by string
    //TODO Sparse matrix-matrix multiplication ?
    //TODO Preconditioning ?


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

      // Check there are no duplicates (really slows things down)
      /*for ( std::size_t k = 0; k < N; k++ )
      {
        if ( ROW_INDEX[ k ] == row && col_index[ k ] == col )
        {
          throw Error( "SparseMatrix: duplicate entry in triplet list." );
        }
      }*/
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
    Vector<std::size_t> col_start( COLS + 1, 0 );
    // Computer number of elements in each column
    for ( std::size_t n = 0; n < N; ++n )
    {
      col_start[ col_index[ n ] ]++;
    }
    std::size_t sum( 0 );
    std::size_t k;
    // Cumulative sum
    for ( k = 0; k < COLS; ++k )
    {
      std::size_t ck = col_start[ k ];
      col_start[ k ] = sum;
      sum += ck;
    }
    col_start[ COLS ] = sum;
    return col_start;
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

    std::size_t i, j;
    T xj;

    for ( j = 0; j < COLS; ++j )
    {
      xj = x[ j ];
      for ( i = COL_START[ j ]; i < COL_START[ j + 1 ]; ++i )
      {
        y[ ROW_INDEX[ i ] ] += VAL[ i ] * xj;
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

    std::size_t i, j;

    for ( i = 0; i < COLS; i++ )
    {
      for ( j = COL_START[ i ]; j < COL_START[ i + 1 ]; j++ )
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

  template <typename T>
  inline int SparseMatrix<T>::solve_BiCG( const Vector<T>& b, Vector<T>& x,
                                          int& max_iter, double& tol, int itol )
  {
    if ( ROWS != b.size() )
    {
      throw Error( "solve_BiCG error: ROWS != b.size() " );
    }
    if ( ROWS != COLS )
    {
      throw Error( "solve_BiCG error: Matrix is not square.");
    }
    if ( b.size() != x.size() )
    {
      throw Error( "solve_BiCG error: b.size() != x.size() " );
    }

    T alpha, beta, rho_1, rho_2 = 1.0;
    double bnrm, znrm, err;
    int j, n = b.size();
    Vector<T> p( n ), pp( n ), r( n ), rr( n ), z( n ), zz( n );

    r = this->multiply( x );
    r = b - r;
    rr = r;

    if ( itol == 1 ) {
      bnrm = b.norm_2();
      //this->diagonal_preconditioner( r, z );           // replaces asolve in NR
      this->identity_preconditioner( r, z );
    }
    else if ( itol == 2 ) {
      //this->diagonal_preconditioner( b, z );
      this->identity_preconditioner( b, z );
      bnrm = z.norm_2();
      //this->diagonal_preconditioner( r, z );
      this->identity_preconditioner( r, z );
    }
    else {
      throw Error( "SparseMatrix solve_BiCG error: illegal itol." );
    }

    // Main loop
    int iter( 0 );
    while ( iter < max_iter )
    {
      ++iter;
      //this->diagonal_preconditioner( rr, zz );
      this->identity_preconditioner( rr, zz );

      rho_1 = z.dot( rr );

      if ( iter == 1 ) {
        p = z;
        pp = zz;
      } else {
        beta = rho_1 / rho_2;
        p = beta * p + z;
        pp = beta * pp + zz;
      }

      z = this->multiply( p );
      alpha = rho_1 / z.dot( pp );
      zz = this->transpose_multiply( pp );
      x += alpha * p;
      r -= alpha * z;
      rr -= alpha * zz;

      //this->diagonal_preconditioner( r, z );
      this->identity_preconditioner( r, z );

      rho_2 = rho_1;

      if ( itol == 1 ) {
        err = r.norm_2() / bnrm;
      } else if ( itol == 2 ) {
        err = z.norm_2() / bnrm;
      }
      if ( err <= tol )
      {
        max_iter = iter;
        tol = err;
        return 0;
      }
    } // end of while loop
    return 1;
  }

  template <typename T>
  inline void SparseMatrix<T>::diagonal_preconditioner( const Vector<T>& b,
                                                       Vector<T>& x )
  {
    Vector<std::size_t> col_index = this->col_index();
    for ( std::size_t i = 0; i < ROWS; i++ )
    {
      for ( std::size_t k = 0; k < N; k++)
      {
        if ( ROW_INDEX[ k ] == i && col_index[ k ] == i ) {
          x[ i ] = b[ i ] / VAL[ k ];
        }
        else {
          x[ i ] = b[ i ];
        }
      }
    }
  }

  template <typename T>
  inline void SparseMatrix<T>::identity_preconditioner( const Vector<T>& b,
                                                        Vector<T>& x )
  {
    for ( std::size_t i = 0; i < ROWS; i++ )
    {
      x[ i ] =  b[ i ];
    }
  }

  template <typename T>
  inline int SparseMatrix<T>::solve_BiCGSTAB( const Vector<T>& b, Vector<T>& x,
                                              int& max_iter, double& tol )
  {
    // https://math.nist.gov/iml++
    double resid;
    Vector<T> p( x.size(), 0.0 ), phat( x.size(), 0.0 );
    Vector<T> s( x.size(), 0.0 ), shat( x.size(), 0.0 );
    Vector<T> v( x.size(), 0.0 );
    Vector<T> temp, t;

    T rho_1 = 1.0;
    T rho_2 = 1.0;
    T alpha = 1.0;
    T beta;
    T omega = 1.0;

    double normb = b.norm_2();
    Vector<T> r;
    r = this->multiply( x );
    r = b - r;
    Vector<T> rtilde = r;

    if ( normb == 0.0 )
    {
      normb = 1.;
    }
    if ( ( resid = r.norm_2() / normb ) <= tol )
    {
      tol = resid;
      max_iter = 0;
      return 0;
    }

    for ( int i = 1; i <= max_iter; i++ )
    {
      rho_1 = rtilde.dot( r );
      if ( rho_1 == 0.0 ) {
        tol = r.norm_2() / normb;
        return 2;
      }
      if ( i == 1 ) {
        p = r;
      } else {
        beta = ( rho_1 / rho_2 ) * ( alpha / omega );
        temp = p - omega * v;
        p = r + beta * temp;
        //p = r + beta * ( p - omega * v );
      }
      //phat = p; //could have preconditioner here phat = M.solve(p);
      this->identity_preconditioner( p, phat );
      v = this->multiply( phat );
      alpha = rho_1 / rtilde.dot( v );
      s = r - alpha * v;
      if ( ( resid = s.norm_2() / normb ) <= tol ) {
        x += alpha * phat;
        tol = resid;
        max_iter = i;
        return 0;
      }
      //shat = s; //could have preconditioner here shat = M.solve(s);
      this->identity_preconditioner( s, shat );
      t = this->multiply( shat );
      omega = t.dot( s ) / t.dot( t );
      x += alpha * phat;
      x += omega * shat;
      r = s - omega * t;

      rho_2 = rho_1;

      if ( ( resid = r.norm_2() / normb ) < tol ) {
        tol = resid;
        max_iter = i;
        return 0;
      }
      if ( omega == 0.0 ) {
        tol = r.norm_2() / normb;
        return 3;
      }
    }

    tol = resid;
    return 1;
  }

  template <typename T>
  inline int SparseMatrix<T>::solve_CG( const Vector<T>& b, Vector<T>& x,
                                        int& max_iter, double& tol )
  {
    double resid;
    Vector<T> p( x.size(), 0.0 ), z( x.size(), 0.0 ), q( x.size(), 0.0 );
    Vector<T> temp;

    T rho = 1.0;
    T rho_1 = 1.0;
    T alpha = 1.0;
    T beta;

    double normb = b.norm_2();
    Vector<T> r;
    r = this->multiply( x );
    r = b - r;

    if ( normb == 0.0 )
    {
      normb = 1.;
    }
    if ( ( resid = r.norm_2() / normb ) <= tol )
    {
      tol = resid;
      max_iter = 0;
      return 0;
    }

    for ( int i = 1; i <= max_iter; i++ )
    {
      //this->diagonal_preconditioner( r, z );
      //this->identity_preconditioner( r, z );
      z = r; // Could have preconditioner here?
      rho = r.dot( z );

      if ( i == 1 ) {
        p = z;
      } else {
        beta = rho / rho_1;
        p = z + beta * p;
      }

      q = this->multiply( p );
      alpha = rho / p.dot( q );

      x += alpha * p;
      r -= alpha * q;

      if ( ( resid = r.norm_2() / normb ) <= tol )
      {
        tol = resid;
        max_iter = i;
        return 0;
      }

      rho_1 = rho;
    }

    tol = resid;
    return 1;
  }

  template <typename T>
  inline int SparseMatrix<T>::solve_QMR( const Vector<T>& b, Vector<T>& x,
                                         int& max_iter, double& tol )
  {
    double resid;

    T rho( 1.0 ), rho_1( 1.0 ), xi( 1.0 );
    T delta( 1.0 ), ep( 1.0 );
    T beta;

    double gamma, gamma_1, eta, theta, theta_1;

    Vector<T> r, v, y, w, z;
    Vector<T> v_tld, w_tld, y_tld, z_tld;
    Vector<T> p, q, p_tld, d, s;

    double normb = b.norm_2();
    r = this->multiply( x );
    r = b - r;

    if ( normb == 0.0 )
    {
      normb = 1.0;
    }
    if ( ( resid = r.norm_2() / normb ) <= tol )
    {
      tol = resid;
      max_iter = 0;
      return 0;
    }

    v_tld = r;
    y = v_tld; // Could have preconditioner here
    rho = y.norm_2();

    w_tld = r;
    z = w_tld; // Could have preconditioner here
    xi = z.norm_2();

    gamma = 1.0;
    eta = - 1.0;
    theta = 0.0;

    for ( int i = 1; i <= max_iter; i++ )
    {
      if ( rho == 0.0 )
        return 2;                     // Return on breakdown

      if ( xi == 0.0 )
        return 7;                     // Return on breakdown

      v = ( 1. / rho ) * v_tld;
      y = ( 1. / rho ) * y;

      w = ( 1. / xi ) * w_tld;
      z = ( 1. / xi ) * z;

      delta = z.dot( y );
      if ( delta == 0.0 )
        return 5;                     // Return on breakdown

      y_tld = y;
      z_tld = z;  // Could have preconditioners here

      if (i > 1) {
        p = y_tld - ( xi * delta / ep ) * p;
        q = z_tld - ( rho * delta / ep ) * q;
      } else {
        p = y_tld;
        q = z_tld;
      }

      p_tld = this->multiply( p );
      ep = q.dot( p_tld );
      if ( ep == 0.0 )
        return 6;                       // Return on breakdown

      beta = ep / delta;
      if ( beta == 0.0 )
        return 3;                       // Return on breakdown

      v_tld = p_tld - beta * v;
      y = v_tld; // Could have preconditioner here

      rho_1 = rho;
      rho = y.norm_2();
      w_tld = this->transpose_multiply( q );
      w_tld -= beta * w;
      z = w_tld; // Could have preconditioner here

      xi = z.norm_2();

      gamma_1 = gamma;
      theta_1 = theta;

      theta = rho / ( gamma_1 * beta );
      gamma = 1.0 / std::sqrt( 1.0 + theta * theta );

      if ( gamma == 0.0 )
        return 4;                       // Return on breakdown

      eta = - eta * rho_1 * gamma * gamma / ( beta * gamma_1 * gamma_1 );

      if (i > 1) {
        d = eta * p + ( theta_1 * theta_1 * gamma * gamma ) * d;
        s = eta * p_tld + ( theta_1 * theta_1 * gamma * gamma ) * s;
      } else {
        d = eta * p;
        s = eta * p_tld;
      }

      x += d;
      r -= s;

      if ( ( resid = r.norm_2() / normb ) <= tol )
      {
        tol = resid;
        max_iter = i;
        return 0;
      }

    }

    tol = resid;
    return 1;
  }

}  // End of namespace Luna

#endif
