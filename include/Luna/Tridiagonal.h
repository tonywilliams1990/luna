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

namespace Luna
{

  // Templated Tridiagonal matrix class
  template <class T>

  /// A Tridiagonal matrix class for use with double and std::complex<double>
  class Tridiagonal
  {
  private:
    Vector<T> a;                // Subdiagonal of the matrix
    Vector<T> b;                // Main diagonal of the matrix
    Vector<T> c;                // Superdiagonal of the matrix
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
    /// \param super The subdiagonal of the matrix as a Vector
    Tridiagonal( const Vector<T>& sub, const Vector<T>& main,
                 const Vector<T>& super );

    //TODO constructor just using elements rather than Vectors (like fill_tridiag)

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



    //TODO more operator overloading

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
    /// \param super The subdiagonal of the matrix as a Vector
    void set_super( const Vector<T>& super );

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

    //TODO transpose_in_place, transpose, conjugate, resize, swap_rows, random
    //TODO various fill methods -> sub, main, super 

    /* ----- Norms ----- */


    /* ----- Solve linear systems ----- */

    //TODO solve tridiagonal system (Vector and Matrix inputs)

    /* ----- Determinant ----- */


    /* ----- Inverse ----- */


  };	// End of class Tridiagonal

  template <typename T>
  inline Tridiagonal<T>::Tridiagonal( const std::size_t& size )
  {
    a.resize( size );
    b.resize( size );
    c.resize( size );
    N = size;
  }

  template <typename T>
  inline Tridiagonal<T>::Tridiagonal( const Vector<T>& sub,
                                      const Vector<T>& main,
                                      const Vector<T>& super )
  {
    a = sub;
    b = main;
    c = super;
    N = main.size();
  }

  /* ----- Operator overloading ----- */

  template <class Type>
  inline std::ostream& operator<<( std::ostream& os,
                                   const Tridiagonal<Type>& m )
  {
    os << std::endl;
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
    return os;
  }

  template <typename T>
  inline const T& Tridiagonal<T>::operator() ( const std::size_t& i,
                                               const std::size_t& j ) const
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

}  // End of namespace Luna

#endif
