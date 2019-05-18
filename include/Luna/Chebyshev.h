/// \file Chebyshev.h
/// A class specifying Chebyshev basis functions to be used in spectral methods.


#ifndef CHEBYSHEV_H
#define CHEBYSHEV_H

#include <iostream>
#include <string>

#include "Matrix.h"
#include "Vector.h"
#include "Error.h"
#include "Basis.h"
#include "Special.h"

namespace Luna
{
  /// A templated class for Chebyshev basis functions
	template <class T>

	class Chebyshev : public Basis<T>
	{
		private:

    public:

      /// Constructor
			Chebyshev();

      /// Destructor ( virtual since we have virtual functions )
      virtual ~Chebyshev();

			/* ----- Operator overloading ----- */

			/// Evaluation operator for the nth Chebyshev polynomial at point x
			/// \param x The point where the Chebyshev polynomial is evaluated
			/// \param n The degree of the Chebyshev polynomial
			/// \return The value of the nth degree Chebyshev polynomial at x
			T operator()( const T& x, const int& n );

			/// Evaluation operator for the nth Chebyshev polynomial at a Vector
			/// of points
			/// \param x A Vector of points where the Chebyshev polynomial is
			/// evaluated
			/// \param n The degree of the Chebyshev polynomial
			/// \return A Vector of values of the nth degree Chebyshev polynomial
			Vector<T> operator()( const Vector<T>& x, const int& n );

			/// Evaluation operator for the dth derivative of the nth Chebyshev
			/// polynomial at point x
			/// \param x The point where the Chebyshev polynomial is evaluated
			/// \param n The degree of the Chebyshev polynomial
			/// \param d The derivative to return
			/// \return The dth derivative of the nth degree Chebyshev polynomial at
			/// the point x
			T operator()( const T& x, const int& n, const int& d );

			/// Evaluation operator for the dth derivative of the nth Chebyshev
			/// polynomial at a Vector of points
			/// \param x A Vector of points where the Chebyshev polynomial is
			/// evaluated
			/// \param n The degree of the Chebyshev polynomial
			/// \param d The derivative to return
			/// \return A Vector of the dth derivative of the nth degree Chebyshev
			/// polynomial at the Vector of points
			Vector<T> operator()( const Vector<T>& x, const int& n,
														const int& d );

			/// Copy assignment
	    /// \param original Chebyshev polynomial to be copied
	    /// \return The new Chebyshev polynomial
	    Chebyshev<T>& operator=( const Chebyshev<T>& original );

      /* ----- Methods ----- */

			/// Return the nth Chebyshev polynomial of the first kind at point x
			/// \param x The point where the Chebyshev polynomial is evaluated
			/// \param n The degree of the Chebyshev polynomial
			/// \return The value of the nth degree Chebyshev polynomial at x
			T first_kind( const T& x, const int& n );

			/// Return the nth Chebyshev polynomial at a Vector of points
			/// \param x A Vector of points where the Chebyshev polynomial is
			/// evaluated
			/// \param n The degree of the Chebyshev polynomial
			/// \return A Vector of values of the nth degree Chebyshev polynomial
			Vector<T> first_kind( const Vector<T>& x, const int& n );

			/// Return the Gegenbauer or ultraspherical polynomial at a point x
			/// \param x The point where the Gegenbauer polynomial is evaluated
			/// \param n The degree of the Gegenbauer polynomial
			/// \param m The m parameter in the Gegenbauer polynomial
			/// \return The value of the Gegenbauer polynomial at x
			T gegenbauer( const T& x, const int& n, const int& m );

			/// Approximate a given function by a Chebyshev polynomial in the
			/// interval [-1, 1]
			/// \param func The function to be approximated
			/// \param N The number of Chebyshev coefficients to calculate
			/// \return A Vector containing the Chebyshev coefficients
			Vector<T> approximate( T func( const T& ), const int& N );

			//TODO fast cosine transform to improve speed ??? NR

			/// Approximate a given function by an even Chebyshev polynomial in the
			/// interval [-1, 1]
			/// \param func The function to be approximated
			/// \param N The number of even Chebyshev coefficients to calculate
			/// \return A Vector containing the even Chebyshev coefficients
			Vector<T> approximate_even( T func( const T& ), const int& N );

			/// Approximate a given function by an odd Chebyshev polynomial in the
			/// interval [-1, 1]
			/// \param func The function to be approximated
			/// \param N The number of odd Chebyshev coefficients to calculate
			/// \return A Vector containing the odd Chebyshev coefficients
			Vector<T> approximate_odd( T func( const T& ), const int& N );


  }; // End of class Chebyshev

	template <typename T>
	Chebyshev<T>::Chebyshev()
	{
		this->set_identifier( "Chebyshev" );
	}

	template <typename T>
	Chebyshev<T>::~Chebyshev()
	{}

	template <typename T>
	inline T Chebyshev<T>::operator()( const T& x, const int& n )
	{
		return first_kind( x, n );
	}

	template <typename T>
	inline Vector<T> Chebyshev<T>::operator()( const Vector<T>& x, const int& n )
	{
		return first_kind( x, n );
	}

	template <typename T>
	inline T Chebyshev<T>::operator()( const T& x, const int& n, const int& d )
	{
		if ( n < d ){
			return 0;
		} else {
			return n * std::pow( 2.0, d - 1 ) * factorial( d - 1 )
						 	 * gegenbauer( x, n - d, d );
		}
	}

	template <typename T>
	inline Vector<T> Chebyshev<T>::operator()( const Vector<T>& x, const int& n,
																						 const int& d )
	{
		Vector<T> temp( x.size(), 0 );
		if ( n < d ){
			return temp;
		} else {
			for ( std::size_t i = 0; i != x.size(); i++) {
				temp[ i ] = n * std::pow( 2.0, d - 1 ) * factorial( d - 1 )
											* gegenbauer( x[ i ], n - d, d );
			}
			return temp;
		}
	}

	template <typename T>
  inline Chebyshev<T>& Chebyshev<T>::operator=( const Chebyshev<T>& original )
  {
		if ( this == &original ){ return *this; }
		this->set_identifier( "Chebyshev" );
    return *this;
  }

	/* ----- Methods ----- */

	template <typename T>
	inline T Chebyshev<T>::first_kind( const T& x, const int& n )
	{
		if ( n == 0 ){
			return 1;
		}
		if ( n == 1 ){
			return x;
		}
		else {
			return 2 * x * first_kind( x, n - 1 ) - first_kind( x, n - 2 );
		}
	}

	template <typename T>
	inline Vector<T> Chebyshev<T>::first_kind( const Vector<T>& x, const int& n )
	{
		Vector<T> temp( x.size() );
		for ( std::size_t i = 0; i != x.size(); i++) {
			temp[ i ] = first_kind( x[ i ], n );
		}
		return temp;
	}

	template <typename T>
	inline T Chebyshev<T>::gegenbauer( const T& x, const int& n,
																		 const int& m )
	{
		if ( n <= 0 ){
			return 1;
		}
		if ( n == 1 ){
			return 2 * m * x;
		}
		else {
			return ( 2 * x * ( n + m - 1 ) * gegenbauer( x, n - 1, m )
							- ( n + 2 * m - 2 ) * gegenbauer( x, n - 2, m ) )
							/ ( 1.0 * n );
		}
	}

	template <typename T>
	inline Vector<T> Chebyshev<T>::approximate( T func( const T& ), const int& N )
	{
		Vector<T> c( N );
		Vector<double> x;
		x.chebyshev_grid( N );
		for ( std::size_t j = 0; j < N; j++ )
		{
			T sum( 0.0 );
			for ( std::size_t k = 0; k < N; k++ )
			{
				sum += func( x[ k ] ) * std::cos( ( j * M_PI * ( k + 0.5 ) ) / N );
			}
			c[ j ] = ( 2.0 / N ) * sum;
		}
		return c;
	}

	template <typename T>
	inline Vector<T> Chebyshev<T>::approximate_even( T func( const T& ),
																									 const int& N )
	{
		Vector<T> c( N );
		Vector<double> x;
		x.chebyshev_grid( 2 * N );
		for ( std::size_t j = 0; j < N; j++ )
		{
			T sum( 0.0 );
			for ( std::size_t k = 0; k < 2 * N; k++ )
			{
				sum += func( x[ k ] ) * std::cos( ( 2 * j * M_PI * ( k + 0.5 ) )
																				/ ( 2. * N ) );
			}
			c[ j ] = ( 2.0 / ( 2. * N ) ) * sum;
		}
		return c;
	}

	template <typename T>
	inline Vector<T> Chebyshev<T>::approximate_odd( T func( const T& ),
																									const int& N )
	{
		Vector<T> c( N );
		Vector<double> x;
		x.chebyshev_grid( 2 * N );
		for ( std::size_t j = 0; j < N; j++ )
		{
			T sum( 0.0 );
			for ( std::size_t k = 0; k < 2 * N; k++ )
			{
				sum += func( x[ k ] ) * std::cos( ( ( 2 * j + 1 ) * M_PI *
																					( k + 0.5 ) ) / ( 2. * N ) );
			}
			c[ j ] = ( 2.0 / ( 2. * N ) ) * sum;
		}
		return c;
	}


}  // End of namespace Luna

#endif