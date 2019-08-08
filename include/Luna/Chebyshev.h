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
			// Derivatives of the Chebyshev polynomials of the first kind
			T first_kind_derivative( const T& x, const int& n );
			T first_kind_second_derivative( const T& x, const int& n );
			T first_kind_third_derivative( const T& x, const int& n );

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

			/// Return the nth Chebyshev polynomial of the second kind at point x
			/// \param x The point where the Chebyshev polynomial is evaluated
			/// \param n The degree of the Chebyshev polynomial
			/// \return The value of the nth degree Chebyshev polynomial at x
			T second_kind( const T& x, const int& n );

			/// Return the nth Chebyshev polynomial of the second kind at a Vector of
			/// points
			/// \param x A Vector of points where the Chebyshev polynomial is
			/// evaluated
			/// \param n The degree of the Chebyshev polynomial
			/// \return A Vector of values of the nth degree Chebyshev polynomial
			Vector<T> second_kind( const Vector<T>& x, const int& n );

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
			Vector<T> approximate( T func ( const T& ), const int& N );

			/// \todo TODO fast cosine transform to improve speed ??? NR

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

			/// Return a Vector containing the the nth Chebyshev polynomial at point x
			/// and its first two derivatives
			/// \param x The point where the polynomial is evaluated
			/// \param n The degree of the polynomial
			/// \return A Vector containing the nth degree polynomial at x, the first
			/// derivative and the second derivative.
			Vector<T> eval_2( const T& x, const int& n );


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
		} if ( d == 0 ) {
			return first_kind( x, n );
		}	if ( d == 1 ) {
			return first_kind_derivative( x, n );
		} if ( d == 2 ) {
				/*if ( x == 1.0 ) {
					return 1. * n * n * ( n * n - 1 ) / 3.0;
				}
				if ( x == -1.0 ) {
					int p;
					p = n % 2 ? -1 : 1;
					return 1. * p * n * n * ( n * n - 1 ) / 3.0;
				} else {
					return 1. * n * ( ( n + 1 ) * first_kind( x, n )
																			- second_kind( x, n ) ) / ( x * x - 1. );
				}*/
				return first_kind_second_derivative( x, n );
		} if ( d == 3 ) {
			return first_kind_third_derivative( x, n );
			/// \todo TODO d = 4,... (really speeds things up)
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
		} if ( d == 0 ) {
			return first_kind( x, n );
		} if ( d == 1 ) {
			for ( std::size_t i = 0; i != x.size(); i++) {
				temp[ i ] = first_kind_derivative( x[ i ], n );
			}
			return temp;
		} if ( d == 2 ) {
			for ( std::size_t i = 0; i != x.size(); i++) {
				temp[ i ] = first_kind_second_derivative( x[ i ], n );
			}
			return temp;
		} if ( d == 3 ) {
			for ( std::size_t i = 0; i != x.size(); i++) {
				temp[ i ] = first_kind_third_derivative( x[ i ], n );
			}
			return temp;
		/// \todo TODO d = 4,... (really speeds things up)
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
		if ( n == 2 ){
			return 2 * x * x - 1.0;
		}
		else {
			// Could use the recursuve definition
			// return 2 * x * first_kind( x, n - 1 ) - first_kind( x, n - 2 );
			// but this is a bit slow.

			T tnm1( 2 * x * x - 1.0 );
			T tnm2( x );
			T tn( tnm1 );

			for ( unsigned l = 3; l <= n; l++ )
			{
				tn = 2. * x * tnm1 - tnm2;
				tnm2 = tnm1;
				tnm1 = tn;
			}
			return tn;
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
	inline T Chebyshev<T>::second_kind( const T& x, const int& n )
	{
		if ( n == 0 ){
			return 1;
		}
		if ( n == 1 ){
			return 2 * x;
		}
		if ( n == 2 ){
			return 4 * x * x - 1.0;
		}
		else {
			T unm1( 4 * x * x - 1.0 );
			T unm2( 2 * x );
			T un( unm1 );

			for ( unsigned l = 3; l <= n; l++ )
			{
				un = 2. * x * unm1 - unm2;
				unm2 = unm1;
				unm1 = un;
			}
			return un;
		}
	}

	template <typename T>
	inline Vector<T> Chebyshev<T>::second_kind( const Vector<T>& x, const int& n )
	{
		Vector<T> temp( x.size() );
		for ( std::size_t i = 0; i != x.size(); i++) {
			temp[ i ] = second_kind( x[ i ], n );
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

	template <typename T>
	inline Vector<T> Chebyshev<T>::eval_2( const T& x, const int& n )
	{
		Vector<T> temp( 3 );

		if ( std::abs( x ) < 1. )
		{
			T t( std::acos( x ) );
			T c( std::cos( t ) );
			T s( std::sin( t ) );
			//T tn( std::cos( n * t ) );
			//T tnt( - n * std::sin( n * t ) );
			//T tntt( - n * n * tn );

			//T tnx( - tnt / s );
			//T tnxx( ( tntt / ( s * s ) ) + ( c * tnx / ( s * s ) ) );

			temp[ 0 ] = std::cos( n * t );

			//temp[ 1 ] = - tnt / s;
			temp[ 1 ]  = std::sin( n * t );
			temp[ 1 ] *= n;
			temp[ 1 ] /= s;

			//temp[ 2 ] = ( tntt / ( s * s ) ) - ( c * tnt / ( s * s * s ) );
			temp[ 2 ]  = - temp[ 0 ];
			temp[ 2 ] *= n * n;
			temp[ 2 ] += c * temp[ 1 ];
			temp[ 2 ] /= s * s;

		} else {
			/// \todo TODO could this be done more efficiently x = +/- 1?
			temp[ 0 ] = std::pow( x, n );
			temp[ 1 ] = std::pow( x, n + 1 ) * n * n;
			temp[ 2 ] = std::pow( x, n ) * n * n * ( n * n - 1. ) / 3.;
		}

		return temp;
	}

	/* ----- Private ----- */

	template <typename T>
	inline T Chebyshev<T>::first_kind_derivative( const T& x, const int& n )
	{
		return second_kind( x, n - 1 ) * n;
	}

	template <typename T>
	inline T Chebyshev<T>::first_kind_second_derivative( const T& x,
																											 const int& n )
	{
		if ( n == 0 ){
			return 0;
		}
		if ( n == 1 ){
			return 0;
		}
		if ( n == 2 ){
			return 4.;
		}
		if ( n == 3 ){
			return 24 * x;
		}
		if ( n == 4 ){
			return 96 * x * x - 16. ;
		}
		else {
			T tnm1( 96 * x * x - 16. );
			T tnm2( 24 * x );
			T tn( tnm1 );

			for ( unsigned l = 5; l <= n; l++ )
			{
				tn = 2. * x * tnm1 + 4 * first_kind_derivative( x, l - 1 ) - tnm2;
				tnm2 = tnm1;
				tnm1 = tn;
			}
			return tn;
		}
	}

	template <typename T>
	inline T Chebyshev<T>::first_kind_third_derivative( const T& x, const int& n )
	{
		if ( n == 0 ){
			return 0;
		}
		if ( n == 1 ){
			return 0;
		}
		if ( n == 2 ){
			return 0;
		}
		if ( n == 3 ){
			return 24.;
		}
		if ( n == 4 ){
			return 192 * x;
		}
		if ( n == 5 ){
			return 960 * x * x - 120.;
		}
		else {
			T tnm1( 960 * x * x - 120. );
			T tnm2( 192 * x );
			T tn( tnm1 );

			for ( unsigned l = 6; l <= n; l++ )
			{
				tn = 2. * x * tnm1 + 6 * first_kind_second_derivative( x, l - 1 )
						 - tnm2;
				tnm2 = tnm1;
				tnm1 = tn;
			}
			return tn;
		}
	}


}  // End of namespace Luna

#endif
