/// \file RationalSemi.h
/// A class specifying Rational Chebyshev basis functions on a semi-infinite
/// interval to be used in spectral methods. TODO \todo more of a definition of
/// rational chebyshev functions.


#ifndef RATIONALSEMI_H
#define RATIONALSEMI_H

#include <iostream>
#include <string>

#include "Matrix.h"
#include "Vector.h"
#include "Error.h"
#include "Basis.h"
#include "Special.h"
#include "Chebyshev.h"

namespace Luna
{
  /// A templated class for Rational Chebyshev basis functions TL_n(y)
	template <class T>

	class RationalSemi : public Basis<T>
	{
		private:
			double L; 						// Map parameter

    public:

      /// Constructor
			/// \param l The map parameter
			RationalSemi( const double& l );

      /// Destructor ( virtual since we have virtual functions )
      virtual ~RationalSemi();

			/* ----- Operator overloading ----- */

			/// Evaluation operator for the nth Semi-infinite Rational Chebyshev
			/// polynomial at point y
			/// \param y The point where the polynomial is evaluated
			/// \param n The degree of the polynomial
			/// \return The value of the nth degree polynomial at y
			T operator()( const T& y, const int& n );

			/// Evaluation operator for the nth Semi-infinite Rational Chebyshev
			/// polynomial at a Vector of points
			/// \param y A Vector of points where the polynomial is evaluated
			/// \param n The degree of the polynomial
			/// \return A Vector of values of the nth degree polynomial
			Vector<T> operator()( const Vector<T>& y, const int& n );

			/// Evaluation operator for the dth derivative of the nth Semi-infinite
			/// Rational Chebyshev polynomial at point y
			/// \param y The point where the polynomial is evaluated
			/// \param n The degree of the polynomial
			/// \param d The derivative to return
			/// \return The dth derivative of the nth degree polynomial at the point y
			T operator()( const T& y, const int& n, const int& d );

			/// Evaluation operator for the dth derivative of the nth Semi-infinite
			/// Rational Chebyshev polynomial at a Vector of points
			/// \param y A Vector of points where the polynomial is evaluated
			/// \param n The degree of the polynomial
			/// \param d The derivative to return
			/// \return A Vector of the dth derivative of the nth degree polynomial at
			/// the Vector of points
			Vector<T> operator()( const Vector<T>& y, const int& n,
														const int& d );

			/// Copy assignment
	    /// \param orig Semi-infinite Rational Chebyshev polynomial to be copied
	    /// \return The new Semi-infinite Rational Chebyshev polynomial
	    RationalSemi<T>& operator=( const RationalSemi<T>& orig );

      /* ----- Methods ----- */

			/// Return the nth Semi-infinite Rational Chebyshev polynomial at point y
			/// \param y The point where the polynomial is evaluated
			/// \param n The degree of the polynomial
			/// \return The value of the nth degree polynomial at y
			T evaluate( const T& y, const int& n );

			/// Return the nth Semi-infinite Rational Chebyshev polynomial at a Vector
			/// of points
			/// \param y A Vector of points where the polynomial is evaluated
			/// \param n The degree of the polynomial
			/// \return A Vector of values of the nth degree polynomial
			Vector<T> evaluate( const Vector<T>& y, const int& n );

			/// Evaluation of the dth derivative of the nth Semi-infinite Rational
			/// Chebyshev polynomial at point y
			/// \param y The point where the polynomial is evaluated
			/// \param n The degree of the polynomial
			/// \param d The derivative to return
			/// \return The dth derivative of the nth degree polynomial at the point y
			T derivative( const T& y, const int& n, const int& d );

			/// Evaluation of the dth derivative of the nth Semi-infinite Rational
			/// Chebyshev polynomial at point y
			/// \param y A Vector of points where the polynomial is evaluated
			/// \param n The degree of the polynomial
			/// \param d The derivative to return
			/// \return The dth derivative of the nth degree polynomial at the point y
			Vector<T> derivative( const Vector<T>& y, const int& n, const int& d );

			//TODO \todo approximate function


  }; // End of class Chebyshev

	template <typename T>
	RationalSemi<T>::RationalSemi( const double& l )
	{
		L = l;
		this->set_identifier( "RationalSemi" );
	}

	template <typename T>
	RationalSemi<T>::~RationalSemi()
	{}

	template <typename T>
	inline T RationalSemi<T>::operator()( const T& y, const int& n )
	{
		return evaluate( y, n );
	}

	template <typename T>
	inline Vector<T> RationalSemi<T>::operator()( const Vector<T>& y,
																								const int& n )
	{
		return evaluate( y, n );
	}

	template <typename T>
	inline T RationalSemi<T>::operator()( const T& y, const int& n, const int& d )
	{
		return derivative( y, n, d );
	}

	template <typename T>
	inline Vector<T> RationalSemi<T>::operator()( const Vector<T>& y,
																								const int& n, const int& d )
	{
		return derivative( y, n, d );
	}

	template <typename T>
  inline RationalSemi<T>& RationalSemi<T>::operator=(
																								const RationalSemi<T>& orig )
  {
		if ( this == &orig ){ return *this; }
		this->set_identifier( "RationalSemi" );
    return *this;
  }

	/* ----- Methods ----- */

	template <typename T>
	inline T RationalSemi<T>::evaluate( const T& y, const int& n )
	{
		Chebyshev<T> cheby;
		T x( ( y - L ) / ( y + L ) );
		return cheby( x, n );

		/*T t;
		t = 2 * std::atan( std::sqrt( L / y ) );
		return std::cos( n * t );*/
	}

	template <typename T>
	inline Vector<T> RationalSemi<T>::evaluate( const Vector<T>& y, const int& n )
	{
		Vector<T> temp( y.size() );
		for ( std::size_t i = 0; i != y.size(); i++) {
			temp[ i ] = evaluate( y[ i ], n );
		}
		return temp;
	}

	template <typename T>
	inline T RationalSemi<T>::derivative( const T& y, const int& n, const int& d )
	{
		if ( d == 1 )
		{
			Chebyshev<T> cheby;
			T x( ( y - L ) / ( y + L ) );
			T dxdy( 2 * L / std::pow( y + L, 2 ) );
			return dxdy * cheby( x, n, 1 );

			/*T t, s, c;
			t = 2 * std::atan( std::sqrt( L / y ) );
			s = std::sin( t / 2 );
			c = std::cos( t / 2 );
			return n * s * s * s * std::sin( n * t ) / ( L * c );*/
		}
		if ( d == 2 )
		{
			Chebyshev<T> cheby;
			T x( ( y - L ) / ( y + L ) );
			T dxdy( 2 * L / std::pow( y + L, 2 ) );
			T d2xdy2( - 4 * L / std::pow( y + L, 3 ) );
			return dxdy * dxdy * cheby( x, n, 2 ) + d2xdy2 * cheby( x, n, 1 );

			/*T t, s, c;
		 	t = 2 * std::atan( std::sqrt( L / y ) );
		 	s = std::sin( t / 2 );
		 	c = std::cos( t / 2 );
			return std::pow( s, 5 ) * ( 2 * c * s * ( - n * n * std::cos( n * t ) )
						 + ( 3 - 2 * s * s ) * ( - n * std::sin( n * t ) ) )
						 / ( 2 * L * L * std::pow( c, 3 ) );*/
		}
		if ( d == 3 )
		{
			Chebyshev<T> cheby;
			T x( ( y - L ) / ( y + L ) );
			T dxdy( 2 * L / std::pow( y + L, 2 ) );
			T d2xdy2( - 4 * L / std::pow( y + L, 3 ) );
			T d3xdy3( 12 * L / std::pow( y + L, 4 ) );
			return std::pow( dxdy, 3 ) * cheby( x, n, 3 ) +
						 3 * dxdy * d2xdy2 * cheby( x, n, 2 ) + d3xdy3 * cheby( x, n, 1 );
		}
		if ( d == 4 )
		{
			Chebyshev<T> cheby;
			T x( ( y - L ) / ( y + L ) );
			T dxdy( 2 * L / std::pow( y + L, 2 ) );
			T d2xdy2( - 4 * L / std::pow( y + L, 3 ) );
			T d3xdy3( 12 * L / std::pow( y + L, 4 ) );
			T d4xdy4( - 48 * L / std::pow( y + L, 5 ) );
			return std::pow( dxdy, 4 ) * cheby( x, n, 4 ) +
						 6 * dxdy * dxdy * d2xdy2 * cheby( x, n, 3 ) +
						 ( 4 * dxdy * d3xdy3 + 3 * d2xdy2 * d2xdy2 ) * cheby( x, n, 2 ) +
						 d4xdy4 * cheby( x, n, 1 );
		}
		else {
			std::string problem;
			problem += "RationalSemi: Only upto the 4th derivative has been";
			problem += "implemented so far.";
			throw Error( problem );
		}
	}

	template <typename T>
	inline Vector<T> RationalSemi<T>::derivative( const Vector<T>& y,
																								const int& n, const int& d )
	{
		Vector<T> temp( y.size() );
		for ( std::size_t i = 0; i != y.size(); i++) {
			temp[ i ] = derivative( y[ i ], n, d );
		}
		return temp;
	}

}  // End of namespace Luna

#endif