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
			T operator()( const T& x, const std::size_t& n );

			/// Evaluation operator for the nth Chebyshev polynomial at a Vector
			/// of points
			/// \param x A Vector of points where the Chebyshev polynomial is
			/// evaluated
			/// \param n The degree of the Chebyshev polynomial
			/// \return A Vector of values of the nth degree Chebyshev polynomial
			Vector<T> operator()( const Vector<T>& x, const std::size_t& n );

			/// Evaluation operator for the dth derivative of the nth Chebyshev
			///polynomial at point x
			/// \param x The point where the Chebyshev polynomial is evaluated
			/// \param n The degree of the Chebyshev polynomial
			/// \param d The derivative to return
			/// \return The dth derivative of the nth degree Chebyshev polynomial at
			/// the point x
			T operator()( const T& x, const std::size_t& n, const std::size_t& d );

			/// Evaluation operator for the dth derivative of the nth Chebyshev
			/// polynomial at a Vector of points
			/// \param x A Vector of points where the Chebyshev polynomial is
			/// evaluated
			/// \param n The degree of the Chebyshev polynomial
			/// \param d The derivative to return
			/// \return A Vector of the dth derivative of the nth degree Chebyshev
			/// polynomial at the Vector of points
			Vector<T> operator()( const Vector<T>& x, const std::size_t& n,
														const std::size_t& d );

      /* ----- Methods ----- */

			/// Return the nth Chebyshev polynomial of the first kind at point x
			/// \param x The point where the Chebyshev polynomial is evaluated
			/// \param n The degree of the Chebyshev polynomial
			/// \return The value of the nth degree Chebyshev polynomial at x
			T first_kind( const T& x, const std::size_t& n );

			/// Return the nth Chebyshev polynomial at a Vector of points
			/// \param x A Vector of points where the Chebyshev polynomial is
			/// evaluated
			/// \param n The degree of the Chebyshev polynomial
			/// \return A Vector of values of the nth degree Chebyshev polynomial
			Vector<T> first_kind( const Vector<T>& x, const std::size_t& n );

			/// Return the Gegenbauer or ultraspherical polynomial at a point x
			/// \param x The point where the Gegenbauer polynomial is evaluated
			/// \param n The degree of the Gegenbauer polynomial
			/// \param m The m parameter in the Gegenbauer polynomial
			/// \return The value of the Gegenbauer polynomial at x
			T gegenbauer( const T& x, const std::size_t& n, const std::size_t& m );


  }; // End of class Chebyshev

	template <typename T>
	Chebyshev<T>::Chebyshev()
	{}

	template <typename T>
	Chebyshev<T>::~Chebyshev()
	{}

	template <typename T>
	inline T Chebyshev<T>::operator()( const T& x, const std::size_t& n )
	{
		return first_kind( x, n );
	}

	template <typename T>
	inline Vector<T> Chebyshev<T>::operator()( const Vector<T>& x,
																						 const std::size_t& n )
	{
		return first_kind( x, n );
	}

	template <typename T>
	inline T Chebyshev<T>::operator()( const T& x, const std::size_t& n,
											 							 const std::size_t& d )
	{
		return n * std::pow( 2.0, d - 1 ) * factorial( d - 1 )
						 * gegenbauer( x, n - d, d );
	}

	template <typename T>
	inline Vector<T> Chebyshev<T>::operator()( const Vector<T>& x,
																						 const std::size_t& n,
																						 const std::size_t& d )
	{
		Vector<T> temp( x.size() );
		for ( std::size_t i = 0; i != x.size(); i++) {
			temp[ i ] = n * std::pow( 2.0, d - 1 ) * factorial( d - 1 )
										* gegenbauer( x[ i ], n - d, d );
		}
		return temp;
	}

	/* ----- Methods ----- */

	template <typename T>
	inline T Chebyshev<T>::first_kind( const T& x, const std::size_t& n )
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
	inline Vector<T> Chebyshev<T>::first_kind( const Vector<T>& x,
																						 const std::size_t& n )
	{
		Vector<T> temp( x.size() );
		for ( std::size_t i = 0; i != x.size(); i++) {
			temp[ i ] = first_kind( x[ i ], n );
		}
		return temp;
	}

	template <typename T>
	inline T Chebyshev<T>::gegenbauer( const T& x, const std::size_t& n,
																		 const std::size_t& m )
	{
		if ( n == 0 ){
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


}  // End of namespace Luna

#endif
