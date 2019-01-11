/// \file Polynomial.h


#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <iostream>
#include <string>
#include <cmath>

#include "Vector.h"
#include "Error.h"

namespace Luna
{
  /// A templated base class to be inherited by objects that define residuals
	template <class T>

	class Polynomial
	{
		private:
      Vector<T> COEFFS;	// Vector for storing the coefficients of the polynomial
			std::size_t N; 		// The degree of the polynomial

    public:

			/// Constructor for a unspecified polynomial
			Polynomial(){}

			/// Constructor for a polynomial with specified coefficients
			Polynomial( const Vector<T>& coeffs ) : COEFFS( coeffs )
			{
				N = coeffs.size() - 1;
			}

			/// Constructor for a polynomial of given coefficients (using std::vector)
			Polynomial( const std::vector<T>& coeffs )
			{
				COEFFS.VECTOR = coeffs;
				N = coeffs.size() - 1;
			}

			/// Destructor
			~Polynomial(){}

			/* ----- Operator overloading ----- */

			/// Indexing operator ( read only )
	    /// \param i Index
	    /// \return A constant reference to the coefficient at the given index
	    const T& operator[] ( const std::size_t& i ) const;

			/// Indexing operator
	    /// \param i Index
	    /// \return A constant reference to the coefficient at the given index
	    T& operator[] ( const std::size_t& i );

			/// Evaluation operator
			/// \param x Arguement of the Polynomial function
			/// \return Value of the Polynomial function
			T operator() ( const T& x );

			/* ----- Methods ----- */

			/// Evaluation of the Polynomial
			/// \param x Arguement of the Polynomial function
			/// \return Value of the Polynomial function
			T evaluate( const T& x );

			/// Quadratic initialisation \f$ ax^2 + bx + c \f$
			/// \param a Coefficient of \f$ x^2 \f$
			/// \param b Coefficient of \f$ x \f$
			/// \param c Constant value in the quadratic
			void quadratic( const T& a, const T& b, const T& c );

			/// Solve the quadratic equation \f$ ax^2 + bx + c = 0 \f$
			/// \param a Coefficient of \f$ x^2 \f$
			/// \param b Coefficient of \f$ x \f$
			/// \param c Constant value in the quadratic
			/// \return A Vector< std::complex<double> > containing the two roots
			Vector< std::complex<double> > quadratic_solve( const T& a, const T& b,
																											const T& c );

			/// Cubic initialisation \f$ ax^3 + bx^2 + cx + d \f$
			/// \param a Coefficient of \f$ x^3 \f$
			/// \param b Coefficient of \f$ x^2 \f$
			/// \param c Coefficient of \f$ x \f$
			/// \param d Constant value in the quadratic
			void cubic( const T& a, const T& b, const T& c, const T& d );

			/// Solve the cubic equation \f$ ax^3 + bx^2 + cx + d = 0 \f$
			/// \param a Coefficient of \f$ x^3 \f$
			/// \param b Coefficient of \f$ x^2 \f$
			/// \param c Coefficient of \f$ x \f$
			/// \param d Constant value in the quadratic
			/// \return A Vector< std::complex<double> > containing the three roots
			Vector< std::complex<double> > cubic_solve( const T& a, const T& b,
																									const T& c, const T& d );

			//TODO solve method with specialisation for quadratic and cubic (and linear/constant N=1/0)


			//TODO polynomial division

			//TODO rational functions

			//TODO derivatives of the polynomial at a given point




  }; // End of class Polynomial

	template <typename T>
  inline const T& Polynomial<T>::operator[]( const std::size_t& i ) const
  {
    if ( i < 0 || N < i )
    {
      throw Error( "Polynomial operator [] range error" );
    }
    return COEFFS[ i ];
  }

	template <typename T>
  inline T& Polynomial<T>::operator[]( const std::size_t& i )
  {
    if ( i < 0 || N < i )
    {
      throw Error( "Polynomial operator [] range error" );
    }
    return COEFFS[ i ];
  }

	template <typename T>
	inline T Polynomial<T>::operator() ( const T& x )
	{
		return evaluate( x );
	}

	template <typename T>
	inline T Polynomial<T>::evaluate( const T& x )
	{
		T p( COEFFS[ N ] );
		for ( std::size_t j = 0; j < N; j++ )
		{
			p = p * x + COEFFS[ N - j - 1 ];
		}
		return p;
	}

	template <typename T>
	inline void Polynomial<T>::quadratic( const T& a, const T& b, const T& c )
	{
		N = 2;
		COEFFS.resize( 3 );
		COEFFS[ 0 ] = c;
		COEFFS[ 1 ] = b;
		COEFFS[ 2 ] = a;
	}

	template <typename T>
	inline Vector< std::complex<double> > Polynomial<T>::quadratic_solve(
																						const T& a, const T& b, const T& c )
	{
		this->quadratic( a, b, c );
		Vector< std::complex<double> > result( 2 );

		std::complex<double> disc( b * b - 4. * a * c );
		double sgn;
		sgn = std::real( std::conj( b ) * std::sqrt( disc ) );

		if ( sgn >= 0.0 ){
			sgn = 1.0;
		} else {
			sgn = - 1.0;
		}

		std::complex<double> q;
		q = - 0.5 * ( b + sgn * std::sqrt( disc ) );
		result[ 0 ] = q / a;
		result[ 1 ] = c / q;

		return result;
	}

	template <typename T>
	inline void Polynomial<T>::cubic( const T& a, const T& b, const T& c,
																		const T& d )
	{
		N = 3;
		COEFFS.resize( 4 );
		COEFFS[ 0 ] = d;
		COEFFS[ 1 ] = c;
		COEFFS[ 2 ] = b;
		COEFFS[ 3 ] = a;
	}

	template <typename T>
	inline Vector< std::complex<double> > Polynomial<T>::cubic_solve( const T& a,
																						const T& b, const T& c, const T& d )
	{
		this->cubic( a, b, c, d );
		Vector< std::complex<double> > result( 3 );

		T e, f, g;
		e = b / a;
		f = c / a;
		g = d / a;

		T Q, R;
		Q = ( e * e - 3. * f ) / 9.;
		R = ( 2. * e * e * e - 9. * e * f + 27. * g ) / 54.;

		std::complex<double> A, B;
		A = - std::pow( R + std::sqrt( R * R - Q * Q * Q ) , 1. / 3. );
		if ( A != 0. ) {
			B = Q / A;
		} else {
			B = 0.;
		}

		result[ 0 ] = A + B - ( e / 3. );
		std::complex<double> first, second;
		std::complex<double> i( 0.0, 1.0 );
		first  = - 0.5 * ( A + B ) - ( e / 3. );
		second = i * 0.5 * std::sqrt( 3. ) * ( A - B );
		result[ 1 ] = first + second;
		result[ 2 ] = first - second;

		return result;
	}

}  // End of namespace Luna

#endif
