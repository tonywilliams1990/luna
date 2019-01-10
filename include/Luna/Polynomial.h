/// \file Polynomial.h


#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <iostream>
#include <string>

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

			//TODO quadratic initialisation ax^2 + bx + c -> c[0] + c[1]*x + c[2]*x^2

			//TODO quadratic solve
			// return complex values always  (in Vector)

			//TODO cubic initialisation

			//TODO cubic solve

			//TODO solve method with specialisation for quadratic and cubic (and linear/constant N=1/0)


			//TODO polynomial division

			//TODO rational functions

			//TODO derivatives of the polynomial at a given point




  }; // End of class Polynomial

	template <typename T>
  inline const T& Polynomial<T>::operator[]( const std::size_t& i ) const
  {
    if ( i < 0 || N <= i )
    {
      throw Error( "Polynomial operator [] range error" );
    }
    return COEFFS[ i ];
  }

	template <typename T>
  inline T& Polynomial<T>::operator[]( const std::size_t& i )
  {
    if ( i < 0 || N <= i )
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

}  // End of namespace Luna

#endif
