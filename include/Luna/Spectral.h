/// \file Spectral.h
/// A class for defining Spectral solutions using specified basis functions and
/// coefficients. The spectral function \f$ u_N(x) \f$ is defined by
/// \f[ u_N(x) = \sum_{n=0}^{N} a_n \phi_n(x),  \f] where \f$ a_n \f$ are
/// coefficients and \f$ \phi_n(x) \f$ are known basis functions. 

#ifndef SPECTRAL_H
#define SPECTRAL_H

#include <iostream>
#include <string>
#include <cmath>

#include "Vector.h"
#include "Error.h"
#include "Chebyshev.h"

namespace Luna
{
  /// A templated class for spectral solutions
	template <class T>

	class Spectral
	{
		private:
			Vector<T> COEFFICIENTS; 	// Coefficients multiplying the basis functions
			std::string BASIS;				// Basis function specifier string

    public:

			/// Constructor ( default = no coefficients and Chebyshev basis )
			Spectral();

			/// Constructor with specified basis and coefficients
			/// \param coefficients The Vector of coefficients
			/// \param basis The basis specifier string (i.e. Chebyshev)
			Spectral( const Vector<T>& coefficients, const std::string& basis );

			/// Destructor
			~Spectral(){}

			/* ----- Operator overloading ----- */

			/// Evaluation operator at a point x
			/// \param x The point where the spectral function is evaluated
			T operator()( const T& x );

			/// Evaluation operator at a Vector of points
			/// \param x A Vector of points where the spectral function is evaluated
			Vector<T> operator()( const Vector<T>& x );

			//TODO derivatives operator, Coefficient access operator []

			/* ----- Methods ----- */

			/// Return the basis specifier string
			/// \return The string containing the basis specifier
			std::string get_basis();

			/// Set the basis specifier string
			/// \param basis The basis specifier string (i.e. Chebyshev)
			void set_basis( const std::string& basis );

			//TODO set coefficients using Vector, get Vector of coefficients


  }; // End of class Spectral

	template <typename T>
	Spectral<T>::Spectral()
	{
		BASIS = "Chebyshev";
	}

	template <typename T>
	Spectral<T>::Spectral( const Vector<T>& coefficients,
												 const std::string& basis )
	{
		if ( basis != "Chebyshev" )//&& basis != "Fourier" )
		{
			std::string problem;
			problem += "Spectral set_basis error: The basis function " + basis + "\n";
			problem += "is not supported.";
			throw Error( problem );
		}
		COEFFICIENTS = coefficients;
		BASIS = basis;
	}

	/* ----- Operator overloading ----- */

	template <typename T>
	inline T Spectral<T>::operator()( const T& x )
	{
		T temp( 0 );
		std::size_t N( COEFFICIENTS.size() );

		if ( BASIS == "Chebyshev" )
		{
			Chebyshev<T> basis;
			for ( std::size_t n = 0; n < N; ++n )
			{
				temp += basis( x, n ) * COEFFICIENTS[ n ];
			}
		}

		return temp;
	}

	template <typename T>
	inline Vector<T> Spectral<T>::operator()( const Vector<T>& x )
	{
		Vector<T> temp( x.size() );
		std::size_t N( COEFFICIENTS.size() );

		if ( BASIS == "Chebyshev" )
		{
			Chebyshev<T> basis;
			for ( std::size_t n = 0; n < N; ++n )
			{
				temp += basis( x, n ) * COEFFICIENTS[ n ];
			}
		}

		return temp;
	}

	/* ----- Methods ----- */

	template <typename T>
	std::string Spectral<T>::get_basis()
	{
		return BASIS;
	}

	template <typename T>
	void Spectral<T>::set_basis( const std::string& basis )
	{
		if ( basis != "Chebyshev" )//&& basis != "Fourier" )
		{
			std::string problem;
			problem += "Spectral set_basis error: The basis function " + basis + "\n";
			problem += "is not supported.";
			throw Error( problem );
		}
		BASIS = basis;
	}

}  // End of namespace Luna

#endif
