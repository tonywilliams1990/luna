/// \file Spectral2D.h
/// \todo TODO fix description for 2D
/// A class for defining Spectral solutions using specified basis functions and
/// coefficients. The spectral function \f$ u_N(x) \f$ is defined by
/// \f[ u_N(x) = \sum_{n=0}^{N} a_n \phi_n(x),  \f] where \f$ a_n \f$ are
/// coefficients and \f$ \phi_n(x) \f$ are known basis functions.

#ifndef SPECTRAL2D_H
#define SPECTRAL2D_H

#include <iostream>
#include <string>
#include <cmath>

#include "Vector.h"
#include "Error.h"
#include "Chebyshev.h"
#include "Special.h"
#include "RationalSemi.h"

namespace Luna
{
  /// A templated class for spectral solutions
	template <class T>

	class Spectral2D
	{
		private:
			Vector<T> COEFFICIENTS; 	// Coefficients multiplying the basis functions
			std::string BASIS;				// Basis function specifier string
			double PARAM;							// Parameter for rational Chebyshev basis
			int I;										// Number of x collocation points
			int J;										// Number of y collocation points

    public:

			/// Constructor ( default = no coefficients and Chebyshev basis )
			Spectral2D();

			/// Constructor with specified basis and coefficients
			/// \param coefficients The Vector of coefficients
			/// \param i The number of collocation points in the x-direction
			/// \param j The number of collocation points in the y-direction
			/// \param basis The basis specifier string (i.e. Chebyshev)
			Spectral2D( const Vector<T>& coefficients, const int& i, const int& j,
									const std::string& basis );

			/// Constructor with specified basis, coefficients and extra parameter
			/// \param coefficients The Vector of coefficients
			/// \param basis The basis specifier string (i.e. Chebyshev)
			/// \param param An extra parameter for use with rational Chebyshev basis
			Spectral2D( const Vector<T>& coefficients, const std::string& basis,
			 						const double& param );
			/// \todo TODO modify this constructor to include number of collocation points

			/// Destructor
			~Spectral2D(){}

			/* ----- Operator overloading ----- */

			/// Evaluation operator at a point (x,y)
			/// \param x The x-coordinate where the spectral function is evaluated
			/// \param y The y-coordinate where the spectral function is evaluated
			/// \return The value of the spectral function at (x,y)
			T operator()( const T& x, const T& y );

			/// Evaluation operator at a Vector of points
			/// \param x A Vector of points where the spectral function is evaluated
			/// \return A Vector of values of the spectral function
			Vector<T> operator()( const Vector<T>& x );
			/// \todo TODO change this to passing a Mesh2D object and a variable number

			/// Evaluation operator for the dth derivative of the spectral function
			/// at point x
			/// \param x The point where the spectral function is evaluated
			/// \param d The derivative to return
			/// \return The dth derivative of the spectral function at the point x
			T operator()( const T& x, const int& d );
			/// \todo TODO modify this for x and y derivatives

			/// Evaluation operator for the dth derivative of the spectral function
			/// at a Vector of points
			/// \param x A Vector of points where the spectral function is evaluated
			/// \param d The derivative to return
			/// \return A Vector of the dth derivative of the spectral function
			/// at the Vector of points
			Vector<T> operator()( const Vector<T>& x, const int& d );
			/// \todo TODO pass Mesh2D object and modify for x and y derivatives

			/// Indexing operator ( read only )
	    /// \param i Index
	    /// \return A constant reference to the coefficient located at the given
			/// index
	    const T& operator[] ( const std::size_t& i ) const;

	    /// Indexing operator ( read / write )
	    /// \param i Index
	    /// \return A reference to the coefficient located at the given index
	    T& operator[] ( const std::size_t& i );

			/* ----- Methods ----- */

			/// Return the basis specifier string
			/// \return The string containing the basis specifier
			std::string get_basis();

			/// Set the basis specifier string
			/// \param basis The basis specifier string (i.e. Chebyshev)
			void set_basis( const std::string& basis );

			/// Return the spectral coefficients
			/// \return The Vector of spectral coefficients
			Vector<T> get_coefficients();

			/// Set the spectral coefficients
			/// \param coeffs The Vector of coefficients
			void set_coefficients( const Vector<T>& coeffs );

			/// Update the spectral coefficients (add a Vector to the Vector of
			/// coefficients)
			/// \param update The Vector to be added to the Vector of coefficients
			void update_coefficients( const Vector<T>& update );

			/// Add a new coefficient to the back of the Vector of coefficients
			/// \param new_elem New element to be appended to the end of the Vector
			/// of coefficients
			void push_back( const T& new_elem );

			/// Add a new coefficient to the front of the Vector of coefficients
			/// \param new_elem New element to be added to the front of the Vector
			/// of coefficients
			void push_front( const T& new_elem );

			/// Return a pointer to the extra parameter
			/// \return A handle to the parameter
			double& parameter();

			/// Return a pointer to the number of x collocation points
			/// \return A handle to the number of x collocation points
			int& x_points();

			/// Return a pointer to the number of y collocation points
			/// \return A handle to the number of y collocation points
			int& y_points();

  }; // End of class Spectral2D

	template <typename T>
	Spectral2D<T>::Spectral2D()
	{
		BASIS = "Chebyshev";
	}

	template <typename T>
	Spectral2D<T>::Spectral2D( const Vector<T>& coefficients, const int& i,
														 const int& j, const std::string& basis )
	{
		if ( basis != "Chebyshev" && basis != "EvenChebyshev"
		  && basis != "OddChebyshev" && basis != "RationalSemi" )
			//&& basis != "Fourier" )
		{
			std::string problem;
			problem += "Spectral set_basis error: The basis function " + basis + "\n";
			problem += "is not supported.";
			throw Error( problem );
		}
		if ( coefficients.size() != i * j )
		{
			std::string problem;
			problem += "Spectral2D error: Number of coefficients and collocation \n";
			problem += "points do not agree.";
			throw Error( problem );
		}
		COEFFICIENTS = coefficients;
		BASIS = basis;
		PARAM = 1.0;
		I = i;
		J = j;
	}

	template <typename T>
	Spectral2D<T>::Spectral2D( const Vector<T>& coefficients,
												 const std::string& basis, const double& param )
	{
		if ( basis != "Chebyshev" && basis != "EvenChebyshev"
		  && basis != "OddChebyshev" && basis != "RationalSemi" )
			//&& basis != "Fourier" )
		{
			std::string problem;
			problem += "Spectral set_basis error: The basis function " + basis + "\n";
			problem += "is not supported.";
			throw Error( problem );
		}
		COEFFICIENTS = coefficients;
		BASIS = basis;
		PARAM = param;
	}

	/* ----- Operator overloading ----- */

	template <typename T>
	inline T Spectral2D<T>::operator()( const T& x, const T& y )
	{
		T temp( 0 );
		int size( I * J );
		int f, g;

		if ( BASIS == "Chebyshev" )
		{
			Chebyshev<T> basis;
			for ( std::size_t n = 0; n < size; n++ )
      {
        f = n / J;
        g = n % J;
        temp += COEFFICIENTS[ n ] * basis( x, f ) * basis( y, g );
      }
		}

		return temp;
	}

	template <typename T>
	inline Vector<T> Spectral2D<T>::operator()( const Vector<T>& x )
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

		if ( BASIS == "EvenChebyshev" )
		{
			Chebyshev<T> basis;
			for ( std::size_t n = 0; n < N; ++n )
			{
				temp += basis( x, 2 * n ) * COEFFICIENTS[ n ];
			}
		}

		if ( BASIS == "OddChebyshev" )
		{
			Chebyshev<T> basis;
			for ( std::size_t n = 0; n < N; ++n )
			{
				temp += basis( x, 2 * n + 1 ) * COEFFICIENTS[ n ];
			}
		}

		if ( BASIS == "RationalSemi" )
		{
			RationalSemi<T> basis( PARAM );
			for ( std::size_t n = 0; n < N; ++n )
			{
				temp += basis( x, n ) * COEFFICIENTS[ n ];
			}
		}

		return temp;
	}

	template <typename T>
	inline T Spectral2D<T>::operator()( const T& x, const int& d )
	{
		T temp( 0.0 );
		std::size_t N( COEFFICIENTS.size() );

		if ( BASIS == "Chebyshev" )
		{
			Chebyshev<T> basis;
			for ( int n = 0; n < N; ++n )
			{
				temp += basis( x, n, d ) * COEFFICIENTS[ n ];
			}
		}

		if ( BASIS == "EvenChebyshev" )
		{
			Chebyshev<T> basis;
			for ( int n = 0; n < N; ++n )
			{
				temp += basis( x, 2 * n, d ) * COEFFICIENTS[ n ];
			}
		}

		if ( BASIS == "OddChebyshev" )
		{
			Chebyshev<T> basis;
			for ( int n = 0; n < N; ++n )
			{
				temp += basis( x, 2 * n + 1, d ) * COEFFICIENTS[ n ];
			}
		}

		if ( BASIS == "RationalSemi" )
		{
			RationalSemi<T> basis( PARAM );
			for ( int n = 0; n < N; ++n )
			{
				temp += basis( x, n, d ) * COEFFICIENTS[ n ];
			}
		}

		return temp;
	}

	template <typename T>
	inline Vector<T> Spectral2D<T>::operator()( const Vector<T>& x, const int& d )
	{
		Vector<T> temp( x.size() );
		std::size_t N( COEFFICIENTS.size() );

		if ( BASIS == "Chebyshev" )
		{
			Chebyshev<T> basis;
			for ( std::size_t n = 0; n < N; ++n )
			{
				temp += basis( x, n, d ) * COEFFICIENTS[ n ];
			}
		}

		if ( BASIS == "EvenChebyshev" )
		{
			Chebyshev<T> basis;
			for ( std::size_t n = 0; n < N; ++n )
			{
				temp += basis( x, 2 * n, d ) * COEFFICIENTS[ n ];
			}
		}

		if ( BASIS == "OddChebyshev" )
		{
			Chebyshev<T> basis;
			for ( std::size_t n = 0; n < N; ++n )
			{
				temp += basis( x, 2 * n + 1, d ) * COEFFICIENTS[ n ];
			}
		}

		if ( BASIS == "RationalSemi" )
		{
			RationalSemi<T> basis( PARAM );
			for ( std::size_t n = 0; n < N; ++n )
			{
				temp += basis( x, n, d ) * COEFFICIENTS[ n ];
			}
		}

		return temp;
	}

	template <typename T>
  inline const T& Spectral2D<T>::operator[]( const std::size_t& i ) const
  {
    if ( i < 0 || COEFFICIENTS.size() <= i )
    {
      throw Error( "Spectral operator [] range error" );
    }
    return COEFFICIENTS[ i ];
  }

  template <typename T>
  inline T& Spectral2D<T>::operator[] ( const std::size_t& i )
  {
    if ( i < 0 || COEFFICIENTS.size() <= i )
    {
      throw Error( "Spectral operator [] range error" );
    }
    return COEFFICIENTS[ i ];
  }

	/* ----- Methods ----- */

	template <typename T>
	std::string Spectral2D<T>::get_basis()
	{
		return BASIS;
	}

	template <typename T>
	void Spectral2D<T>::set_basis( const std::string& basis )
	{
		if ( basis != "Chebyshev" && basis != "EvenChebyshev"
		  && basis != "OddChebyshev" && basis != "RationalSemi" )
			//&& basis != "Fourier" )
		{
			std::string problem;
			problem += "Spectral set_basis error: The basis function " + basis + "\n";
			problem += "is not supported.";
			throw Error( problem );
		}
		BASIS = basis;
	}

	template <typename T>
	Vector<T> Spectral2D<T>::get_coefficients()
	{
		return COEFFICIENTS;
	}

	template <typename T>
	void Spectral2D<T>::set_coefficients( const Vector<T>& coeffs )
	{
		COEFFICIENTS = coeffs;
	}

	template <typename T>
	void Spectral2D<T>::update_coefficients( const Vector<T>& update )
	{
		if ( update.size() != COEFFICIENTS.size() )
		{
			throw Error( "Spectral: update_coefficients vector size incorrect." );
		}
		COEFFICIENTS += update;
	}

	template <typename T>
	void Spectral2D<T>::push_back( const T& new_elem )
	{
		COEFFICIENTS.push_back( new_elem );
	}

	template <typename T>
	void Spectral2D<T>::push_front( const T& new_elem )
	{
		COEFFICIENTS.push_front( new_elem );
	}

	template <typename T>
	double& Spectral2D<T>::parameter()
	{
		return PARAM;
	}

	template <typename T>
	int& Spectral2D<T>::x_points()
	{
		return I;
	}

	template <typename T>
	int& Spectral2D<T>::y_points()
	{
		return J;
	}

}  // End of namespace Luna

#endif
