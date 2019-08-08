/// \file Spectral2D.h
/// A templated class for 2D spectral solutions

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
	template <class T>

	/// A class for defining 2D spectral solutions using specified basis functions
	/// and  coefficients. The spectral function \f$ u_{IJ}(x,y) \f$ is defined by
	/// \f[ u_{IJ}(x,y) = \sum_{N=0}^{IJ-1} a_N \Phi_N(x,y),  \f] where \f$ a_n
	/// \f$ are coefficients and \f$ \Phi_N(x,y) = \phi_{f(N)}(x) \phi_{g(N)}(y)
	/// \f$ are known 2D product basis functions.
	class Spectral2D
	{
		private:
			Vector<T> COEFFICIENTS; 	// Coefficients multiplying the basis functions
			std::string BASIS;				// Basis function specifier string
			double PARAM;							// Parameter for rational Chebyshev basis
			int I;										// Number of x collocation points
			int J;										// Number of y collocation points
			int CASE;

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
			/// \param i The number of collocation points in the x-direction
			/// \param j The number of collocation points in the y-direction
			/// \param basis The basis specifier string (i.e. Chebyshev)
			/// \param param An extra parameter for use with rational Chebyshev basis
			Spectral2D( const Vector<T>& coefficients, const int& i, const int& j,
									const std::string& basis, const double& param );
			/// \todo TODO modify this constructor to include number of collocation points

			/// Destructor
			~Spectral2D(){}

			/* ----- Operator overloading ----- */

			/// Evaluation operator at a point (x,y)
			/// \param x The x-coordinate where the spectral function is evaluated
			/// \param y The y-coordinate where the spectral function is evaluated
			/// \return The value of the spectral function at (x,y)
			T operator()( const T& x, const T& y );

			/// Evaluation operator on a 2D mesh
			/// \param mesh The Mesh2D object to be evaluated
			/// \param var The variable in the mesh to be filled
			void operator()( Mesh2D<T>& mesh, const int& var );

			/// Evaluation operator for the derivatives of the spectral function
			/// at a point (x,y)
			/// \param x The x-coordinate
			/// \param y The y-coordinate
			/// \param dx The x-derivative to return
			/// \param dy The y-derivative to return
			/// \return The dth derivative of the spectral function at the point x
			T operator()( const T& x, const T& y, const int& dx, const int& dy );

			/// Evaluation operator on a 2D mesh for the derivatives of the spectral
			/// function
			/// \param mesh The Mesh2D object to be evaluated
			/// \param var The variable in the mesh to be filled
			/// \param dx The x-derivative to return
			/// \param dy The y-derivative to return
			void operator()( Mesh2D<T>& mesh, const int& var, const int& dx,
				               const int& dy );

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

			/// Return a mesh containing the function on a given grid and its 1st and
			/// 2nd derivatives
			/// \param xnodes The mesh x nodes
			/// \param ynodes The mesh y nodes
			/// \return The spectral solution on a grid and its derivatives: u = 0,
			/// u_x = 1, u_y = 2, u_xx = 3, u_yy = 4, u_xy = 5
			Mesh2D<T> mesh_derivatives( const Vector<double>& xnodes,
																	const Vector<double>& ynodes );

  }; // End of class Spectral2D

	template <typename T>
	Spectral2D<T>::Spectral2D()
	{
		BASIS = "Chebyshev";
		CASE = 1;
	}

	template <typename T>
	Spectral2D<T>::Spectral2D( const Vector<T>& coefficients, const int& i,
														 const int& j, const std::string& basis )
	{
		if ( basis != "Chebyshev" && basis != "EvenChebyshev"
		  && basis != "OddChebyshev" && basis != "RationalSemi"
		  && basis != "EvenRationalSemi" && basis != "OddRationalSemi" )
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
		if ( BASIS == "Chebyshev" ){ CASE = 1; }
		if ( BASIS == "EvenChebyshev" ){ CASE = 2; }
		if ( BASIS == "OddChebyshev" ){ CASE = 3; }
		if ( BASIS == "RationalSemi" ){ CASE = 4; }
		if ( BASIS == "EvenRationalSemi" ){ CASE = 5; }
		if ( BASIS == "OddRationalSemi" ){ CASE = 6; }
	}

	template <typename T>
	Spectral2D<T>::Spectral2D( const Vector<T>& coefficients, const int& i,
							 const int& j, const std::string& basis, const double& param )
	{
		if ( basis != "Chebyshev" && basis != "EvenChebyshev"
		  && basis != "OddChebyshev" && basis != "RationalSemi"
			&& basis != "EvenRationalSemi" && basis != "OddRationalSemi" )
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
		I = i;
		J = j;
		if ( BASIS == "Chebyshev" ){ CASE = 1; }
		if ( BASIS == "EvenChebyshev" ){ CASE = 2; }
		if ( BASIS == "OddChebyshev" ){ CASE = 3; }
		if ( BASIS == "RationalSemi" ){ CASE = 4; }
		if ( BASIS == "EvenRationalSemi" ){ CASE = 5; }
		if ( BASIS == "OddRationalSemi" ){ CASE = 6; }
	}

	/* ----- Operator overloading ----- */

	template <typename T>
	inline T Spectral2D<T>::operator()( const T& x, const T& y )
	{
		T temp( 0 );
		int size( I * J );
		int f, g;

		// Chebyshev
		if ( CASE == 1 )
		{
			Chebyshev<T> basis;
			for ( std::size_t n = 0; n < size; n++ )
      {
        f = n / J;
        g = n % J;
        temp += COEFFICIENTS[ n ] * basis( x, f ) * basis( y, g );
      }
		}

		// EvenChebyshev
		if ( CASE == 2 )
		{
			Chebyshev<T> basis;
			for ( std::size_t n = 0; n < size; n++ )
      {
        f = n / J;
        g = n % J;
        temp += COEFFICIENTS[ n ] * basis( x, 2 * f ) * basis( y, 2 * g );
      }
		}

		// OddChebyshev
		if ( CASE == 3 )
		{
			Chebyshev<T> basis;
			for ( std::size_t n = 0; n < size; n++ )
      {
        f = n / J;
        g = n % J;
        temp += COEFFICIENTS[ n ] * basis( x, 2 * f + 1 )
																	* basis( y, 2 * g + 1 );
      }
		}

		// RationalSemi
		if ( CASE == 4 )
		{
			RationalSemi<T> basis( PARAM );
			for ( std::size_t n = 0; n < size; n++ )
      {
        f = n / J;
        g = n % J;
        temp += COEFFICIENTS[ n ] * basis( x, f )	* basis( y, g );
      }
		}

		// EvenRationalSemi
		if ( CASE == 5 )
		{
			RationalSemi<T> basis( PARAM );
			for ( std::size_t n = 0; n < size; n++ )
      {
        f = n / J;
        g = n % J;
        temp += COEFFICIENTS[ n ] * basis( x, 2 * f ) * basis( y, 2 * g );
      }
		}

		// OddRationalSemi
		if ( CASE == 6 )
		{
			RationalSemi<T> basis( PARAM );
			for ( std::size_t n = 0; n < size; n++ )
      {
        f = n / J;
        g = n % J;
        temp += COEFFICIENTS[ n ] * basis( x, 2 * f + 1 )
																	* basis( y, 2 * g + 1 );
      }
		}

		return temp;
	}

	template <typename T>
	inline void Spectral2D<T>::operator()( Mesh2D<T>& mesh, const int& var )
	{
		if ( var >= mesh.get_nvars() )
		{
			throw Error( "Spectral2D error: Incorrect variable number." );
		}
		int nx( mesh.xnodes().size() );
		int ny( mesh.ynodes().size() );
		int f, g;
		int size( I * J );

		if ( CASE == 1 )
		{
			Chebyshev<T> basis;
			for ( std::size_t i = 0; i < nx; i++ )
			{
				double x( mesh.xnodes()[ i ] );
				for ( std::size_t j = 0; j < ny; j++ )
				{
					double y( mesh.ynodes()[ j ] );
					for ( std::size_t n = 0; n < size; n++ )
					{
						f = n / J;
						g = n % J;
						mesh( i, j, var ) += COEFFICIENTS[ n ] * basis( x, f )
																									 * basis( y, g );
					}
				}
			}
		}

		if ( CASE == 2 )
		{
			Chebyshev<T> basis;
			for ( std::size_t i = 0; i < nx; i++ )
			{
				double x( mesh.xnodes()[ i ] );
				for ( std::size_t j = 0; j < ny; j++ )
				{
					double y( mesh.ynodes()[ j ] );
					for ( std::size_t n = 0; n < size; n++ )
					{
						f = n / J;
						g = n % J;
						mesh( i, j, var ) += COEFFICIENTS[ n ] * basis( x, 2 * f )
																									 * basis( y, 2 * g );
					}
				}
			}
		}

		if ( CASE == 3 )
		{
			Chebyshev<T> basis;
			for ( std::size_t i = 0; i < nx; i++ )
			{
				double x( mesh.xnodes()[ i ] );
				for ( std::size_t j = 0; j < ny; j++ )
				{
					double y( mesh.ynodes()[ j ] );
					for ( std::size_t n = 0; n < size; n++ )
					{
						f = n / J;
						g = n % J;
						mesh( i, j, var ) += COEFFICIENTS[ n ] * basis( x, 2 * f + 1 )
																									 * basis( y, 2 * g + 1 );
					}
				}
			}
		}

		if ( CASE == 4 )
		{
			RationalSemi<T> basis( PARAM );
			for ( std::size_t i = 0; i < nx; i++ )
			{
				double x( mesh.xnodes()[ i ] );
				for ( std::size_t j = 0; j < ny; j++ )
				{
					double y( mesh.ynodes()[ j ] );
					for ( std::size_t n = 0; n < size; n++ )
					{
						f = n / J;
						g = n % J;
						mesh( i, j, var ) += COEFFICIENTS[ n ] * basis( x, f )
																									 * basis( y, g );
					}
				}
			}
		}

		if ( CASE == 5 )
		{
			RationalSemi<T> basis( PARAM );
			for ( std::size_t i = 0; i < nx; i++ )
			{
				double x( mesh.xnodes()[ i ] );
				for ( std::size_t j = 0; j < ny; j++ )
				{
					double y( mesh.ynodes()[ j ] );
					for ( std::size_t n = 0; n < size; n++ )
					{
						f = n / J;
						g = n % J;
						mesh( i, j, var ) += COEFFICIENTS[ n ] * basis( x, 2 * f )
																									 * basis( y, 2 * g );
					}
				}
			}
		}

		if ( CASE == 6 )
		{
			RationalSemi<T> basis( PARAM );
			for ( std::size_t i = 0; i < nx; i++ )
			{
				double x( mesh.xnodes()[ i ] );
				for ( std::size_t j = 0; j < ny; j++ )
				{
					double y( mesh.ynodes()[ j ] );
					for ( std::size_t n = 0; n < size; n++ )
					{
						f = n / J;
						g = n % J;
						mesh( i, j, var ) += COEFFICIENTS[ n ] * basis( x, 2 * f + 1 )
																									 * basis( y, 2 * g + 1 );
					}
				}
			}
		}

	}

	template <typename T>
	inline T Spectral2D<T>::operator()( const T& x, const T& y, const int& dx,
																			const int& dy )
	{
		T temp( 0 );
		int size( I * J );
		int f, g;

		if ( CASE == 1 )
		{
			Chebyshev<T> basis;
			for ( std::size_t n = 0; n < size; n++ )
			{
				f = n / J;
				g = n % J;
				temp += COEFFICIENTS[ n ] * basis( x, f, dx ) * basis( y, g, dy );
			}
		}

		if ( CASE == 2 )
		{
			Chebyshev<T> basis;
			for ( std::size_t n = 0; n < size; n++ )
			{
				f = n / J;
				g = n % J;
				temp += COEFFICIENTS[ n ] * basis( x, 2 * f, dx )
																	* basis( y, 2 * g, dy );
			}
		}

		if ( CASE == 3 )
		{
			Chebyshev<T> basis;
			for ( std::size_t n = 0; n < size; n++ )
			{
				f = n / J;
				g = n % J;
				temp += COEFFICIENTS[ n ] * basis( x, 2 * f + 1, dx )
																	* basis( y, 2 * g + 1, dy );
			}
		}

		if ( CASE == 4 )
		{
			RationalSemi<T> basis( PARAM );
			for ( std::size_t n = 0; n < size; n++ )
			{
				f = n / J;
				g = n % J;
				temp += COEFFICIENTS[ n ] * basis( x, f, dx ) * basis( y, g, dy );
			}
		}

		if ( CASE == 5 )
		{
			RationalSemi<T> basis( PARAM );
			for ( std::size_t n = 0; n < size; n++ )
			{
				f = n / J;
				g = n % J;
				temp += COEFFICIENTS[ n ] * basis( x, 2 * f, dx )
																	* basis( y, 2 * g, dy );
			}
		}

		if ( CASE == 6 )
		{
			RationalSemi<T> basis( PARAM );
			for ( std::size_t n = 0; n < size; n++ )
			{
				f = n / J;
				g = n % J;
				temp += COEFFICIENTS[ n ] * basis( x, 2 * f + 1, dx )
																	* basis( y, 2 * g + 1, dy );
			}
		}

		return temp;
	}

	template <typename T>
	inline void Spectral2D<T>::operator()( Mesh2D<T>& mesh, const int& var,
																				 const int& dx, const int& dy )
	{
		if ( var >= mesh.get_nvars() )
		{
			throw Error( "Spectral2D error: Incorrect variable number." );
		}
		int nx( mesh.xnodes().size() );
		int ny( mesh.ynodes().size() );
		int f, g;
		int size( I * J );

		if ( CASE == 1 )
		{
			Chebyshev<T> basis;
			for ( std::size_t i = 0; i < nx; i++ )
			{
				double x( mesh.xnodes()[ i ] );
				for ( std::size_t j = 0; j < ny; j++ )
				{
					double y( mesh.ynodes()[ j ] );
					for ( std::size_t n = 0; n < size; n++ )
					{
						f = n / J;
						g = n % J;
						mesh( i, j, var ) += COEFFICIENTS[ n ] * basis( x, f, dx )
																									 * basis( y, g, dy );
					}
				}
			}
		}

		if ( CASE == 2 )
		{
			Chebyshev<T> basis;
			for ( std::size_t i = 0; i < nx; i++ )
			{
				double x( mesh.xnodes()[ i ] );
				for ( std::size_t j = 0; j < ny; j++ )
				{
					double y( mesh.ynodes()[ j ] );
					for ( std::size_t n = 0; n < size; n++ )
					{
						f = n / J;
						g = n % J;
						mesh( i, j, var ) += COEFFICIENTS[ n ] * basis( x, 2 * f, dx )
																									 * basis( y, 2 * g, dy );
					}
				}
			}
		}

		if ( CASE == 3 )
		{
			Chebyshev<T> basis;
			for ( std::size_t i = 0; i < nx; i++ )
			{
				double x( mesh.xnodes()[ i ] );
				for ( std::size_t j = 0; j < ny; j++ )
				{
					double y( mesh.ynodes()[ j ] );
					for ( std::size_t n = 0; n < size; n++ )
					{
						f = n / J;
						g = n % J;
						mesh( i, j, var ) += COEFFICIENTS[ n ] * basis( x, 2 * f + 1, dx )
																									 * basis( y, 2 * g + 1, dy );
					}
				}
			}
		}

		if ( CASE == 4 )
		{
			RationalSemi<T> basis( PARAM );
			for ( std::size_t i = 0; i < nx; i++ )
			{
				double x( mesh.xnodes()[ i ] );
				for ( std::size_t j = 0; j < ny; j++ )
				{
					double y( mesh.ynodes()[ j ] );
					for ( std::size_t n = 0; n < size; n++ )
					{
						f = n / J;
						g = n % J;
						mesh( i, j, var ) += COEFFICIENTS[ n ] * basis( x, f, dx )
																									 * basis( y, g, dy );
					}
				}
			}
		}

		if ( CASE == 5 )
		{
			RationalSemi<T> basis( PARAM );
			for ( std::size_t i = 0; i < nx; i++ )
			{
				double x( mesh.xnodes()[ i ] );
				for ( std::size_t j = 0; j < ny; j++ )
				{
					double y( mesh.ynodes()[ j ] );
					for ( std::size_t n = 0; n < size; n++ )
					{
						f = n / J;
						g = n % J;
						mesh( i, j, var ) += COEFFICIENTS[ n ] * basis( x, 2 * f, dx )
																									 * basis( y, 2 * g, dy );
					}
				}
			}
		}

		if ( CASE == 6 )
		{
			RationalSemi<T> basis( PARAM );
			for ( std::size_t i = 0; i < nx; i++ )
			{
				double x( mesh.xnodes()[ i ] );
				for ( std::size_t j = 0; j < ny; j++ )
				{
					double y( mesh.ynodes()[ j ] );
					for ( std::size_t n = 0; n < size; n++ )
					{
						f = n / J;
						g = n % J;
						mesh( i, j, var ) += COEFFICIENTS[ n ] * basis( x, 2 * f + 1, dx )
																									 * basis( y, 2 * g + 1, dy );
					}
				}
			}
		}

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
		  && basis != "OddChebyshev" && basis != "RationalSemi"
			&& basis != "EvenRationalSemi" && basis != "OddRationalSemi" )
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

	template <typename T>
	Mesh2D<T> Spectral2D<T>::mesh_derivatives( const Vector<double>& xnodes,
																						 const Vector<double>& ynodes )
	{
		Mesh2D<double> mesh( xnodes, ynodes, 6 );

		int nx( xnodes.size() );
		int ny( ynodes.size() );
		int f, g;
		int size( I * J );

		if ( CASE == 1 )
		{
			Chebyshev<T> basis;
			Vector<double> phix, phiy;
			for ( std::size_t i = 0; i < nx; i++ )
			{
				double x( xnodes[ i ] );
				for ( std::size_t j = 0; j < ny; j++ )
				{
					double y( ynodes[ j ] );
					for ( std::size_t n = 0; n < size; n++ )
					{
						f = n / J;
						g = n % J;
						phix = basis.eval_2( x, f );  // phix, phix' and phix''
		        phiy = basis.eval_2( y, g );  // phiy, phiy' and phiy''
						mesh( i, j, 0 ) += COEFFICIENTS[ n ] * phix[ 0 ]
																								 * phiy[ 0 ];	// u
						mesh( i, j, 1 ) += COEFFICIENTS[ n ] * phix[ 1 ]
																								 * phiy[ 0 ];	// u_x
						mesh( i, j, 2 ) += COEFFICIENTS[ n ] * phix[ 0 ]
																								 * phiy[ 1 ]; // u_y
						mesh( i, j, 3 ) += COEFFICIENTS[ n ] * phix[ 2 ]
																								 * phiy[ 0 ]; // u_xx
						mesh( i, j, 4 ) += COEFFICIENTS[ n ] * phix[ 0 ]
																								 * phiy[ 2 ]; // u_yy
						mesh( i, j, 5 ) += COEFFICIENTS[ n ] * phix[ 1 ]
																								 * phiy[ 1 ]; // u_xy
						/*mesh( i, j, 0 ) += COEFFICIENTS[ n ] * basis( x, f )
																								 * basis( y, g );	   // u
						mesh( i, j, 1 ) += COEFFICIENTS[ n ] * basis( x, f, 1 )
																								 * basis( y, g );	   // u_x
						mesh( i, j, 2 ) += COEFFICIENTS[ n ] * basis( x, f )
																								 * basis( y, g, 1 ); // u_y
						mesh( i, j, 3 ) += COEFFICIENTS[ n ] * basis( x, f, 2 )
																								 * basis( y, g );		 // u_xx
						mesh( i, j, 4 ) += COEFFICIENTS[ n ] * basis( x, f )
																								 * basis( y, g, 2 ); // u_yy
						mesh( i, j, 5 ) += COEFFICIENTS[ n ] * basis( x, f, 1 )
																								 * basis( y, g, 1 ); // u_xy*/
					}
				}
			}
		}

		if ( CASE == 2 )
		{
			Chebyshev<T> basis;
			Vector<double> phix, phiy;
			for ( std::size_t i = 0; i < nx; i++ )
			{
				double x( xnodes[ i ] );
				for ( std::size_t j = 0; j < ny; j++ )
				{
					double y( ynodes[ j ] );
					for ( std::size_t n = 0; n < size; n++ )
					{
						f = n / J;
						g = n % J;
						phix = basis.eval_2( x, 2 * f );  // phix, phix' and phix''
		        phiy = basis.eval_2( y, 2 * g );  // phiy, phiy' and phiy''
						mesh( i, j, 0 ) += COEFFICIENTS[ n ] * phix[ 0 ]
																								 * phiy[ 0 ];	// u
						mesh( i, j, 1 ) += COEFFICIENTS[ n ] * phix[ 1 ]
																								 * phiy[ 0 ];	// u_x
						mesh( i, j, 2 ) += COEFFICIENTS[ n ] * phix[ 0 ]
																								 * phiy[ 1 ]; // u_y
						mesh( i, j, 3 ) += COEFFICIENTS[ n ] * phix[ 2 ]
																								 * phiy[ 0 ]; // u_xx
						mesh( i, j, 4 ) += COEFFICIENTS[ n ] * phix[ 0 ]
																								 * phiy[ 2 ]; // u_yy
						mesh( i, j, 5 ) += COEFFICIENTS[ n ] * phix[ 1 ]
																								 * phiy[ 1 ]; // u_xy
					}
				}
			}
		}

		if ( CASE == 3 )
		{
			Chebyshev<T> basis;
			Vector<double> phix, phiy;
			for ( std::size_t i = 0; i < nx; i++ )
			{
				double x( xnodes[ i ] );
				for ( std::size_t j = 0; j < ny; j++ )
				{
					double y( ynodes[ j ] );
					for ( std::size_t n = 0; n < size; n++ )
					{
						f = n / J;
						g = n % J;
						phix = basis.eval_2( x, 2 * f + 1 );  // phix, phix' and phix''
		        phiy = basis.eval_2( y, 2 * g + 1 );  // phiy, phiy' and phiy''
						mesh( i, j, 0 ) += COEFFICIENTS[ n ] * phix[ 0 ]
																								 * phiy[ 0 ];	// u
						mesh( i, j, 1 ) += COEFFICIENTS[ n ] * phix[ 1 ]
																								 * phiy[ 0 ];	// u_x
						mesh( i, j, 2 ) += COEFFICIENTS[ n ] * phix[ 0 ]
																								 * phiy[ 1 ]; // u_y
						mesh( i, j, 3 ) += COEFFICIENTS[ n ] * phix[ 2 ]
																								 * phiy[ 0 ]; // u_xx
						mesh( i, j, 4 ) += COEFFICIENTS[ n ] * phix[ 0 ]
																								 * phiy[ 2 ]; // u_yy
						mesh( i, j, 5 ) += COEFFICIENTS[ n ] * phix[ 1 ]
																								 * phiy[ 1 ]; // u_xy
					}
				}
			}
		}

		if ( CASE == 4 )
		{
			RationalSemi<T> basis( PARAM );
			Vector<double> phix, phiy;
			for ( std::size_t i = 0; i < nx; i++ )
			{
				double x( xnodes[ i ] );
				for ( std::size_t j = 0; j < ny; j++ )
				{
					double y( ynodes[ j ] );
					for ( std::size_t n = 0; n < size; n++ )
					{
						f = n / J;
						g = n % J;
						phix = basis.eval_2( x, f );  // phix, phix' and phix''
		        phiy = basis.eval_2( y, g );  // phiy, phiy' and phiy''
						mesh( i, j, 0 ) += COEFFICIENTS[ n ] * phix[ 0 ]
																								 * phiy[ 0 ];	// u
						mesh( i, j, 1 ) += COEFFICIENTS[ n ] * phix[ 1 ]
																								 * phiy[ 0 ];	// u_x
						mesh( i, j, 2 ) += COEFFICIENTS[ n ] * phix[ 0 ]
																								 * phiy[ 1 ]; // u_y
						mesh( i, j, 3 ) += COEFFICIENTS[ n ] * phix[ 2 ]
																								 * phiy[ 0 ]; // u_xx
						mesh( i, j, 4 ) += COEFFICIENTS[ n ] * phix[ 0 ]
																								 * phiy[ 2 ]; // u_yy
						mesh( i, j, 5 ) += COEFFICIENTS[ n ] * phix[ 1 ]
																								 * phiy[ 1 ]; // u_xy
					}
				}
			}
		}

		if ( CASE == 5 )
		{
			RationalSemi<T> basis( PARAM );
			Vector<double> phix, phiy;
			for ( std::size_t i = 0; i < nx; i++ )
			{
				double x( xnodes[ i ] );
				for ( std::size_t j = 0; j < ny; j++ )
				{
					double y( ynodes[ j ] );
					for ( std::size_t n = 0; n < size; n++ )
					{
						f = n / J;
						g = n % J;
						phix = basis.eval_2( x, 2 * f );  // phix, phix' and phix''
		        phiy = basis.eval_2( y, 2 * g );  // phiy, phiy' and phiy''
						mesh( i, j, 0 ) += COEFFICIENTS[ n ] * phix[ 0 ]
																								 * phiy[ 0 ];	// u
						mesh( i, j, 1 ) += COEFFICIENTS[ n ] * phix[ 1 ]
																								 * phiy[ 0 ];	// u_x
						mesh( i, j, 2 ) += COEFFICIENTS[ n ] * phix[ 0 ]
																								 * phiy[ 1 ]; // u_y
						mesh( i, j, 3 ) += COEFFICIENTS[ n ] * phix[ 2 ]
																								 * phiy[ 0 ]; // u_xx
						mesh( i, j, 4 ) += COEFFICIENTS[ n ] * phix[ 0 ]
																								 * phiy[ 2 ]; // u_yy
						mesh( i, j, 5 ) += COEFFICIENTS[ n ] * phix[ 1 ]
																								 * phiy[ 1 ]; // u_xy
					}
				}
			}
		}

		if ( CASE == 6 )
		{
			RationalSemi<T> basis( PARAM );
			Vector<double> phix, phiy;
			for ( std::size_t i = 0; i < nx; i++ )
			{
				double x( xnodes[ i ] );
				for ( std::size_t j = 0; j < ny; j++ )
				{
					double y( ynodes[ j ] );
					for ( std::size_t n = 0; n < size; n++ )
					{
						f = n / J;
						g = n % J;
						phix = basis.eval_2( x, 2 * f + 1 );  // phix, phix' and phix''
		        phiy = basis.eval_2( y, 2 * g + 1 );  // phiy, phiy' and phiy''
						mesh( i, j, 0 ) += COEFFICIENTS[ n ] * phix[ 0 ]
																								 * phiy[ 0 ];	// u
						mesh( i, j, 1 ) += COEFFICIENTS[ n ] * phix[ 1 ]
																								 * phiy[ 0 ];	// u_x
						mesh( i, j, 2 ) += COEFFICIENTS[ n ] * phix[ 0 ]
																								 * phiy[ 1 ]; // u_y
						mesh( i, j, 3 ) += COEFFICIENTS[ n ] * phix[ 2 ]
																								 * phiy[ 0 ]; // u_xx
						mesh( i, j, 4 ) += COEFFICIENTS[ n ] * phix[ 0 ]
																								 * phiy[ 2 ]; // u_yy
						mesh( i, j, 5 ) += COEFFICIENTS[ n ] * phix[ 1 ]
																								 * phiy[ 1 ]; // u_xy
					}
				}
			}
		}

		return mesh;
	}

}  // End of namespace Luna

#endif
