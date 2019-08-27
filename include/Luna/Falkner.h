/// \file Falkner.h
/// A class for solving the modified Falkner-Skan equation.

#ifndef FALKNER_H
#define FALKNER_H

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
  // A templated class for spectral solutions
	template <class T>

	/// A class for solving the modified Falkner-Skan equation using a spectral
	/// rational Chebyshev approach. The equation \f[ u''' + (\eta + u)u'' -
	/// \beta u'\left(2 + u' \right) = 0 \f] is solved subject to the boundary
	/// conditions \f$ u(0) = 0\f$, \f$ u'(0) = -1 \f$ and \f$ u' \to 0 \f$ as \f$
	/// \eta \to \infty \f$. To recover the solution of the Falkner-Skan equation
	/// simply add \f$ \eta \f$ to the solution such that \f$ f(\eta) = u(\eta) +
	/// \eta \f$
	class Falkner
	{
		private:
			double BETA;				// Hartree parameter
			/// \todo TODO include parameter for injection BC -> K

    public:

			/// Constructor ( default )
			Falkner();

			/// Constructor with specified basis and coefficients
			/// \param beta The Hartree parameter
			Falkner( const double& beta );

			/// \todo TODO constructor for beta and K

			/// Destructor
			~Falkner(){}

			/* ----- Methods ----- */

			/// Solve the Falkner-Skan equation using a rational Chebyshev spectral
			/// method.
			/// \param n Number of spectral coefficients
			/// \param L Map parameter
			/// \param tol Tolerance for convergence of the solution
			/// \param initial_guess Initial guess function
			/// \param max_iter Maximum number of iterations
			/// \return The spectral solution
			Spectral<T> solve( const int& n, const double& L, const double& tol,
												 T initial_guess( const T& ), const int& max_iter	);



  }; // End of class Falkner

	template <typename T>
	Falkner<T>::Falkner()
	{
		BETA = 0.0;
	}

	template <typename T>
	Falkner<T>::Falkner( const double& beta )
	{
		BETA = beta;
	}

	/* ----- Methods ----- */

	template <typename T>
	Spectral<T> Falkner<T>::solve( const int& n, const double& L,
																 const double& tol, T initial_guess( const T& ),
																 const int& max_iter	)
	{
		Vector<double> y;                 			// Vector for collocation grid
		y.rational_semi_grid( n, L );
		y.push_back( 0.0 );
		RationalSemi<T> rationalsemi( L );
		Vector<T> c( n, 0.0 );                 	// Vector of coefficients

		// Approximate the guess as a rational Chebyshev polynomial
		c = rationalsemi.approximate( initial_guess, n );
		for ( std::size_t i = 0; i < 2; i++ )
		{
			c.push_back( 0.0 ); // Extra coefficients for BCs
		}
		// Make this into a spectral solution
		Spectral<double> u_g( c, "RationalSemi", L );

		double norm;
		int iter( 0 );

		do
		{
			// Setup the system of equations to find the correction
			Matrix<double> M( n + 2, n + 2, 0.0 );
			Vector<double> F( n + 2 ), a_c( n + 2 );

			// Don't need the BC at infinity as this is a "natural" BC

			// Internal nodes
			for ( std::size_t i = 0; i < n; i++ )
			{
				double yi = y[ i ];
				for ( std::size_t j = 0; j < n + 2; j++ )
				{
					// u_c'''
					M( i, j )  = rationalsemi( yi, j, 3 );
					// ( y + u_g ) * u_c''
					M( i, j ) += ( yi + u_g( yi ) ) * rationalsemi( yi, j, 2 );
					// - 2 * beta * ( 1 + u_g' ) * u_c'
					M( i, j ) -= 2 * BETA * ( 1 + u_g( yi, 1 ) ) * rationalsemi( yi, j, 1 );
					// u_g'' * u_c
					M( i, j ) += u_g( yi, 2 ) * rationalsemi( yi, j );
				}
				// - u_g''' - ( y + u_g ) * u_g'' + beta * u_g' ( 2 + u_g' )
				F[ i ]  = - u_g( yi, 3 ) - ( yi + u_g( yi ) ) * u_g( yi, 2 );
				F[ i ] += BETA * u_g( yi, 1 ) * ( 2 + u_g( yi, 1 ) );
			}

			// y = 0 boundary u_c = - u_g
			for ( std::size_t j = 0; j < n + 2; j++ )
			{
				double yi = 0.0;
				M( n, j ) = rationalsemi( yi, j );
			}
			F[ n ] = - u_g( 0.0 );

			// y = 0 boundary u_c' = -1 - u_g'
			for ( std::size_t j = 0; j < n + 2; j++ )
			{
				double yi = 0.0;
				M( n + 1, j ) = rationalsemi( yi, j, 1 );
			}
			F[ n + 1 ] = - 1.0 - u_g( 0.0, 1 );

			// Solve the system for the correction spectral coefficients
			a_c = M.solve( F );
			u_g.update_coefficients( a_c );
			norm = a_c.norm_2();
			++iter;

		}while( norm > tol && iter < max_iter );

		return u_g;
	}



}  // End of namespace Luna

#endif
