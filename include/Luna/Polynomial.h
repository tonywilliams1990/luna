/// \file Polynomial.h
/// A class for defining and evaluating polynomials and finding the roots of
/// polynomial equations.

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

			// Laguerre's method -> given complex coefficients of a polynomial a and a
			// complex value x, this routine improves x until it converges to a root
			// of the polynomial (to within the achievable round off limit).
			void laguer( Vector< std::complex<double> >& a, std::complex<double>& x,
			 						 std::size_t& iterations );

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

			/// Return the Vector of coefficients
			/// \return The Vector of coefficients
			Vector<T> coeffs();

			/// Return the degree of the Polynomial
			/// \return The degree of the Polynomial
			std::size_t degree();

			/// Set the degree of the Polynomial
			/// \param n The degree of the Polynomial
			void set_degree( const std::size_t& n );


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
			Vector<std::complex<double>> quadratic_solve( const T& a, const T& b,
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
			Vector<std::complex<double>> cubic_solve( const T& a, const T& b,
																								const T& c, const T& d );

			/// Divide the polynomial \f$u\f$ by another \f$v\f$ to find the quotient
			/// \f$q\f$ and the remainder \f$r\f$ such that \f$ u = vq + r  \f$
			/// \param v The divisor polynomial
			/// \param q The quotient polynomial
			/// \param r The remainder polynomial
			void polydiv( Polynomial<T>& v, Polynomial<T>& q, Polynomial<T>& r );

			/// Find all the roots of the polynomial using Laguerre's method
			/// \param polish True if you want to refine/polish the roots
			/// \return A Vector<std::complex<double>> containing the roots
			Vector<std::complex<double>> roots( const bool& polish );

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
	inline Vector<T> Polynomial<T>::coeffs(){
		return COEFFS;
	}

	template <typename T>
	inline std::size_t Polynomial<T>::degree()
	{
		return N;
	}

	template <typename T>
	inline void Polynomial<T>::set_degree( const std::size_t& n )
	{
		N = n;
		COEFFS.resize( n + 1 );
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

	template <typename T>
	inline void Polynomial<T>::polydiv( Polynomial<T>& v, Polynomial<T>& q,
																			Polynomial<T>& r )
	{
		int k, j, n=N, nv=v.N;
		if ( nv == 0 )
		{
			 throw Error( "Polynomial.polydiv() divide by zero polynomial");
		}
		while ( nv >= 0 && v[nv] == 0.) nv--;
		if ( nv < 0 )
		{
			 throw Error( "Polynomial.polydiv() divide by zero polynomial");
		}
		r.COEFFS = COEFFS;
		r.N = N;
		q.COEFFS.assign( COEFFS.size(), 0.0 );
		q.N = N;

		for ( k = n - nv; k >= 0; k-- )
		{
			q.COEFFS[ k ] = r.COEFFS[ nv + k ] / v.COEFFS[ nv ];
			for ( j = nv + k - 1; j >= k; j-- )
			{
				 r.COEFFS[ j ] -= q.COEFFS[ k ] * v.COEFFS[ j - k ];
			}
		}
		for ( j = nv; j <= n; j++ )
		{
			 r.COEFFS[ j ] = 0.0;
		}
		// Remove leading zeros from r and q polynomials
		T r_val( r.COEFFS[ n ] );
		while ( std::abs( r_val ) == 0.0 )
		{
			r.COEFFS.pop_back();
			r.N--;
			n--;
			r_val = r.COEFFS[ n ];
		}
		n = N;
		T q_val( q.COEFFS[ n ] );
		while ( std::abs( q_val ) == 0.0 )
		{
			q.COEFFS.pop_back();
			q.N--;
			n--;
			q_val = q.COEFFS[ n ];
		}
	}

	template <typename T>
	inline Vector<std::complex<double>> Polynomial<T>::roots( const bool& polish )
	{
		Vector<std::complex<double>> poly_roots( N );
		if ( N <= 0 )
		{
			std::string problem;
			problem = "Polynomial roots: degree must be at least 1";
			throw Error( problem );
		}
		if ( N == 1 )
		{
			poly_roots[ 0 ] = - COEFFS[ 0 ] / COEFFS[ 1 ];
			return poly_roots;
		}
		if ( N == 2 )
		{
			poly_roots = quadratic_solve( COEFFS[ 2 ], COEFFS[ 1 ], COEFFS[ 0 ] );
			return poly_roots;
		}
		if ( N == 3 )
		{
			poly_roots = cubic_solve( COEFFS[ 3 ], COEFFS[ 2 ], COEFFS[ 1 ],
																COEFFS[ 0 ] );
			return poly_roots;
		}

		const double EPS( 1.0e-14 );
		std::size_t i, its;
		std::complex<double> x, b, c;

		Vector< std::complex<double> > a( N + 1 );
		for ( std::size_t n = 0; n < N + 1; n++ )
		{
			a[ n ] = COEFFS[ n ];
		}

		Vector< std::complex<double> > ad( a );

		for ( int j = N - 1; j >= 0; j-- )
		{
			x = 0.0;
			Vector< std::complex<double> > ad_v( j + 2 );
			for ( int jj = 0; jj < j + 2; jj++ ) ad_v[ jj ] = ad[ jj ];
			laguer( ad_v, x, its );
			if ( std::abs( std::imag(x) ) <=
					 2.0 * EPS * std::abs( std::real( x ) ) )
			{
				x = std::complex<double>( std::real( x ), 0.0 );
			}
			poly_roots[ j ] = x;
			b = ad[ j + 1 ];
			for ( int jj = j; jj >= 0; jj-- ) {
				c = ad[ jj ];
				ad[ jj ] = b ;
				b = x * b + c;
			}
		}

		// Polish the roots
		if ( polish )
		{
			for ( std::size_t j = 0; j < N; j++ )
			{
				laguer( a, poly_roots[ j ], its );
			}
		}

		return poly_roots;
	}

	/* ----- Private functions ----- */
	template <typename T>
	inline void Polynomial<T>::laguer( Vector< std::complex<double> >& a,
																		 std::complex<double>& x,
																		 std::size_t& iterations )
	{
		const int MR=8, MT=10, MAXIT=MT*MR;
		const double EPS=std::numeric_limits<double>::epsilon();
		// Fraction used to break possible limit cycle
		static const double frac[MR+1]=
			{0.0,0.5,0.25,0.75,0.13,0.38,0.62,0.88,1.0};
		std::complex<double> dx,x1,b,d,f,g,h,sq,gp,gm,g2;
		int m=a.size()-1;
		for ( int iter = 1; iter <= MAXIT; iter++ )
		{
			iterations=iter;
			b = a[ m ];
			double err = std::abs( b );
			d = f = 0.0;
			double abx = std::abs( x );
			for (int j = m - 1; j >= 0 ; j-- )
			{
				f = x * f + d;
				d = x * d + b;
				b = x * b + a[ j ];
				err = std::abs( b ) + abx * err;
			}
			err *= EPS;
			if ( std::abs( b ) <= err ){ return; }
			g = d / b;
			g2 = g * g;
			h = g2 - 2.0 * f / b;
			sq = std::sqrt( double( m - 1 ) * ( double( m ) * h - g2 ) );
			gp = g + sq;
			gm = g - sq;
			double abp = std::abs( gp );
			double abm = std::abs( gm );
			if ( abp < abm ){ gp = gm; }
			dx = std::max( abp, abm ) > 0. ? double( m ) / gp :
																			 std::polar( 1 + abx, double( iter ));
			x1 = x - dx;
			if ( x == x1 ){ return; }
			if ( iter % MT != 0 ){ x = x1; }
			else { x -= frac[ iter / MT ] * dx; }
		}
		throw Error( "Polynomial: Too many iterations in laguer method" );
	}

}  // End of namespace Luna

#endif
