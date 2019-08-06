/// \file SVD.h
/// Defines a class for performing singular value decompositions of matrices.

#ifndef SVD_H
#define SVD_H

#include <ctime>
#include <iostream>
#include <limits>
#include <numeric>
#include "Matrix.h"
#include "Vector.h"
#include "Error.h"

namespace Luna
{
	template <class T>

	/// A class for performing singular value decompositions of matrices.
	class SVD
	{

	private:
		int m, n;
		Matrix<T> u, v;
		Vector<T> w;
		double eps, tsh;

		double inv_condition() {
			return ( w[ 0 ] <= 0. || w[ n - 1 ] <= 0. ) ? 0. : w[ n - 1 ] / w[ 0 ];
		}

		/// \todo TODO complex SVD -> algo358
		void decompose();
		void reorder();
		double pythag( const double a, const double b );
		double SIGN( const double &a, const double &b );

	public:

		/* ----- Constructors and Destructor ----- */

		/// Constructor
		SVD( Matrix<T> &a ) : m( a.rows() ), n( a.cols() ), u( a ), v( n, n, 0.0 ),
													w( n )
		{
			eps = std::numeric_limits<double>::epsilon();
			decompose();
			//reorder();
			tsh = 0.5 * std::sqrt( m + n + 1. ) * w[ 0 ] * eps;
		}

		/// Destructor
		~SVD()
		{}

		/* ----- Methods ----- */

		/// Solve the system of equations Ax=b where x and b are Vectors
    /// \param b The right-hand side Vector of the system of equations
    /// \param x The solution Vector
		/// \param thresh Threshold for nonzero singular values
		void solve( Vector<T>& b, Vector<T>& x, double thresh = - 1. );

		/// Solve the system of equations AX=B where X and B are Matrices
    /// \param B The right-hand side Matrix of the system of equations
    /// \param x The solution Matrix
		/// \param thresh Threshold for nonzero singular values
		void solve( Matrix<T>& B, Matrix<T>& X, double thresh = - 1. );

		/// Return the rank of the Matrix
		/// \param thresh Threshold for nonzero singular values ( negative value
		/// leads to default being used based on estimated roundoff )
		int rank( double thresh = - 1. );

		/// Return the nullity of the Matrix
		/// \param thresh Threshold for nonzero singular values ( negative value
		/// leads to default being used based on estimated roundoff )
		int nullity( double thresh = - 1. );

		/// Return an orthonormal basis for the range of the Matrix as the columns
		/// of the returned Matrix
		/// \param thresh Threshold for nonzero singular values ( negative value
		/// leads to default being used based on estimated roundoff )
		Matrix<T> range( double thresh = - 1. );

		/// Return an orthonormal basis for the nullspace of the Matrix as the
		///columns of the returned Matrix
		/// \param thresh Threshold for nonzero singular values ( negative value
		/// leads to default being used based on estimated roundoff )
		Matrix<T> nullspace( double thresh = - 1. );

		/// Return the Vector of singular values
		/// \return The singular values of the Matrix
		Vector<T> singular_values();

		/// Return the left-singular Vectors of the Matrix
		/// \return The left-singular Vectors
		Matrix<T> left_vectors();

		/// Return the right-singular Vectors of the Matrix
		/// \return The right-singular Vectors
		Matrix<T> right_vectors();

	}; // End of SVD class



	/* ----- Methods ----- */

	template <typename T>
	inline void SVD<T>::solve( Vector<T>& b, Vector<T>& x, double thresh )
	{
		int i, j, jj;
		T s;
		if ( b.size() != m || x.size() != n ){
			throw Error( "SVD::solve bad sizes" );
		}
		Vector<T> tmp( n );
		tsh = ( thresh >= 0. ? thresh : 0.5 * sqrt( m + n + 1. ) * w[ 0 ] * eps );
		for ( j = 0; j < n; j++ ) {
			s = 0.0;
			if ( w[ j ] > tsh ) {
				for ( i = 0; i < m; i++ ) {
					s += u[ i ][ j ] * b[ i ];
				}
				s /= w[ j ];
			}
			tmp[ j ] = s;
		}
		for ( j = 0; j < n; j++ ) {
			s = 0.0;
			for ( jj = 0; jj < n; jj++ ) {
				s += v[ j ][ jj ] * tmp[ jj ];
			}
			x[ j ] = s;
		}
	}

	template <typename T>
	inline void SVD<T>::solve( Matrix<T>& b, Matrix<T>& x, double thresh )
	{
		int i, j, p = b.ncols();
		if ( b.nrows() != m || x.nrows() != n || x.ncols() != p )
		{
			throw Error( "SVD::solve bad sizes" );
		}
		Vector<T> xx( n ), bcol( m );
		for ( j = 0; j < p; j++ ) {
			for ( i = 0; i < m; i++ ) {
				bcol[ i ] = b[ i ][ j ];
			}
			solve( bcol, xx, thresh );
			for ( i = 0; i < n; i++ ) {
				x[ i ][ j ] = xx[ i ];
			}
		}
	}

	template <typename T>
	inline int SVD<T>::rank( double thresh )
	{
		int j, nr = 0;
		tsh = ( thresh >= 0. ? thresh : 0.5 * sqrt( m + n + 1. ) * w[ 0 ] * eps );
		for ( j = 0; j < n; j++ ){
			if ( w[ j ] > tsh ) nr++;
		}
		return nr;
	}

	template <typename T>
	inline int SVD<T>::nullity( double thresh )
	{
		int j, nn = 0;
		tsh = ( thresh >= 0. ? thresh : 0.5 * sqrt( m + n + 1. ) * w[ 0 ] * eps );
		for ( j = 0; j < n; j++ ) {
			if ( w[ j ] <= tsh ) nn++;
		}
		return nn;
	}

	template <typename T>
	inline Matrix<T> SVD<T>::range( double thresh )
	{
		int i, j, nr = 0;
		Matrix<T> rnge( m, rank( thresh ) );
		for ( j = 0; j < n; j++ ) {
			if ( w[ j ] > tsh ) {
				for ( i = 0; i < m ; i++ ) {
					rnge[ i ][ nr ] = u[ i ][ j ];
				}
				nr++;
			}
		}
		return rnge;
	}

	template <typename T>
	inline Matrix<T> SVD<T>::nullspace( double thresh )
	{
		int j, jj, nn = 0;
		Matrix<T> nullsp( n, nullity( thresh ), 0.0 );
		for ( j = 0; j < n; j++ ) {
			if ( w[ j ] <= tsh ) {
				for ( jj = 0; jj < n; jj++ ) {
					nullsp[ jj ][ nn ] = v[ jj ][ j ];
				}
				nn++;
			}
		}
		return nullsp;
	}

	template <typename T>
	inline Vector<T> SVD<T>::singular_values()
	{
		return w;
	}

	template <typename T>
	inline Matrix<T> SVD<T>::left_vectors()
	{
		return u;
	}

	template <typename T>
	inline Matrix<T> SVD<T>::right_vectors()
	{
		return v;
	}

	/* ----- Private methods ----- */

	template <>
	inline void SVD<double>::decompose() {
		bool flag;
		int i, its, j, jj, k, l, nm;
		double anorm, c, f, g, h, s, scale, x, y, z;
		Vector<double> rv1( n );
		g = scale = anorm = 0.0;
		for ( i = 0; i < n; i++ ) {
			l = i + 2;
			rv1[ i ] = scale * g;
			g = s = scale = 0.0;
			if ( i < m ) {
				for ( k = i; k < m; k++ ) scale += std::abs( u[ k ][ i ] );
				if ( scale != 0.0 ) {
					for ( k = i; k < m; k++ ) {
						u[ k ][ i ] /= scale;
						s += u[ k ][ i ] * u[ k ][ i ];
					}
					f = u[ i ][ i ];
					g = - SIGN( std::sqrt( s ), f );
					h = f * g - s;
					u[ i ][ i ] = f - g;
					for ( j = l - 1; j < n; j++ ) {
						for ( s = 0.0, k = i; k < m; k++ ){
							s += u[ k ][ i ] * u[ k ][ j ];
						}
						f = s / h;
						for ( k = i; k < m; k++ ){
							u[ k ][ j ] += f * u[ k ][ i ];
						}
					}
					for ( k = i; k < m; k++ ){
						u[ k ][ i ] *= scale;
					}
				}
			}
			w[ i ] = scale * g;
			g = s = scale = 0.0;
			if ( i + 1 <= m && i + 1 != n ) {
				for ( k = l - 1; k < n; k++ ){
					scale += std::abs( u[ i ][ k ] );
				}
				if ( scale != 0.0 ) {
					for ( k = l - 1; k < n; k++ ) {
						u[ i ][ k ] /= scale;
						s += u[ i ][ k ] * u[ i ][ k ];
					}
					f = u[ i ][ l - 1 ];
					g = -SIGN( std::sqrt( s ), f );
					h = f * g - s;
					u[ i ][ l - 1 ] = f - g;
					for ( k = l - 1; k < n; k++ ){
						rv1[ k ] = u[ i ][ k ] / h;
					}
					for ( j = l - 1; j < m; j++ ) {
						for ( s = 0.0, k = l - 1; k < n; k++ ){
							s += u[ j ][ k ] * u[ i ][ k ];
						}
						for ( k = l - 1; k < n; k++ ){
							u[ j ][ k ] += s * rv1[ k ];
						}
					}
					for ( k = l - 1; k < n; k++ ){
						u[ i ][ k ] *= scale;
					}
				}
			}
			anorm = std::max( anorm, ( std::abs( w[ i ] )+ std::abs( rv1[ i ] ) ) );
		}
		for ( i = n - 1; i >= 0; i-- ) {
			if ( i < n - 1 ) {
				if ( g != 0.0 ) {
					for ( j = l; j < n; j++ ) {
						v[ j ][ i ] = ( u[ i ][ j ] / u[ i ][ l ] ) / g;
					}
					for ( j = l; j < n; j++ ) {
						for ( s = 0.0, k = l; k < n; k++ ){
							s += u[ i ][ k ] * v[ k ][ j ];
						}
						for ( k = l; k < n; k++ ){
							v[ k ][ j ] += s * v[ k ][ i ];
						}
					}
				}
				for ( j = l; j < n; j++ ){
					v[ i ][ j ] = v[ j ][ i ] = 0.0;
				}
			}
			v[ i ][ i ] = 1.0;
			g = rv1[ i ];
			l = i;
		}
		for ( i = std::min( m, n ) - 1; i >= 0; i-- ) {
			l = i + 1;
			g = w[ i ];
			for ( j = l; j < n; j++ ) {
				u[ i ][ j ] = 0.0;
			}
			if ( g != 0.0 ) {
				g = 1.0 / g;
				for ( j = l; j < n; j++ ) {
					for ( s = 0.0, k = l; k < m; k++ ){
						s += u[k][i]*u[k][j];
					}
					f = ( s / u[ i ][ i ] ) * g;
					for ( k = i; k < m; k++ ) {
						u[ k ][ j ] += f * u[ k ][ i ];
					}
				}
				for ( j = i; j < m; j++ ) {
					u[j][i] *= g;
				}
			} else for ( j = i; j < m; j++ ){
				u[ j ][ i ] = 0.0;
			}
			++u[ i ][ i ];
		}
		for ( k = n - 1; k >= 0; k-- ) {
			for ( its = 0; its < 30; its++ ) {
				flag = true;
				for ( l = k; l >= 0; l-- ) {
					nm = l - 1;
					if ( l == 0 || std::abs( rv1[ l ] ) <= eps * anorm ) {
						flag = false;
						break;
					}
					if ( std::abs( w[ nm ] ) <= eps * anorm ) break;
				}
				if ( flag ) {
					c = 0.0;
					s = 1.0;
					for ( i = l; i < k + 1; i++ ) {
						f = s * rv1[ i ];
						rv1[ i ] = c * rv1[ i ];
						if ( std::abs( f ) <= eps * anorm ) break;
						g = w[ i ];
						h = pythag( f, g );
						w[ i ] = h;
						h = 1.0 / h;
						c = g * h;
						s = - f * h;
						for ( j = 0; j < m; j++ ) {
							y = u[ j ][ nm ];
							z = u[ j ][ i ];
							u[ j ][ nm ] = y * c + z * s;
							u[ j ][ i ] = z * c - y * s;
						}
					}
				}
				z = w[ k ];
				if ( l == k ) {
					if ( z < 0.0 ) {
						w[ k ] = -z;
						for ( j = 0; j < n; j++ ) {
							v[ j ][ k ] = - v[ j ][ k ];
						}
					}
					break;
				}
				if ( its == 29 ){
					throw Error( "SVD error: no convergence in 30 iterations" );
				}
				x = w[ l ];
				nm = k - 1;
				y = w[ nm ];
				g = rv1[ nm ];
				h = rv1[ k ];
				f = ( ( y - z ) * ( y + z ) + ( g - h ) * ( g + h ) ) / ( 2.0 * h * y );
				g = pythag( f, 1.0 );
				f = ( ( x - z ) * ( x + z ) + h * ( ( y / ( f + SIGN( g, f ) ) ) - h ) )
				 		/ x;
				c = s = 1.0;
				for ( j = l; j <= nm; j++ ) {
					i = j + 1;
					g = rv1[ i ];
					y = w[ i ];
					h = s * g;
					g = c * g;
					z = pythag( f, h );
					rv1[ j ] = z;
					c = f / z;
					s = h / z;
					f = x * c + g * s;
					g = g * c - x * s;
					h = y * s;
					y *= c;
					for ( jj = 0; jj < n; jj++ ) {
						x = v[ jj ][ j ];
						z = v[ jj ][ i ];
						v[ jj ][ j ] = x * c + z * s;
						v[ jj ][ i ] = z * c - x * s;
					}
					z = pythag( f, h );
					w[ j ] = z;
					if ( z ) {
						z = 1.0 / z;
						c = f * z;
						s = h * z;
					}
					f = c * g + s * y;
					x = c * y - s * g;
					for ( jj = 0; jj < m; jj++ ) {
						y = u[ jj ][ j ];
						z = u[ jj ][ i ];
						u[ jj ][ j ] = y * c + z * s;
						u[ jj ][ i ] = z * c - y * s;
					}
				}
				rv1[ l ] = 0.0;
				rv1[ k ] = f;
				w[ k ] = x;
			}
		}
	}

	template <>
	inline void SVD< std::complex<double> >::decompose()
	{

	}

	template <typename T>
	inline void SVD<T>::reorder()
	{
		int i, j, k, s, inc = 1;
		T sw;
		Vector<T> su( m ), sv( n );
		do { inc *= 3; inc++; } while ( inc <= n );
		do {
			inc /= 3;
			for ( i = inc; i < n; i++ ) {
				sw = w[ i ];
				for ( k = 0; k < m; k++ ) {
					su[ k ] = u[ k ][ i ];
				}
				for ( k = 0; k < n; k++ ) {
					sv[ k ] = v[ k ][ i ];
				}
				j = i;
				while ( w[ j - inc ] < sw ) {
					w[ j ] = w[ j - inc ];
					for ( k = 0; k < m; k++ ) {
						u[ k ][ j ] = u[ k ][ j - inc ];
					}
					for ( k = 0; k < n; k++ ) {
						v[ k ][ j ] = v[ k ][ j - inc ];
					}
					j -= inc;
					if ( j < inc ) break;
				}
				w[ j ] = sw;
				for ( k = 0; k < m; k++ ) {
					u[ k ][ j ] = su[ k ];
				}
				for ( k = 0; k < n; k++ ) {
					v[ k ][ j ] = sv[ k ];
				}
			}
		} while ( inc > 1 );
		for ( k = 0; k < n; k++ ) {
			s = 0;
			for ( i = 0; i < m; i++ ) {
				if ( u[ i ][ k ] < 0. ) s++;
			}
			for ( j = 0; j < n; j++ ) {
				if ( v[ j ][ k ] < 0. ) s++;
			}
			if ( s > ( m + n ) / 2 ) {
				for ( i = 0; i < m; i++ ) {
					u[ i ][ k ] = - u[ i ][ k ];
				}
				for ( j = 0; j < n; j++ ) {
					v[ j ][ k ] = - v[ j ][ k ];
				}
			}
		}
	}

	template <typename T>
	inline double SVD<T>::pythag( const double a, const double b )
	{
		double absa = std::abs( a ), absb = std::abs( b );
		double ba = absb / absa, ab = absa / absb;
		return ( absa > absb ? absa * std::sqrt( 1.0 + ba * ba ) :
			( absb == 0.0 ? 0.0 : absb * std::sqrt( 1.0 + ab * ab ) ) );
	}

	template <typename T>
  inline double SVD<T>::SIGN( const double &a, const double &b )
  {
    return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
  }

}	// End of namespace Luna

#endif
