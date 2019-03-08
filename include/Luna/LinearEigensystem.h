/// \file LinearEigensystem.h
/// A templated class for solving linear eigenvalue problems of the form
/// A*v=lambda*v where A is a matrix, lambda is an eigenvalue and v is an
/// eigenvector.

#ifndef LINEAREIGENSYSTEM_H
#define LINEAREIGENSYSTEM_H

#include <complex>
#include <limits>
#include <algorithm>

#include "Error.h"
#include "Vector.h"
#include "Matrix.h"


namespace Luna
{

  // Templated linear eigensystem class
  template <class T>

  /// A LinearEigensystem class for use with double and std::complex<double>
  class LinearEigensystem
  {
  private:
    std::size_t N;                             // Number of eigenvalues
    Vector<std::complex<double>> EIGENVALUES;  // Vector for eigenvalues
    Vector<std::complex<double>> ALPHAS;       // Vector for complex numerators
    Vector<double> BETAS;                      // Vector for real denominators
    Matrix<T> EIGENVECTORS;                    // Matrix for eigenvectors
    bool EIGENVALUES_COMPUTED;                 // Eigenvalues calculated?
    bool EIGENVECTORS_COMPUTED;                // Eigenvectors calculated?

    Matrix<T> A;                               // Eigensystem Matrix
    Vector<double> SCALE;                      // Scale factors Vector
    Vector<int> PERM;                          // Permutations

  public:

    /* ----- Constructors and Destructor ----- */

    /// Constructor for a LinearEigensystem
    LinearEigensystem() : EIGENVALUES_COMPUTED( false ),
                          EIGENVECTORS_COMPUTED( false )
    {}

    /// Destructor
    ~LinearEigensystem() {}

    /* ----- Methods ----- */

    //TODO compute_real_symmetric_tridiagonal, tidy up

    /// Solve a LinearEigensystem A*v=lambda*v where A is a real symmetric
    /// Matrix to give the eigenvalues ( and optionally the eigenvectors )
    /// \param A The Matrix A in the eigensystem
    /// \param compute_evecs True if the eigenvectors are to be computed
    void compute_real_symmetric( const Matrix<double>& a,
                                 bool compute_evecs = true );

    /// Solve a LinearEigensystem A*v=lambda*v where A is a real symmetric
    /// tridiagonal Matrix to give the eigenvalues ( and optionally the
    /// eigenvectors )
    /// \param diag The Vector of main diagonal entries in the A Matrix
    /// \param off_diag The Vector of off-diagonal entries in the A Matrix
    /// \param compute_evecs True if the eigenvectors are to be computed
    void compute_real_symmetric_tridiagonal( const Vector<double>& diag,
                                             const Vector<double>& off_diag,
                                             bool compute_evecs = true );

    /// Solve a LinearEigensystem A*v=lambda*v where A is a real Matrix to give
    /// the eigenvalues ( and optionally the eigenvectors )
    /// \param A The Matrix A in the eigensystem
    /// \param compute_evecs True if the eigenvectors are to be computed too
    void compute_real( const Matrix<double>& a, bool compute_evecs = true );

    //TODO compute_real -> eigenvectors not working properly


    //TODO solve_complex

    //TODO solve_generalised_real

    //TODO solve_generalised_complex

    //TODO general solve

    //TODO eigenvalues computed, eigenvectors computed

    //TODO return alphas, return betas 



    /// Return the eigenvalues
    /// \return The eigenvalues as a Vector of std::complex<double>
    Vector<std::complex<double>> eigenvalues() const;

    /// Return a Matrix containing the eigenvectors
    /// \return A Matrix where each column is an eigenvector
    Matrix<T> eigenvector_matrix() const;

    /// Return an std::vector of eigenvectors
    /// \return An std::vector<Vector<std::complex<double>>> where each entry in
    /// the std::vector is an eigenvector
    std::vector<Vector<std::complex<double>>> eigenvectors() const;

  private:

    // Householder reduction of real symmetric matrix to tridiagonal form
    void tred2( Matrix<double>& z, Vector<double>& d, Vector<double>& e,
                bool compute_evecs );

    // QL algorithm with implicit shift to determine the eigenvalues ( and
    // eigenvalues ) of a real symmetric tridiagonal matrix
    void tqli( Matrix<double>& z, Vector<double>& d, Vector<double>& e,
               bool compute_evecs );

    double pythag(const double a, const double b);

    // Sort the eigenvalues (and eigenvectors) into descending order
    void eigsrt( Vector<double> &d, Matrix<double> *v=NULL);
    void sort( Vector<double>& d, Matrix<double>& z, bool evecs_computed );



    //TODO need to sort out functions below here

    void balance();

    void elmhes();

    void hqr();

    T sign( const T &a, const T &b );

    //void sort();

    void eltran();

    void hqr2();

    void balbak()
    {
      for ( int i = 0; i < N; i++ )
    		for ( int j = 0; j < N; j++ )
    			EIGENVECTORS( i, j ) *= SCALE[ i ];
    }

    void sortvecs()
    {
    	int i;
    	Vector<T> temp( N );
    	for (int j=1;j<N;j++) {
    		std::complex<double> x = EIGENVALUES[ j ];
    		for (int k=0;k<N;k++)
    			temp[k]=EIGENVECTORS( k, j );
    		for (i=j-1;i>=0;i--) {
    			if (std::real(EIGENVALUES[i]) >= std::real(x)) break;
    			EIGENVALUES[i+1]=EIGENVALUES[i];
    			for (int k=0;k<N;k++)
    				EIGENVECTORS( k, i + 1 )=EIGENVECTORS( k, i );
    		}
    		EIGENVALUES[i+1]=x;
    		for ( int k=0;k<N;k++)
    			EIGENVECTORS( k, i + 1 ) = temp[ k ];
    	}
    }

    //TODO elmhes etc (possibly rename) - brief descriptions

    void find_small_subdiagonals( int& l, int& nn, double& s,
                                  const double& anorm )
    {
      const double EPS = std::numeric_limits<double>::epsilon();
      for ( l = nn; l > 0; l-- ) {
        s = std::abs( A( l - 1, l - 1 ) ) + std::abs( A( l, l ) );
        if ( s == 0.0 ) s = anorm;
        if ( std::abs( A( l, l - 1 ) ) <= EPS * s ) {
          A( l, l - 1 ) = 0.0;
          break;
        }
      }
    }

    void form_shift( int& i, int& m, int& nn, int& l, double& p, double& q,
                     double& r, double& s, double& u, double& v, double& w,
                     double& x, double& y, double& z )
    {
      const double EPS = std::numeric_limits<double>::epsilon();
      for ( m = nn - 2; m >= l; m-- ) {
        z = A( m, m );
        r = x - z;
        s = y - z;
        p = ( r * s - w ) / A( m + 1, m ) + A( m, m + 1 );
        q = A( m + 1, m + 1 ) - z - r - s;
        r = A( m + 2, m + 1 );
        s = std::abs( p ) + std::abs( q ) + std::abs( r );
        p /= s;
        q /= s;
        r /= s;
        if (m == l) break;
        u = std::abs( A( m, m - 1 ) ) * ( std::abs( q ) + std::abs( r ) );
        v = std::abs( p ) * ( std::abs( A( m - 1, m - 1 ) ) +
            std::abs( z ) + std::abs( A( m + 1, m + 1 ) ) );
        if ( u <= EPS * v ) break;
      }
      for ( i = m; i < nn - 1; i++ ) {
        A( i + 2, i ) = 0.0;
        if ( i != m ) A( i + 2, i - 1 ) = 0.0;
      }
    }

    void form_exceptional_shift( int& i, int& nn, double& s, double& t,
                                 double& w, double& x, double& y )
    {
      t += x;
      for ( i = 0; i < nn + 1; i++ ) A( i, i ) -= x;
      s = std::abs( A( nn, nn - 1) ) + std::abs( A( nn - 1, nn - 2 ) );
      y = x = 0.75 * s;
      w = -0.4375 * s * s;
    }

    void setup_householder( int& k, int& nn, double& p, double& q, double& r,
                            double& x )
    {
      p = A( k, k - 1 );
      q = A( k + 1, k - 1 );
      r = 0.0;
      if ( k + 1 != nn ) r = A( k + 2, k - 1 );
      if ( ( x = std::abs(p) + std::abs(q) + std::abs(r) ) != 0.0 ) {
        p /= x;
        q /= x;
        r /= x;
      }
    }

    void row_modification( int& j, int& k, int& nn, double& p, double& q,
                   double& r, double& s, double& x, double& y, double& z )
    {
      p += s;
      x = p / s;
      y = q / s;
      z = r / s;
      q /= p;
      r /= p;
      for ( j = k; j < nn + 1; j++ ) {
        p = A( k, j ) + q * A( k + 1, j );
        if ( k + 1 != nn ) {
          p += r * A( k + 2, j );
          A( k + 2, j ) -= p * z;
        }
        A( k + 1, j ) -= p * y;
        A( k, j ) -= p * x;
      }
    }

    void column_modification( int& i, int& l, int& k, int& nn, int& mmin,
            double& p, double& q, double& r, double& x, double& y, double& z )
    {
      for ( i = l; i < mmin + 1; i++ ) {
        p = x * A( i, k ) + y * A( i, k + 1 );
        if ( k + 1 != nn ) {
          p += z * A( i, k + 2 );
          A( i, k + 2 ) -= p * r;
        }
        A( i, k + 1 ) -= p * q;
        A( i, k ) -= p;
      }
    }


  };	// End of class LinearEigensystem

  template <>
  inline void LinearEigensystem<double>::compute_real_symmetric(
                                  const Matrix<double>& a, bool compute_evecs )
  {
    if ( a.rows() != a.cols() )
    {
      throw Error( "compute_real_symmetric error: Matrix is not square.");
    }
    Matrix<double> z( a );
    Vector<double> diag, off_diag;
    tred2( z, diag, off_diag, compute_evecs );
    tqli( z, diag, off_diag, compute_evecs );
    sort( diag, z, compute_evecs );

    std::size_t n( z.rows() );
    EIGENVALUES.resize( n );
    EIGENVECTORS.resize( n, n );
    for ( std::size_t i = 0; i < n; i++ )
    {
      EIGENVALUES[ i ] = diag[ i ];
    }
    EIGENVALUES_COMPUTED = true;
    if ( compute_evecs )
    {
      EIGENVECTORS = z;
      EIGENVECTORS_COMPUTED = true;
    }
  }

  template <>
  inline void LinearEigensystem<double>::compute_real_symmetric_tridiagonal(
              const Vector<double>& diag, const Vector<double>& off_diag,
              bool compute_evecs )
  {
    if ( diag.size() != off_diag.size() + 1 )
    {
      std::string problem;
      problem = "compute_real_symmetric_tridiagonal error: ";
      problem += "size of Vectors are incompatible.";
      throw Error( problem );
    }
    std::size_t n( diag.size() );
    Matrix<double> z( n, n, 0.0 );
    z.eye();

    Vector<double> diag_temp( diag ), off_diag_temp( off_diag );
    off_diag_temp.push_front( 0.0 );
    tqli( z, diag_temp, off_diag_temp, compute_evecs );
    sort( diag_temp, z, compute_evecs );

    EIGENVALUES.resize( n );
    EIGENVECTORS.resize( n, n );
    for ( std::size_t i = 0; i < n; i++ )
    {
      EIGENVALUES[ i ] = diag_temp[ i ];
    }
    EIGENVALUES_COMPUTED = true;
    if ( compute_evecs )
    {
      EIGENVECTORS = z;
      EIGENVECTORS_COMPUTED = true;
    }
  }

  template <>
  inline void LinearEigensystem<double>::compute_real( const Matrix<double>& a,
                                                       bool compute_evecs)
  {
    //TODO check A is square
    A = a;
    N = a.rows();

    balance();
    elmhes();
    std::cout << "a = " << A << std::endl;
    if ( compute_evecs ) {
      EIGENVECTORS.resize( N, N );
      EIGENVECTORS.fill( 1.0 );
      eltran();
      hqr2();
      balbak();
      sortvecs();

      EIGENVALUES_COMPUTED = true;
      EIGENVECTORS_COMPUTED = true;
    } else {
      hqr();
      //sort();
      EIGENVALUES_COMPUTED = true;
    }



  }



  template <typename T>
  inline Vector<std::complex<double>> LinearEigensystem<T>::eigenvalues() const
  {
    if ( EIGENVALUES_COMPUTED ) { return EIGENVALUES; }
    else { throw Error( "Eigensystem: eigenvalues not computed." ); }
  }

  template <typename T>
  inline Matrix<T> LinearEigensystem<T>::eigenvector_matrix() const
  {
    if ( EIGENVECTORS_COMPUTED ) { return EIGENVECTORS; }
    else { throw Error( "Eigensystem: eigenvectors not computed." ); }
  }

  typedef std::complex<double> cmplx;
  template <typename T>
  inline std::vector<Vector<cmplx>> LinearEigensystem<T>::eigenvectors() const
  {
    if ( EIGENVECTORS_COMPUTED )
    {
      std::vector< Vector< std::complex<double> > > evecs;
      std::size_t rows = EIGENVECTORS.rows();
      std::size_t cols = EIGENVECTORS.cols();

      for (std::size_t j = 0; j < cols; ++j )
      {
          Vector< std::complex<double> > evec;
          for (std::size_t i=0; i<rows; ++i)
          {
              evec.push_back( EIGENVECTORS( i, j ) );
          }
          evecs.push_back( evec );
      }
      return evecs;
    }
    else { throw Error( "Eigensystem: eigenvectors not computed." ); }
  }

  /* ----- Private functions ----- */

  template <typename T>
  inline void LinearEigensystem<T>::tred2( Matrix<double>& z, Vector<double>& d,
                                         Vector<double>& e, bool compute_evecs )
  {
    int i, j, k, l, n( z.rows() );
    d.resize( n );
    e.resize( n );
    double scale, hh, h, g, f;
    for ( i = n - 1; i > 0; i-- ) {
      l = i - 1;
      h = scale = 0.0;
      if (l > 0) {
        for ( k = 0; k < i; k++ ) {
          scale += std::abs( z( i, k ) );
        }
        if ( scale == 0.0 )
          e[ i ] = z( i, l );
        else {
          for ( k = 0; k < i; k++ ) {
            z( i, k ) /= scale;
            h += z( i, k ) * z( i, k );
          }
          f = z( i, l );
          g = ( f >= 0.0 ? -std::sqrt( h ) : std::sqrt( h ) );
          e[ i ] = scale * g;
          h -= f * g;
          z( i, l ) = f - g;
          f = 0.0;
          for ( j = 0; j < i; j++ ) {
            if ( compute_evecs ){
              z( j, i ) = z( i, j ) / h;
            }
            g = 0.0;
            for ( k = 0; k < j + 1; k++ )
              g += z( j, k ) * z( i, k );
            for ( k = j + 1; k < i; k++ )
              g += z( k, j ) * z( i, k );
            e[ j ] = g / h;
            f += e[ j ] * z( i, j );
          }
          hh = f / ( h + h );
          for ( j = 0; j < i; j++ ) {
            f = z( i, j );
            e[ j ] = g = e[ j ] - hh * f;
            for ( k = 0; k < j + 1; k++ ){
              z( j, k ) -= ( f * e[ k ] + g * z( i, k ) );
            }
          }
        }
      } else
        e[ i ] = z( i, l );
      d[ i ] = h;
    }
    if ( compute_evecs ){
       d[ 0 ] = 0.0;
    }
    e[ 0 ] = 0.0;
    for ( i = 0; i < n; i++ ) {
      if ( compute_evecs ) {
        if ( d[ i ] != 0.0 ) {
          for ( j = 0; j < i; j++ ) {
            g = 0.0;
            for ( k = 0; k < i; k++ ){
              g += z( i, k ) * z( k, j );
            }
            for (k = 0; k < i; k++ ) {
              z( k, j ) -= g * z( k, i );
            }
          }
        }
        d[ i ] = z( i, i );
        z( i, i ) = 1.0;
        for ( j = 0; j < i; j++ ){
           z( j, i ) = z( i, j ) = 0.0;
        }
      } else {
        d[ i ] = z( i, i );
      }
    }
  }

  template <typename T>
  inline void LinearEigensystem<T>::tqli( Matrix<double>& z, Vector<double>& d,
                                        Vector<double>& e, bool compute_evecs )
  {
    int m, l, iter, i, k, n( z.rows() );
    double s, r, p, g, f, dd, c, b;
    const double EPS = std::numeric_limits<double>::epsilon();
    for ( i = 1; i < n; i++ ){
       e[ i - 1 ] = e[ i ];
    }
    e[ n - 1 ]=0.0;
    for ( l = 0; l < n; l++ ) {
      iter = 0;
      do {
        for ( m = l; m < n - 1; m++ ) {
          dd = std::abs( d[ m ] ) + std::abs( d[ m + 1 ] );
          if ( std::abs( e[ m ] ) <= EPS * dd ) break;
        }
        if (m != l) {
          if ( iter++ == 30 ){ throw Error("Too many iterations in tqli."); }
          g = ( d[ l + 1 ] - d[ l ] ) / ( 2.0 * e[ l ] );
          r = this->pythag( g, 1.0 );
          g = d[ m ] - d[ l ] + e[ l ] / ( g + this->sign( r, g ) );
          s = c = 1.0;
          p = 0.0;
          for ( i = m - 1; i >= l; i-- ) {
            f = s * e[ i ];
            b = c * e[ i ];
            e[ i + 1 ] = ( r = this->pythag( f, g ) );
            if ( r == 0.0 ) {
              d[ i + 1 ] -= p;
              e[ m ] = 0.0;
              break;
            }
            s = f / r;
            c = g / r;
            g = d[ i + 1 ] - p;
            r = ( d[ i ] - g ) * s + 2.0 * c * b;
            d[ i + 1 ] = g + ( p = s * r );
            g = c * r - b;
            if ( compute_evecs ) {
              for ( k = 0; k < n; k++ ) {
                f = z( k, i + 1 );
                z( k, i + 1 ) = s * z( k, i ) + c * f;
                z( k, i ) = c * z( k, i ) - s * f;
              }
            }
          }
          if ( r == 0.0 && i >= l ) continue;
          d[ l ] -= p;
          e[ l ] = g;
          e[ m ] = 0.0;
        }
      } while (m != l);
    }
  }

  template <typename T>
  inline double LinearEigensystem<T>::pythag(const double a, const double b)
  {
    double absa = std::abs(a), absb = std::abs(b);
    return (absa > absb ? absa * std::sqrt( 1. + (absb/absa)*(absb/absa) ) :
      (absb == 0.0 ? 0.0 : absb * std::sqrt( 1. + (absa/absb)*(absa/absb) )));
  }

  template <typename T>
  inline void LinearEigensystem<T>::eigsrt( Vector<double> &d,
                                            Matrix<double> *v )
  {
    int k;
    int n = d.size();
    for ( int i = 0; i < n - 1; i++ ) {
      double p = d[ k = i ];
      for (int j = i; j < n; j++ ){
        if ( d[ j ] >= p ){
           p = d[ k = j ];
        }
      }
      if ( k != i ) {
        d[ k ] = d[ i ];
        d[ i ] = p;
        if ( v != NULL ) {
          for ( int j = 0; j < n; j++ ) {
            p = (*v)( j, i );
            (*v)( j, i ) = (*v)( j, k );
            (*v)( j, k ) = p;
          }
        }
      }
    }
  }

  template <typename T>
  inline void LinearEigensystem<T>::sort( Vector<double>& d, Matrix<double>& z,
                                          bool evecs_computed ) {
    if ( evecs_computed ){
      eigsrt(d,&z);
    } else {
      eigsrt(d);
    }
  }

  template <typename T>
  inline void LinearEigensystem<T>::balance()
  {
    SCALE.assign( N, 1.0 );
    const double RADIX = std::numeric_limits<double>::radix;
    bool done( false );
    double sqrdx( RADIX * RADIX );
    while (!done) {
      done=true;
      for ( int i = 0; i < N; i++ ) {
        double r = 0.0, c = 0.0;
        for ( int j = 0; j < N; j++ )
          if ( j != i ) {
            c += abs( A( j, i ) );
            r += abs( A( i, j ) );
          }
        if ( c != 0.0 && r != 0.0 ) {
          double g = r / RADIX;
          double f = 1.0;
          double s = c + r;
          while ( c < g ) {
            f *= RADIX;
            c *= sqrdx;
          }
          g = r * RADIX;
          while ( c > g ) {
            f /= RADIX;
            c /= sqrdx;
          }
          if ( ( c + r ) / f < 0.95 * s ) {
            done = false;
            g = 1.0 / f;
            SCALE[ i ] *= f;
            for ( int j = 0; j < N; j++ ) A( i, j ) *= g;
            for ( int j = 0; j < N; j++ ) A( j, i ) *= f;
          }
        }
      }
    }
  }

  template <typename T>
  inline void LinearEigensystem<T>::elmhes()
  {
    PERM.resize( N );
    for ( int m = 1; m < N - 1; m++ ) {
      double x = 0.0;
      int i = m;
      for ( int j = m; j < N; j++ ) {
        if ( std::abs( A( j, m - 1 ) ) > std::abs( x ) ) {
          x = A( j, m - 1 );
          i = j;
        }
      }
      PERM[ m ] = i;
      if ( i != m ) {
        for ( int j = m - 1; j < N; j++ ) std::swap( A( i, j ) , A( m, j ) );
        for ( int j = 0; j < N; j++ ) std::swap( A( j, i ), A( j, m ) );
      }
      if ( x != 0.0 ) {
        for ( i = m + 1; i < N; i++ ) {
          double y = A( i, m - 1 );
          if ( y != 0.0 ) {
            y /= x;
            A( i, m - 1 ) = y;
            for ( int j = m; j < N; j++ ) A( i, j ) -= y * A( m, j );
            for ( int j = 0; j < N; j++ ) A( j, m ) += y * A( j, i );
          }
        }
      }
    }
  }

  template <typename T>
  inline void LinearEigensystem<T>::hqr()
  {
  	int nn, m, l, k, j, its, i, mmin;
  	double z, y, x, w, v, u, t, s, r, q, p, anorm = 0.0;
    EIGENVALUES.resize( N );

  	const double EPS = std::numeric_limits<double>::epsilon();
    anorm = A.norm_1();

  	nn = N - 1;
  	t = 0.0;
  	while ( nn >= 0 ) {
  		its = 0;
  		do {
        find_small_subdiagonals( l, nn, s, anorm );
  			x = A( nn, nn );
  			if ( l == nn ) { //one root found
  				EIGENVALUES[ nn-- ] = x + t;
  			} else {
  				y = A( nn - 1, nn - 1 );
  				w = A( nn, nn - 1 ) * A( nn - 1, nn );
  				if ( l == nn - 1 ) { //two roots found
  					p = 0.5 * ( y - x );
  					q = p * p + w;
  					z = std::sqrt( std::abs( q ) );
  					x += t;
  					if ( q >= 0.0 ) {
  						z = p + sign( z, p );
  						EIGENVALUES[ nn - 1 ] = EIGENVALUES[ nn ] = x + z;
  						if ( z != 0.0 ) EIGENVALUES[ nn ] = x - w / z;
  					} else {
  						EIGENVALUES[ nn ] = std::complex<double>( x + p, - z );
  						EIGENVALUES[ nn - 1 ] = std::conj( EIGENVALUES[ nn ] );
  					}
  					nn -= 2;
  				} else { //no roots found
  					if ( its == 30 ){ throw Error( "Too many iterations in hqr." ); }
  					if ( its == 10 || its == 20 ) {
              form_exceptional_shift( i, nn, s, t, w, x, y );
  					}
  					++its;
            form_shift( i, m, nn, l, p, q, r, s, u, v, w, x, y, z );

  					for ( k = m; k < nn; k++ ) {
  						if ( k != m ) {
                setup_householder( k, nn, p, q, r, x );
  						}
  						if ( ( s = sign( std::sqrt( p*p + q*q + r*r ), p ) ) != 0.0 ) {
  							if ( k == m ) {
  								if ( l != m )
  								A( k, k - 1 ) = - A( k, k - 1 );
  							} else {
  								A( k, k - 1 ) = - s * x;
                }

                row_modification( j, k, nn, p, q, r, s, x, y, z );
  							mmin = nn < k + 3 ? nn : k + 3;
                column_modification( i, l, k, nn, mmin, p, q, r, x, y, z );
  						}
  					}
  				}
  			}
  		} while ( l + 1 < nn );
  	}
  }

  template <typename T>
  inline T LinearEigensystem<T>::sign( const T &a, const T &b )
  {
    return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
  }

  /*template <typename T>
  inline void LinearEigensystem<T>::sort()
  {
    int i;
    for ( int j = 1; j < N; j++ ) {
      std::complex<double> x = EIGENVALUES[ j ];
      for ( i = j - 1; i >= 0; i-- ) {
        if ( std::real( EIGENVALUES[ i ] ) >= std::real( x ) ) break;
        EIGENVALUES[ i + 1 ] = EIGENVALUES[ i ];
      }
      EIGENVALUES[ i + 1 ] = x;
    }
  }*/

  template <typename T>
  inline void LinearEigensystem<T>::eltran()
  {
    for ( int mp = N - 2; mp > 0; mp-- ) {
      for ( int k = mp + 1; k < N; k++ )
        EIGENVECTORS( k, mp ) = A( k, mp - 1 );
      int i = PERM[ mp ];
      if ( i != mp ) {
        for ( int j = mp; j < N; j++ ) {
          EIGENVECTORS( mp, j ) = EIGENVECTORS( i, j );
          EIGENVECTORS( i, j ) = 0.0;
        }
        EIGENVECTORS( i, mp ) = 1.0;
      }
    }
  }

  template <typename T>
  inline void LinearEigensystem<T>::hqr2()
  {
    int nn, m, l, k, j, its, i, mmin, na;
  	double z, y, x, w, v, u, t, s, r, q, p, anorm=0.0, ra, sa, vr, vi;
    typedef std::complex<double> cmplx;
    EIGENVALUES.resize( N );

  	const double EPS = std::numeric_limits<double>::epsilon();
    anorm = A.norm_1();

  	nn = N - 1;
  	t = 0.0;
  	while (nn >= 0) {
  		its=0;
  		do {
        find_small_subdiagonals( l, nn, s, anorm );
  			x = A( nn, nn );
  			if (l == nn) {
  				EIGENVALUES[ nn ] = A( nn, nn ) = x + t;
  				nn--;
  			} else {
  				y = A( nn - 1, nn - 1 );
  				w = A( nn, nn - 1 ) * A( nn - 1, nn );
  				if ( l == nn - 1 ) {
  					p = 0.5 * ( y - x );
  					q = p * p + w;
  					z = std::sqrt( std::abs( q ) );
  					x += t;
  					A( nn, nn ) = x;
  					A( nn - 1, nn - 1 ) = y + t;
  					if ( q >= 0.0 ) {
  						z = p + sign( z, p );
  						EIGENVALUES[ nn - 1 ] = EIGENVALUES[ nn ] = x + z;
  						if ( z != 0.0 ) EIGENVALUES[ nn ] = x - w / z;
  						x = A( nn, nn - 1 );
  						s = std::abs( x ) + std::abs( z );
  						p = x / s;
  						q = z / s;
  						r = std::sqrt( p * p + q * q );
  						p /= r;
  						q /= r;
  						for ( j = nn - 1; j < N; j++ ) {
  							z = A( nn - 1, j );
  							A( nn - 1, j ) = q * z + p * A( nn, j );
  							A( nn, j ) = q * A( nn, j ) - p * z;
  						}
  						for ( i = 0; i <= nn; i++ ) {
  							z = A( i, nn - 1 );
  							A( i, nn - 1 ) = q * z + p * A( i, nn );
  							A( i, nn ) = q * A( i, nn ) - p * z;
  						}
  						for ( i = 0; i < N; i++ ) {
  							z = EIGENVECTORS( i, nn - 1 );
  							EIGENVECTORS( i, nn - 1 ) = q  * z + p * EIGENVECTORS( i, nn );
  							EIGENVECTORS( i, nn ) = q * EIGENVECTORS( i, nn ) - p * z;
  						}
  					} else {
  						EIGENVALUES[ nn ] = std::complex<double>( x + p, - z );
  						EIGENVALUES[ nn - 1 ] = std::conj( EIGENVALUES[ nn ] );
  					}
  					nn -= 2;
  				} else {
  					if ( its == 30 ) throw Error( "Too many iterations in hqr2" );
  					if ( its == 10 || its == 20 ) {
  						form_exceptional_shift( i, nn, s, t, w, x, y );
  					}
  					++its;
  					form_shift( i, m, nn, l, p, q, r, s, u, v, w, x, y, z );

  					for ( k = m; k < nn; k++ ) {
  						if ( k != m ) {
  							p = A( k, k - 1 );
  							q = A( k + 1, k - 1 );
  							r = 0.0;
  							if ( k + 1 != nn ) r = A( k + 2, k - 1 );
  							if ( ( x = std::abs(p) + std::abs(q) + std::abs(r) ) != 0.0 ) {
  								p /= x;
  								q /= x;
  								r /= x;
  							}
  						}
  						if ( ( s = sign( std::sqrt( p*p + q*q + r*r ), p ) ) != 0.0 ) {
  							if ( k == m ) {
  								if ( l != m )
  								A( k, k - 1 ) = - A( k, k - 1 );
  							} else
  								A( k, k - 1 ) = - s * x;

                row_modification( j, k, nn, p, q, r, s, x, y, z );
  							mmin = nn < k+3 ? nn : k+3;
  							column_modification( i, l, k, nn, mmin, p, q, r, x, y, z );

                for ( i = 0; i < N; i++ ) {
  								p = x * EIGENVECTORS( i, k ) + y * EIGENVECTORS( i, k + 1 );
  								if ( k + 1 != nn ) {
  									p += z * EIGENVECTORS( i, k + 2 );
  									EIGENVECTORS( i, k + 2 ) -= p * r;
  								}
  								EIGENVECTORS( i, k + 1 ) -= p * q;
  								EIGENVECTORS( i, k ) -= p;
  							}
  						}
  					}
  				}
  			}
  		} while ( l + 1 < nn );
  	}
  	if (anorm != 0.0) {
  		for ( nn = N - 1; nn >= 0; nn-- ) {
  			p = std::real( EIGENVALUES[ nn ] );
  			q = std::imag( EIGENVALUES[ nn ] );
  			na = nn - 1;
  			if ( q == 0.0 ) {
  				m = nn;
  				A( nn, nn ) = 1.0;
  				for ( i = nn - 1; i >= 0; i-- ) {
  					w = A( i, i ) - p;
  					r = 0.0;
  					for ( j = m; j <= nn; j++ )
  						r += A( i, j ) * A( j, nn );
  					if ( std::imag( EIGENVALUES[ i ] ) < 0.0 ) {
  						z = w;
  						s = r;
  					} else {
  						m = i;

  						if ( std::imag( EIGENVALUES[ i ] ) == 0.0 ) {
  							t = w;
  							if ( t == 0.0 )
  								t = EPS * anorm;
  							A( i, nn ) = - r / t;
  						} else {
  							x = A( i, i + 1 );
  							y = A( i + 1, i );
  							q = std::pow( std::real( EIGENVALUES[ i ] ) - p, 2 ) +
                    std::pow( std::imag( EIGENVALUES[ i ] ), 2 );
  							t = ( x * s - z * r ) / q;
  							A( i, nn ) = t;
  							if ( std::abs( x ) > std::abs( z ) )
  								A( i + 1, nn ) = ( - r - w * t ) / x;
  							else
  								A( i + 1, nn ) = ( - s - y * t ) / z;
  						}
  						t = std::abs( A( i, nn ) );
  						if ( EPS * t * t > 1 )
  							for ( j = i; j <= nn; j++ )
  								A( j, nn ) /= t;
  					}
  				}
  			} else if ( q < 0.0 ) {
  				m = na;
  				if ( std::abs( A( nn, na ) ) > std::abs( A( na, nn ) ) ) {
  					A( na, na ) = q / A( nn, na );
  					A( na, nn ) = - ( A( nn, nn ) - p ) / A( nn, na );
  				} else {
  					cmplx temp = cmplx( 0.0, - A( na, nn ) ) /
                         cmplx( A( na, na ) - p , q );
  					A( na, na ) = std::real( temp );
  					A( na, nn ) = std::imag( temp );
  				}
  				A( nn, na ) = 0.0;
  				A( nn, nn ) = 1.0;
  				for ( i = nn - 2; i >= 0; i-- ) {
  					w = A( i, i ) - p;
  					ra = sa = 0.0;
  					for ( j = m; j <= nn; j++ ) {
  						ra += A( i, j ) * A( j, na );
  						sa += A( i, j ) * A( j, nn );
  					}
  					if ( std::imag( EIGENVALUES[ i ] ) < 0.0 ) {
  						z = w;
  						r = ra;
  						s = sa;
  					} else {
  						m = i;
  						if ( std::imag( EIGENVALUES[ i ] ) == 0.0 ) {
  							cmplx temp = cmplx( - ra, - sa ) / cmplx( w, q );
  							A( i, na ) = std::real( temp );
  							A( i, nn ) = std::imag( temp );
  						} else {
  							x = A( i, i + 1 );
  							y = A( i + 1, i );
  							vr = std::pow( std::real( EIGENVALUES[ i ] ) - p, 2 ) +
                     std::pow( std::imag( EIGENVALUES[ i ] ), 2 ) - q * q;
  							vi = 2.0 * q * ( std::real( EIGENVALUES[ i ] ) - p );
  							if ( vr == 0.0 && vi == 0.0 )
  								vr = EPS * anorm * ( std::abs(w) + std::abs(q) + std::abs(x) +
                       std::abs(y) + std::abs(z) );
  							cmplx temp = cmplx( x * r - z * ra + q * sa,
                                    x * s - z * sa - q * ra ) / cmplx( vr, vi );
  							A( i, na ) = std::real( temp );
  							A( i, nn ) = std::imag( temp );
  							if ( std::abs(x) > std::abs(z) + std::abs(q) ) {
  								A( i + 1, na ) = ( - ra - w * A( i, na ) + q * A( i, nn ) )
                                   / x;
  								A( i + 1, nn ) = ( - sa - w * A( i, nn ) - q * A( i, na ) )
                                   / x;
  							} else {
  								cmplx temp = cmplx( - r - y * A( i, na ),
                                      - s - y * A( i, nn ) ) / cmplx( z, q );
  								A( i + 1, na ) = std::real( temp );
  								A( i + 1, nn ) = std::imag( temp );
  							}
  						}
  					}
  					t = std::max( std::abs( A( i, na ) ), std::abs( A( i, nn ) ) );
  					if ( EPS * t * t > 1 )
  						for ( j = i; j <= nn; j++ ) {
  							A( j, na ) /= t;
  							A( j, nn ) /= t;
  						}
  				}
  			}
  		}
  		for ( j = N-1; j >= 0; j-- )
  			for ( i = 0; i < N; i++ ) {
  				z = 0.0;
  				for ( k = 0; k <= j; k++ )
  					z += EIGENVECTORS( i, k ) * A( k, j );
  				EIGENVECTORS( i, j ) = z;
  			}
  	}
  }


}  // End of namespace Luna

#endif
