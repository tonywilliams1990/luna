/// \file Special.h
/// A file specifying some special mathematical functions not contained
/// in the standard library

#ifndef SPECIAL_H
#define SPECIAL_H

#include <vector>
#include <fstream>
#include <iostream>
#include <complex>
#include <limits>

#include "Error.h"
#include "Vector.h"


namespace Luna
{
    /// The gamma function \f[ \Gamma(z) = \int_0^{\infty} x^{z-1}e^{-x} dx \f]
    /// \param z The argument of the function (not a negative integer or zero)
    /// \return The gamma function \f$\Gamma(z)\f$ using Spouge's approximation
    template <typename T>
    T gamma( const T& z )
    {
      const int a = 12;
      Vector<T> c;
      c.resize( a );
      int k;
      T accm;

      double factrl = 1.0; // (k - 1)!*(-1)^k with 0!==1
      c[ 0 ] = std::sqrt( 2.0 * M_PI );
      for(k=1; k < a; k++) {
        c[ k ] =    std::exp( 1. * a - k )
                  * std::pow( 1. * a - k, k - 0.5 ) / factrl;
  	    factrl *= -k;
      }
      accm = c[ 0 ];
      for( k = 1; k < a; k++ )
      {
        accm += c[ k ] / ( z + 1. * k );
      }
      accm *= std::exp( - ( z + 1. * a ) ) * std::pow( z + 1. * a, z + 0.5 );
      return accm / z; // Gamma(z) = Gamma(z+1) / z
    }

    /// The log-gamma function
    /// \param z The argument of the function \f$ (z > 0) \f$
    /// \return The log of the gamma function \f$\ln \Gamma(z) \f$
    template <typename T>
    T lngamma(const T& z) {
    	T x, tmp, y, ser;
    	static const double cof[ 14 ]={ 57.1562356658629235, -59.5979603554754912,
    	14.1360979747417471, -0.491913816097620199, .339946499848118887e-4,
    	.465236289270485756e-4, -.983744753048795646e-4, .158088703224912494e-3,
    	-.210264441724104883e-3, .217439618115212643e-3, -.164318106536763890e-3,
    	.844182239838527433e-4, -.261908384015814087e-4, .368991826595316234e-5 };
    	if ( std::real( z ) <= 0) throw Error( " lngamma: invalid argument " );
    	y = x = z;
    	tmp = x + 5.24218750000000000;
    	tmp = ( x + 0.5 ) * log( tmp ) - tmp;
    	ser = 0.999999999999997092;
    	for ( std::size_t j = 0;j < 14; j++ )
      {
        y += 1.0;
        ser += cof[ j ] / y;
      }
    	return tmp + std::log( 2.5066282746310005 * ser / x );
    }

    /// The factorial function \f$ n! \f$
    /// \param n The arguement of the function
    /// \return The factorial of n as a double precision floating-point number
    double factorial( const std::size_t& n )
    {
      // Create static table on the first call then lookup on subsequent calls
      static std::vector<double> a( 171 );
      static bool init( true );
      if ( init ) {
        init = false;
        a[ 0 ] = 1.;
        for ( std::size_t i = 1; i < 171; i++ )
        {
          a[ i ] = i * a[ i - 1 ];
        }
      }
      if ( n < 0 || n > 170 ){ throw Error( " Range error in factorial" ); }
      return a[ n ];
    }

    /// The Log-Factorial function \f$ \ln(n!) \f$
    /// \param n The arguement of the function
    /// \return The log-factorial as a double precision floating-point number
    double lnfactorial( const std::size_t& n )
    {
      static const std::size_t NTOP( 2000 );
      static std::vector<double> a( NTOP );
      static bool init( true );
      if ( init ) {
        init = false;
        for ( std::size_t i = 0; i < NTOP; i++ )
        {
          a[ i ] = lngamma( i + 1. );
        }
      }
      if ( n < 0 ){ throw Error( " Negative arguement in lnfactorial" ); }
      if ( n < NTOP ){ return a[ n ]; }
      return lngamma( n + 1. ); // Out of range of table
    }

    /// The beta function \f[ B(z,w)=\frac{\Gamma(z)\Gamma(w)}{\Gamma(z+w)} \f]
    /// \param z The first argument of the function \f$ (Re(z) > 0) \f$
    /// \param w The second argument of the function \f$ (Re(w) > 0) \f$
    /// \return The beta function \f$ B(z,w) \f$
    template <typename T>
    T beta( const T& z, const T& w )
    {
      return std::exp( lngamma(z) + lngamma(w) - lngamma(z+w) );
    }

    /// Complementary error function Chebyshev approximation
    /// \param z The arguement of the function
    /// \return The Chebyshev approximation to the complementary error function
    template <typename T>
    T erfccheb( const T& z )
    {
      int ncof( 28 );
      const double cof[28] = {-1.3026537197817094, 6.4196979235649026e-1,
        1.9476473204185836e-2,-9.561514786808631e-3,-9.46595344482036e-4,
        3.66839497852761e-4,4.2523324806907e-5,-2.0278578112534e-5,
        -1.624290004647e-6,1.303655835580e-6,1.5626441722e-8,-8.5238095915e-8,
        6.529054439e-9,5.059343495e-9,-9.91364156e-10,-2.27365122e-10,
        9.6467911e-11, 2.394038e-12,-6.886027e-12,8.94487e-13, 3.13092e-13,
        -1.12708e-13,3.81e-16,7.106e-15,-1.523e-15,-9.4e-17,1.21e-16,-2.8e-17};
      T t, ty, tmp, d = 0., dd = 0.;
      if ( std::real( z ) < 0. )
      {
        throw Error(" Nonnegative argument required in erfccheb");
      }
      t = 2. / ( 2. + z );
      ty = 4. * t - 2.;
      for ( std::size_t j = ncof - 1; j > 0; j-- )
      {
        tmp = d;
        d = ty * d - dd + cof[ j ];
        dd = tmp;
      }
      return t * exp( - z * z + 0.5 * ( cof[ 0 ] + ty * d ) - dd );
    }

    /// The error function \f[ erf(z) = \frac{1}{\sqrt{\pi}}
    /// \int_{-z}^{z} e^{-t^2} dt \f]
    /// \param z The arguement of the function
    /// \return The error function \f$ erf(z) \f$
    template <typename T>
    T erf( const T& z )
    {
      if ( std::real( z ) >= 0. ) return 1.0 - erfccheb( z );
      else return erfccheb( - z ) - 1.0;
    }

    /// The complementary error function \f[ erfc(z) = 1 - erf(z) \f]
    /// \param z The arguement of the function
    /// \return The complementary error function \f$ erfc(z) \f$
    template <typename T>
    T erfc( const T& z )
    {
      if ( std::real( z ) >= 0. ) return erfccheb( z );
      else return 2.0 - erfccheb( - z );
    }

    template <typename T>
    struct Bessel {
      const double xj00,xj10,xj01,xj11,twoopi,pio4;
      const double j0r[7]={1.682397144220462e-4,2.058861258868952e-5,
        5.288947320067750e-7,5.557173907680151e-9,2.865540042042604e-11,
        7.398972674152181e-14,7.925088479679688e-17};
      const double j0s[7]={1.0,1.019685405805929e-2,5.130296867064666e-5,
        1.659702063950243e-7,3.728997574317067e-10,
        5.709292619977798e-13,4.932979170744996e-16};
      const double j0pn[5]={9.999999999999999e-1,1.039698629715637,
      	2.576910172633398e-1,1.504152485749669e-2,1.052598413585270e-4};
      const double j0pd[5]={1.0,1.040797262528109,2.588070904043728e-1,
      	1.529954477721284e-2,1.168931211650012e-4};
      const double j0qn[5]={-1.562499999999992e-2,-1.920039317065641e-2,
      	-5.827951791963418e-3,-4.372674978482726e-4,-3.895839560412374e-6};
      const double j0qd[5]={1.0,1.237980436358390,3.838793938147116e-1,
      	3.100323481550864e-2,4.165515825072393e-4};

      const double j1r[7]={7.309637831891357e-5,3.551248884746503e-6,
      	5.820673901730427e-8,4.500650342170622e-10,1.831596352149641e-12,
      	3.891583573305035e-15,3.524978592527982e-18};
      const double j1s[7]={1.0,9.398354768446072e-3,4.328946737100230e-5,
      	1.271526296341915e-7,2.566305357932989e-10,
      	3.477378203574266e-13,2.593535427519985e-16};
      const double j1pn[5]={1.0,1.014039111045313,2.426762348629863e-1,
      	1.350308200342000e-2,9.516522033988099e-5};
      const double j1pd[5]={1.0,1.012208056357845,2.408580305488938e-1,
      	1.309511056184273e-2,7.746422941504713e-5};
      const double j1qn[5]={4.687499999999991e-2,5.652407388406023e-2,
      	1.676531273460512e-2,1.231216817715814e-3,1.178364381441801e-5};
      const double j1qd[5]={1.0,1.210119370463693,3.626494789275638e-1,
      	2.761695824829316e-2,3.240517192670181e-4};

      const double y0r[9]={-7.653778457189104e-3,-5.854760129990403e-2,
      	3.720671300654721e-4,3.313722284628089e-5,4.247761237036536e-8,
      	-4.134562661019613e-9,-3.382190331837473e-11,
      	-1.017764126587862e-13,-1.107646382675456e-16};
      const double y0s[9]={1.0,1.125494540257841e-2,6.427210537081400e-5,
      	2.462520624294959e-7,7.029372432344291e-10,1.560784108184928e-12,
      	2.702374957564761e-15,3.468496737915257e-18,2.716600180811817e-21};
      const double y0pn[5]={9.999999999999999e-1,1.039698629715637,
      	2.576910172633398e-1,1.504152485749669e-2,1.052598413585270e-4};
      const double y0pd[5]={1.0,1.040797262528109,2.588070904043728e-1,
      	1.529954477721284e-2,1.168931211650012e-4};
      const double y0qn[5]={-1.562499999999992e-2,-1.920039317065641e-2,
      	-5.827951791963418e-3,-4.372674978482726e-4,-3.895839560412374e-6};
      const double y0qd[5]={1.0,1.237980436358390,3.838793938147116e-1,
      	3.100323481550864e-2,4.165515825072393e-4};

      const double y1r[8]={-1.041835425863234e-1,-1.135093963908952e-5,
      	2.212118520638132e-4,1.270981874287763e-6,
      	-3.982892100836748e-8,-4.820712110115943e-10,
      	-1.929392690596969e-12,-2.725259514545605e-15};
      const double y1s[8]={1.0,1.186694184425838e-2,7.121205411175519e-5,
      	2.847142454085055e-7,8.364240962784899e-10,1.858128283833724e-12,
      	3.018846060781846e-15,3.015798735815980e-18};
      const double y1pn[5]={1.0,1.014039111045313,2.426762348629863e-1,
      	1.350308200342000e-2,9.516522033988099e-5};
      const double y1pd[5]={1.0,1.012208056357845,2.408580305488938e-1,
      	1.309511056184273e-2,7.746422941504713e-5};
      const double y1qn[5]={4.687499999999991e-2,5.652407388406023e-2,
      	1.676531273460512e-2,1.231216817715814e-3,1.178364381441801e-5};
      const double y1qd[5]={1.0,1.210119370463693,3.626494789275638e-1,
      	2.761695824829316e-2,3.240517192670181e-4};

      double ax;
    	T nump,denp,numq,denq,xx,y,z;

      Bessel() : xj00( 5.783185962946785 ),
                 xj10( 3.047126234366209e1 ),
                 xj01( 1.468197064212389e1 ),
                 xj11( 4.921845632169460e1 ),
                 twoopi( 0.6366197723675813 ),
                 pio4( 0.7853981633974483 )
      {
      }

      void rat(const T x, const double *r, const double *s, const int n)
      {
    		y = x * x;
    		z = 64.0 - y;
    		nump = r[ n ];
    		denp = s[ n ];
    		for ( int i = n-1; i >= 0; i-- )
        {
    			nump = nump * z + r[ i ];
    			denp = denp * y + s[ i ];
    		}
    	}

      void asp( const double *pn, const double *pd, const double *qn,
                const double *qd, const double fac)
      {
    		z = 8.0 / ax;
    		y = z * z;
    		xx = ax - fac * pio4;
    		nump = pn[ 4 ];
    		denp = pd[ 4 ];
    		numq = qn[ 4 ];
    		denq = qd[ 4 ];
    		for ( int i = 3; i >= 0; i-- )
        {
    			nump = nump * y + pn[ i ];
    			denp = denp * y + pd[ i ];
    			numq = numq * y + qn[ i ];
    			denq = denq * y + qd[ i ];
    		}
    	}

      /// The Bessel function \f$ J_0(x) \f$
      /// \param x The arguement of the function
      /// \return The Bessel function of the first kind \f$ J_0(x) \f$
      T J0( const T& x )
      {
    		if ( ( ax = abs( std::real(x) ) ) < 8.0 ) {
    			rat( x, j0r, j0s, 6 );
    			return nump * ( y - xj00 ) * ( y - xj10 ) / denp;
    		} else {
    			asp( j0pn, j0pd, j0qn, j0qd, 1. );
    			return sqrt( twoopi / ax ) * ( cos( xx ) * nump / denp - z * sin( xx )
                 * numq / denq );
    		}
    	}

      /// The Bessel function \f$ J_1(x) \f$
      /// \param x The arguement of the function
      /// \return The Bessel function of the first kind \f$ J_1(x) \f$
      T J1( const T& x ) {
    		if ( ( ax = abs( std::real(x) ) ) < 8.0 ) {
    			rat( x, j1r, j1s, 6 );
    			return x * nump * ( y - xj01 ) * ( y - xj11 ) / denp;
    		} else {
    			asp( j1pn, j1pd, j1qn, j1qd, 3. );
    			T ans = sqrt( twoopi / ax ) * ( cos( xx ) * nump / denp
                       - z * sin( xx ) * numq / denq );
    			return std::real( x ) > 0.0 ? ans : -ans;
    		}
    	}

      /// The Bessel function \f$ Y_0(x) \f$
      /// \param x The arguement of the function
      /// \return The Bessel function of the second kind \f$ Y_0(x) \f$
      T Y0( const T& x ) {
    		if ( std::real( x ) < 8.0) {
    			T j0x = J0( x );
    			rat( x, y0r, y0s, 8 );
    			return nump / denp + twoopi * j0x * log( x );
    		} else {
    			ax = std::real( x );
    			asp( y0pn, y0pd, y0qn, y0qd, 1. );
    			return sqrt( twoopi / x ) * ( sin( xx ) * nump / denp
                 + z * cos( xx ) * numq / denq );
    		}
    	}

      /// The Bessel function \f$ Y_1(x) \f$
      /// \param x The arguement of the function
      /// \return The Bessel function of the second kind \f$ Y_1(x) \f$
      T Y1( const T& x ) {
    		if ( std::real( x ) < 8.0 ) {
    			T j1x = J1( x );
    			rat( x, y1r, y1s, 7 );
    			return x * nump / denp + twoopi * ( j1x * log( x ) - 1.0 / x );
    		} else {
    			ax = std::real( x );
    			asp( y1pn, y1pd, y1qn, y1qd, 3. );
    			return sqrt( twoopi / x ) * ( sin( xx ) * nump / denp
                 + z * cos( xx ) * numq / denq );
    		}
    	}

      /// The Bessel function \f$ J_n(x) \f$
      /// \param n The Bessel function number (integer)
      /// \param x The arguement of the function
      /// \return The Bessel function of the first kind \f$ J_n(x) \f$
      T Jn( const int n, const T& x)
      {
        const double ACC = 160.0;
      	const int IEXP = std::numeric_limits<T>::max_exponent / 2;
      	bool jsum;
      	int j , k, m;
      	T bj , bjm, bjp, dum, sum, tox, ans;
        double ax;

      	if (n==0) return J0(x);
      	if (n==1) return J1(x);
        // ax=abs(x);
      	ax = abs( std::real( x ) );
      	if ( ax * ax <= 8.0 * std::numeric_limits<double>::min() ) return 0.0;
      	else if ( ax > double( n ) ) {
      		//tox=2.0/ax;
      		//bjm=J0(ax);
      		//bj=J1(ax);
          tox = 2.0 / x;
          bjm = J0( x );
          bj = J1( x );
      		for ( j = 1; j < n; j++ ) {
      			bjp = 1. * j * tox * bj - bjm;
      			bjm = bj;
      			bj = bjp;
      		}
      		ans = bj;
      	} else {
      		//tox=2.0/ax;
          tox = 2.0 / x;
      		m = 2 * ( ( n + int( sqrt( ACC * n ) ) ) / 2 );
      		jsum=false;
      		bjp=ans=sum=0.0;
      		bj=1.0;
      		for ( j = m; j > 0; j-- ) {
      			bjm = 1. * j * tox * bj - bjp;
      			bjp = bj;
      			bj = bjm;
      			dum.real( frexp( std::real( bj ), &k ) );
      			if (k > IEXP) {
      				bj.real( ldexp( std::real( bj ), -IEXP ) );
      				bjp.real( ldexp( std::real( bjp ), -IEXP ) );
      				ans.real( ldexp( std::real( ans ), -IEXP ) );
      				sum.real( ldexp( std::real( sum ), -IEXP ) );
      			}
            dum.imag( frexp( std::imag( bj ), &k ) );
      			if (k > IEXP) {
      				bj.imag( ldexp( std::imag( bj ), -IEXP ) );
      				bjp.imag( ldexp( std::imag( bjp ), -IEXP ) );
      				ans.imag( ldexp( std::imag( ans ), -IEXP ) );
      				sum.imag( ldexp( std::imag( sum ), -IEXP ) );
      			}
      			if ( jsum ) sum += bj;
      			jsum =! jsum;
      			if (j == n) ans = bjp;
      		}
      		sum = 2.0 * sum - bj;
      		ans /= sum;
      	}
      	return std::real( x ) < 0.0 && (n & 1) ? -ans : ans;
      }

      /// The Bessel function \f$ Y_n(x) \f$
      /// \param n The Bessel function number (integer)
      /// \param x The arguement of the function
      /// \return The Bessel function of the second kind \f$ Y_n(x) \f$
      T Yn( const int n, const T& x)
      {
        int j;
      	T by, bym, byp, tox;
      	if ( n == 0 ) return Y0( x );
      	if ( n == 1 ) return Y1( x );
      	tox = 2.0 / x;
      	by = Y1( x );
      	bym = Y0( x );
      	for ( j = 1; j < n; j++ ) {
      		byp =1. * j * tox * by - bym;
      		bym = by;
      		by = byp;
      	}
      	return by;
      }

    }; // End struct Bessel

    //TODO could do modified Bessel functions too? 

    //TODO Airy

    //TODO Elliptic integrals

    //TODO Hypergeometric functions

    //TODO Statistical functions

}  // End of namespace Luna

#endif
