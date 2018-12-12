/// \file Special.h
/// A file specifying some special mathematical functions not contained
/// in the standard library

#ifndef SPECIAL_H
#define SPECIAL_H

#include <vector>
#include <fstream>
#include <iostream>
#include <complex>

#include "Error.h"
#include "Vector.h"


namespace Luna
{
    //TODO Gamma (double in std but want complex too) -> better to do log-gamma
    // gamma.h + p256(280) in Numerical recipes
    // do proper documentation

    /// The Gamma function \f[ \Gamma(z) = \int_0^{\infty} x^{z-1}e^{-x} dx \f]
    /// \param z The argument of the function (not a negative integer or zero)
    /// \return The Gamma function \f$\Gamma(z)\f$ using Spouge's approximation
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

    /// The Log-Gamma function
    /// \param z The argument of the function (z > 0)
    /// \return The log of the Gamma function \f$\ln \Gamma(z) \f$
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


    //TODO factorial - for ints / std::size_t - do bounds checking ??

    //TODO Beta

    //TODO Error

    //TODO Airy

    //TODO Bessel

    //TODO Elliptic integrals

    //TODO Hypergeometric functions

    //TODO Statistical functions

}  // End of namespace Luna

#endif
