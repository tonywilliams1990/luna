/// \file Residual.h
/// A residual class for use with Newton for solving systems of
///	non-linear equations of the form F(x) = 0 where x is a vector
///	and F is vector valued function. Given a current approximation x_k
///	of x we can calculate the residual and approximate the Jacobian.

#ifndef RESIDUAL_H
#define RESIDUAL_H

#include <iostream>
#include <string>

#include "Matrix.h"
#include "Vector.h"
#include "Error.h"

namespace Luna
{
  /// A templated base class to be inherited by objects that define residuals
	template <class T>

	class Residual
	{
    protected:
      /// Because the residual evaluation at the current state is assumed
      /// to have already been done by the 'update' method, this routine is
      /// protected. This default uses a finite-differenced Jacobian.
      /// You can overload this to provide an analytic Jacobian if you wish
      virtual void jacobian( const Vector<T>& state, Matrix<T>& jac ) const;

      Matrix<T> JAC_AT_LAST_STATE; // Jacobian for the last state vector
      Vector<T> FN_AT_LAST_STATE;  // Residual for the last state vector
      Vector<T> LAST_STATE;				 // The last state vector
      T DELTA;										 // Step for FD computation of the Jacobian
      unsigned ORDER_OF_SYSTEM; // The order of the system of equations
      unsigned NUMBER_OF_VARS;	 // The number of elements in the state vector

    public:
      /// Constructor for a 'square' residual (N residuals for N unknowns)
      /// \param order The order of the system of equations
      Residual( const unsigned& order ) : DELTA( 1.e-8 )
      {
        ORDER_OF_SYSTEM = order;
        NUMBER_OF_VARS = order;
        LAST_STATE = Vector<T>( NUMBER_OF_VARS, 0.0 );
        FN_AT_LAST_STATE = Vector<T>( ORDER_OF_SYSTEM, 0.0 );
        JAC_AT_LAST_STATE = Matrix<T>( ORDER_OF_SYSTEM, NUMBER_OF_VARS, 0.0 );
      }

      /// Constructor for a 'non-square' residual object
      /// \param order The order of the system of equations
      /// \param nvars The number of elements in the state vector
      Residual( const unsigned& order,
                const unsigned& nvars ) : DELTA( 1.e-8 )
      {
        ORDER_OF_SYSTEM = order;
        NUMBER_OF_VARS = nvars;
        LAST_STATE = Vector<T>( NUMBER_OF_VARS, 0.0 );
        FN_AT_LAST_STATE = Vector<T>( ORDER_OF_SYSTEM, 0.0 );
        JAC_AT_LAST_STATE = Matrix<T>( ORDER_OF_SYSTEM, NUMBER_OF_VARS, 0.0 );
      }

      /// Destructor ( virtual since we have virtual functions )
      virtual ~Residual()
      {}


      /* ----- Methods ----- */

      /// Get the order of the residual vector
      /// \return The order of the system of equations
      unsigned get_order() const { return ORDER_OF_SYSTEM; }

      /// Get the number of variables
      /// \return The number of elements in the state vector
      unsigned get_number_of_vars() const { return NUMBER_OF_VARS; }

      /// Update the Residual object for the current set of state variables
      /// \param The last state vector
      void update( const Vector<T>& state )
      {
        LAST_STATE = state;
        residual_fn( LAST_STATE, FN_AT_LAST_STATE );
        jacobian( LAST_STATE, JAC_AT_LAST_STATE );
      }

      /// A handle to the residuals corresponding to the last state
      /// \return Residual for the last state vector
      const Vector<T>& residual() const { return FN_AT_LAST_STATE; }

      /// A handle to the Jacobian of the residual
      /// \return Jacobian for the last state vector
      const Matrix<T>& jacobian() const { return JAC_AT_LAST_STATE; }

      /// A handle to the state vector
      /// \return The last state vector
      const Vector<T>& state() const { return LAST_STATE; }

      /// A handle to the step size used when finite-differencing
      /// \return The step for FD computation of the Jacobian
      T& delta() { return DELTA; }

      /// A handle to the step size used when finite-differencing
      /// \return The step for FD computation of the Jacobian
      const T& delta() const { return DELTA; }

      /// The residual function to be defined later
      /// \param state The current state Vector
      /// \param f The function Vector to be updated
      virtual void residual_fn( const Vector<T>& state, Vector<T>& f ) const
      {
        std::string problem;
        problem = "The residual_fn method has not been implemented.\n";
        problem += "You must implement this method to define the residual.\n";
        throw Error( problem );
      }

  }; // End of class Residual

  template <typename T>
  inline void Residual<T>::jacobian( const Vector<T>& state,
                                     Matrix<T>& jac ) const
  {
    Vector<T> new_state( state );
    // evaluation of the function
    Vector<T> f_at_new_state( ORDER_OF_SYSTEM, 0.0 );
    // default is to FD the Jacobian
    for ( std::size_t i = 0; i < NUMBER_OF_VARS; ++i )
    {
      new_state[ i ] += DELTA;
      residual_fn( new_state, f_at_new_state );
      new_state[ i ] -= DELTA;
      jac.set_col( i, ( f_at_new_state - FN_AT_LAST_STATE ) / DELTA );
    }
  }
}  // End of namespace Luna

#endif
