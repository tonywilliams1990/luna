/// \file Newton.h
/// A class for using Newton's method for solving systems of nonlinear equations
/// of the form F(x) = 0 where x is a vector and F is vector valued function.

#ifndef NEWTON_H
#define NEWTON_H

#include <iostream>
#include <string>

#include "Matrix.h"
#include "Vector.h"
#include "Error.h"
#include "Timer.h"
#include "Residual.h"
#include "Arclength.h"

namespace Luna
{
	/// A templated Newton iteration class
	template <class T>

	class Newton : public Arclength< T >
	{
    protected:
      Residual<T>* ptr_RESIDUAL;     // Pointer to the residual object
      std::size_t MAX_ITER;          // Maximum number of iterations
      double TOL;                    // Convergence tolerance
      double DELTA;                  // Derivative step
      std::size_t ORDER;             // Order of the system

    public:
      /// Constructor
			/// \param ptr_residual Pointer to the residual object
			/// \param max_iter Maximum number of iterations
			/// \param tolerance Convergence tolerance
			/// \param delta_step Finite-difference derivative step
      Newton( Residual<T>* ptr_residual, std::size_t max_iter = 20,
              double tolerance = 1.0e-8, double delta_step = 1.0e-8 ) :
      ptr_RESIDUAL( ptr_residual ),
      MAX_ITER( max_iter ),
      TOL( tolerance )
      {
        ptr_RESIDUAL -> delta() = delta_step;
        DELTA = delta_step;
        ORDER = ptr_RESIDUAL -> get_order();
      }

      /// Destructor
      ~Newton() {}

      /* ----- Methods ----- */

      /// Newton iterate using intial guess x_0 (stores the solution in x_0)
			/// \param x_0 The initial guess
      void iterate( Vector<T>& x_0 );

      /// Solve the system using Newton iteration (points to iterate)
			/// \param x_0 The initial guess
      void solve( Vector<T>& x_0 )
      {
        iterate( x_0 );
      }

      /// Arc-length solve the system (init_arc must be called first)
			/// \param x The current state Vector
      void arclength_solve( Vector<T>& x );

	}; // End of class Newton

	template <typename T>
	inline void Newton<T>::iterate( Vector<T>& x_0 )
	{
		//TODO determinant monitoring etc
		if ( x_0.size() != ORDER )
		{
			throw Error( "Newton error: size does not agree." );
		}
		Vector<T> F( ORDER, 0.0 );               // Residual function evaluation
		Matrix<T> J( ORDER, ORDER, 0.0 );        // Jacobian matrix
		Vector<T> dx( ORDER, 0.0 );              // Increment vector
		std::size_t iter( 0 );
		double max_residual;

		do
		{
			++iter;
			ptr_RESIDUAL -> update( x_0 );
			F = ptr_RESIDUAL -> residual();
			max_residual = F.norm_inf();
			J = ptr_RESIDUAL -> jacobian();

			dx = J.solve( -F );
			x_0 += dx;
		}while( iter < MAX_ITER && max_residual > TOL );

		if ( iter == MAX_ITER )
		{
			std::string warning;
			warning = "NEWTON WARNING! MAXIMUM NUMBER OF ITERATIONS REACHED.";
			std::cout << warning << std::endl;
		}
	}

	template <typename T>
	inline void Newton<T>::arclength_solve( Vector<T>& x )
	{
		if ( !this->INITIALISED ) // Error if init_arc has not been called
		{
			throw Error( "Newton arclength_solve error: not initialised." );
		}
		Vector<T> backup_state( x );
		T backup_parameter( *(this->ptr_PARAM) );
		bool step_succeeded( false );
		std::size_t itn( 0 );
		// Guess the next solution
		*( this->ptr_PARAM ) = this->LAST_PARAM + this->PARAM_DERIV_S * this->DS;
		x = this->LAST_X + this->X_DERIV_S * this->DS;

		Matrix<T> J( ORDER, ORDER, 0.0 );
		Vector<T> R1( ORDER, 0.0 );
		Vector<T> R2( ORDER, 0.0 );
		Vector<T> dR_dp( ORDER, 0.0 );
		Vector<T> y( ORDER, 0.0 );
		Vector<T> z( ORDER, 0.0 );
		Vector<T> JacE; 

		do
		{
			this->ptr_RESIDUAL -> update( x );              // Update the residual
			J = this->ptr_RESIDUAL -> jacobian();           // Get the Jacobian
			R1 = this->ptr_RESIDUAL -> residual();          // Get the residual
			double E1 = this -> arclength_residual( x );    // Arclength residual
			// Compute derivatives wrt the parameter
			*(this->ptr_PARAM) += this->DELTA;
			this->ptr_RESIDUAL -> residual_fn( x, R2 );
			double E2 = this -> arclength_residual( x );
			*(this->ptr_PARAM) -= this->DELTA;
			dR_dp = ( R2 - R1 ) / this->DELTA;
			T dE_dp = ( E2 - E1 ) / this->DELTA;
			// Bordering algorithm
			y = J.solve( -R1 );
			z = J.solve( -dR_dp );
			JacE = this->Jac_arclength_residual( x );
			T delta_p = -( E1 + JacE.dot( y ) ) / ( dE_dp + JacE.dot( z ) );
			Vector<T> delta_x = y + z * delta_p;
			double max_correction = std::max( delta_x.norm_inf(),
																				std::abs( delta_p ) );
			if ( max_correction < this->TOL )
			{
				step_succeeded = true;
				break;
			}
			// Add the corrections to the state variables
			x += delta_x;
			*(this->ptr_PARAM) += delta_p;
			++itn;
			if ( itn > MAX_ITER )
			{
				step_succeeded = false;
				break;
			}
		}while( true );

		if ( !step_succeeded )
		{
			x = backup_state;                               // Restore state
			*(this->ptr_PARAM) = backup_parameter;          // Restore the parameter
			this->ptr_RESIDUAL->update( x );                // Restore the residual
			this->DS /= this->ARCSTEP_MULTIPLIER;           // Reduce our step length

		}
		else
		{
			this->update( x );
			if ( itn >= 7 )                                 // Converging to slowly
			{
				this->DS /= this->ARCSTEP_MULTIPLIER;
			}
			if ( itn <= 2 )                                 // Converging to quickly
			{
				if ( std::abs( this->DS * this->ARCSTEP_MULTIPLIER ) < this->MAX_DS )
				{
					this->DS *= this->ARCSTEP_MULTIPLIER;
				}
			}
		}
	}

}  // End of namespace Luna

#endif
