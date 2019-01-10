/// \file Arclength.h
/// A base class for arc length solvers.

#ifndef ARCLENGTH_H
#define ARCLENGTH_H

#include <iostream>

#include "Vector.h"
#include "Error.h"

namespace Luna
{
	/// A templated arc length solver class
	template <class T>

	class Arclength
	{
    protected:
      T *ptr_PARAM;              // Pointer to the arc-length parameter
      Vector<T> LAST_X;          // State variable at the last computed solution
      Vector<T> X_DERIV_S;       // Deriv of the state variable wrt arc-length
      T LAST_PARAM;              // Value at the last computed solution
      T PARAM_DERIV_S;           // Deriv of the parameter wrt arc-length
      double DS;                 // Size of the arc-length step
      double MAX_DS;             // Maximum arc length step to be taken
      double ARCSTEP_MULTIPLIER; // Step change multiplier
      bool INITIALISED;          // True if the arc-length solver is initialised
      double THETA;              // Theta weighting paramter (usually = 1/2)

      /* ----- Methods ----- */

      /// Store the current converged state and parameter & compute derivatives
			/// \param x The current state Vector
      void update( const Vector<T>& x );

      /// Extra constraint that is to be used to replace the unknown arc-length
			/// \param x The current state Vector
      double arclength_residual( const Vector<T>& x ) const;

      /// Return the derivative of the arclength_residual function
			/// with respect to each state variable
			/// \param x The current state Vector
			/// \return The derivative of the arclength_residual function
      Vector<T> Jac_arclength_residual( Vector<T>& x ) const;

      /// Automatically update theta value only if RESCALE_THETA = true
			/// \param x The current state Vector
      void update_theta( const Vector<T>& x );

    private:
      double DESIRED_ARC_PROPORTION;   // Desired proportion of the arc-length
      bool RESCALE_THETA;              // If true theta is rescaled at each step

    public:

      /// Constructor
      Arclength() : ARCSTEP_MULTIPLIER( 2.0 ),
                    INITIALISED( false ),
                    THETA( 0.5 ),
                    DESIRED_ARC_PROPORTION( 0.5 ),
                    RESCALE_THETA( false )
      {}

      /// Destructor
      virtual ~Arclength() {}

      /* ----- Methods ----- */

      /// Initialise the arc-length continuation class
			/// \param x The current state Vector
			/// \param ptr_param Pointer to the arc-length parameter
			/// \param ds The size of the arc-length step
			/// \param max_ds The maximum arc-length step to be taken
      void init_arc( Vector<T> x, T* ptr_param, const double& ds,
										 const double& max_ds );

      /// Solve the system for a given initial guess
			/// \param x The current state Vector
      virtual void solve( Vector<T>& x ) = 0;

      /// A handle to the arc-length step
		  /// \return The arc-length step
      double& ds() { return DS; }

      /// A handle to the arc step multiplier
			/// \return The arc-length step multiplier
      double& arcstep_multiplier() { return ARCSTEP_MULTIPLIER; }

      /// A handle to the RESCALE_THETA flag
			/// \return The theta rescale flag
      bool& rescale_theta() { return RESCALE_THETA; }

      /// A handle to the theta parameter
			/// \return The Theta weighting parameter
      double& theta() { return THETA; }

      /// A handle to the desired proportion of the parameter to be used
			/// \return The desired proportion of the arc-length
      double& desired_arc_proportion() { return DESIRED_ARC_PROPORTION; }

	}; // End of class Arclength

  template <typename T>
  void Arclength<T>::update( const Vector<T>& x )
  {
    if ( RESCALE_THETA ) { update_theta( x ); }
    X_DERIV_S = ( x - LAST_X ) / DS;
    PARAM_DERIV_S = ( *ptr_PARAM - LAST_PARAM ) / DS;
    LAST_X = x;
    LAST_PARAM = *ptr_PARAM;
  }

  template <typename T>
  double Arclength<T>::arclength_residual( const Vector<T>& x ) const
  {
    return THETA * ( x - LAST_X ).norm_2() / x.size() + ( 1.0 - THETA )
		       * std::pow( std::abs( *ptr_PARAM - LAST_PARAM ), 2 ) - DS * DS;
  }

	template <typename T>
	Vector<T> Arclength<T>::Jac_arclength_residual( Vector<T>& x ) const
	{
		Vector<T> Jx( x - LAST_X );
		Jx = Jx * ( THETA / ( x.size() * ( x - LAST_X ).norm_2() ) );
		return Jx;
	}

  template <class T>
  void Arclength<T>::update_theta( const Vector<T>& x )
  {
    if ( RESCALE_THETA )
    {
      double Delta_p2 = std::pow( std::abs( *ptr_PARAM - LAST_PARAM ), 2);
      double Delta_x2 = ( x - LAST_X ).norm_2() / x.size();
      THETA = Delta_p2 * ( DESIRED_ARC_PROPORTION - 1.0 ) /
             ( Delta_p2 * ( DESIRED_ARC_PROPORTION - 1.0 )
              - DESIRED_ARC_PROPORTION * Delta_x2 );
    }
  }

  template <class T>
  void Arclength<T>::init_arc( Vector<T> x, T* ptr_param,
                               const double& ds, const double& max_ds )
  {
    ptr_PARAM = ptr_param;
    solve( x );
    LAST_X = x;
    LAST_PARAM = *ptr_PARAM;
    // We now have one solution which we can use to arc-length continue
    DS = ds;
    *ptr_PARAM += DS;
    // Recompute the state at this new parameter value
    solve( x );
    update( x );
    INITIALISED = true;
    MAX_DS = max_ds;
  }

  template <class T>
  void Arclength<T>::solve( Vector<T>& x ) {}

}  // End of namespace Luna

#endif
