/// \file ODE_BVP.h
/// A class for solving a system of \f$ n \f$ first-order ordinary differential
/// equations subject to \f$ n \f$ boundary conditions defined at the ends of
/// the domain \f$ x = x_{left} \f$ or \f$ x_{right} \f$. The system is defined
/// by \f[ M_0 \frac{d \mathbf{f}}{dx} = \mathbf{R} \f] where
/// \f$ \mathbf{f}=\mathbf{f}(x) \f$ is the vector on unknowns,
/// \f$ \mathbf{R} = \mathbf{R}(\mathbf{f}, x) \f$ is the vector of residuals
/// and \f$ M_0 = M_0(\mathbf{f}, x) \f$ is a matrix.
///
/// The system is solved via Newton iteration by splitting the solution into a
/// known part \f$ \mathbf{f}^g \f$ and a correction \f$ \mathbf{f}^c \f$ such
/// that \f$ \mathbf{f} =\mathbf{f}^g + \mathbf{f}^c \f$. The matrix \f$ M_0 \f$
/// and the residual \f$ \mathbf{R} \f$ are exapanded about the known solution
/// such that \f[ \mathbf{R}(\mathbf{f}, x) = \mathbf{R}(\mathbf{f}^g, x) +
/// \frac{\partial \mathbf{R}}{\partial \mathbf{f}} \bigg\vert_{\mathbf{f}^g}
/// \mathbf{f}^c + O((\mathbf{f}^c)^2), \hspace{0.5cm} M_0(\mathbf{f}, x) =
/// M_0(\mathbf{f}^g, x) + \frac{\partial M_0}{\partial \mathbf{f}}
/// \bigg\vert_{\mathbf{f}^g} \mathbf{f}^c + O((\mathbf{f}^c)^2).\f]
/// The system is linearised for the corrections such that \f[ M_0(\mathbf{f}^g,
/// x) \frac{d \mathbf{f}^c}{d x} + \left(\frac{\partial M_0}{\partial
/// \mathbf{f}} \bigg\vert_{\mathbf{f}^g} - \frac{\partial \mathbf{R}}{\partial
/// \mathbf{f}} \bigg\vert_{\mathbf{f}^g} \right) \mathbf{f}^c = \mathbf{R}
/// - M_0(\mathbf{f}^g, x) \frac{d \mathbf{f}^g}{d x}. \f] The system is then
/// solved for the correction \f$ \mathbf{f}^c \f$ which is used to update the
/// known solution \f$ \mathbf{f}^g \f$, this is repeated until convergence is
/// achieved.

#ifndef ODE_BVP_H
#define ODE_BVP_H

#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "Vector.h"
#include "Error.h"
#include "Matrix.h"
#include "BandedMatrix.h"
#include "Equation.h"
#include "Mesh1D.h"
#include "Arclength.h"

namespace Luna
{
	// Doubly templated class
	template <class T, class X = double>

	/// A templated object for boundary-value problems as systems of first-order
	/// ordinary differential equations.
  class ODE_BVP : public Arclength<T>
  {
  private:
    Mesh1D<T, X> SOLUTION;                // Solution mesh
    int MAX_ITER;			                    // Maximum number of iterations
    double TOL;						                // Tolerance for convergence
    Equation_1matrix<T, X> *ptr_EQUATION; // Pointer to the ODE equation
    Residual<T> *ptr_LEFT_RESIDUAL;       // Pointer to the left residual
    Residual<T> *ptr_RIGHT_RESIDUAL;      // Pointer to the right residual
    int LAST_DET_SIGN;                    // Determinant sign of the Jacobian
    bool MONITOR_DET;                     // Monitor the determinant if true

    // Solve the system for an initial guess by Newton iteration. This method is
    // inherited from Arclength and points to solve_bvp
    void solve( Vector<T>& state );

    // Assemble the Jacobian matrix and the residual vector
    void assemble_matrix_problem( BandedMatrix<T>& a, Vector<T>& b );

  public:

    /// Constructor
    ODE_BVP( Equation_1matrix<T, X>* ptr_equation, const Vector<X>& nodes,
             Residual<T>* ptr_left_residual, Residual<T>* ptr_right_residual );

    //TODO copy constructor

    /// Destructor
    virtual ~ODE_BVP();

    /* ----- Methods ----- */

    /// Return a pointer to the current solution mesh
    Mesh1D<T, X>& solution();

    /// Access the tolerance for convergence
    /// \return a pointer to the tolerance
    double& tolerance();

    /// Access the maximum number of iterations
    /// \return a pointer to the maximum number of iterations
    int& max_iterations();

    /// Set whether the determinant will be monitored
    /// \param monitor The boolean value determining whether the determinant
    /// will be monitored
    void set_monitor_det( bool monitor );

    /// Solve the boundary value problem
    void solve_bvp();

    /// Adapt the computational mesh
    /// \param adapt_tol The residual tolerance at a nodal point to determine
    /// whether mesh refinements will occur
    /// \return A pair of values indicating the number of refinements
    /// and unrefinements.
    std::pair<unsigned, unsigned> adapt( const double& adapt_tol );

    /// Solve the system with mesh refinements or unrefinements
    /// \param adapt_tol The residual tolerance at a nodal point to determine
    /// whether mesh refinements will occur
    void adapt_until( const double& adapt_tol );

		/// Initialise so that we can perform arc-length continuation
		void init_arc( T* ptr_param, const double& ds, const double& max_ds );

		/// Arc-length solve the system (init_arc must be called first)
		double arclength_solve( const double& step );

  }; // End of class ODE_BVP

  template <typename T, typename X>
  ODE_BVP<T, X>::ODE_BVP( Equation_1matrix<T, X>* ptr_equation,
           const Vector<X>& nodes, Residual<T>* ptr_left_residual,
           Residual<T>* ptr_right_residual ) :
           MAX_ITER( 20 ), TOL( 1e-8 ), ptr_EQUATION( ptr_equation ),
           ptr_LEFT_RESIDUAL( ptr_left_residual ),
           ptr_RIGHT_RESIDUAL( ptr_right_residual ),
           LAST_DET_SIGN( 0 ), MONITOR_DET( true )
  {
    SOLUTION = Mesh1D<T, X>( nodes, ptr_EQUATION -> get_order() );

    unsigned eqn_order( ptr_EQUATION -> get_order() );
    unsigned left_vars( ptr_LEFT_RESIDUAL -> get_number_of_vars() );
    unsigned right_vars( ptr_RIGHT_RESIDUAL -> get_number_of_vars() );
    unsigned left_order( ptr_LEFT_RESIDUAL -> get_order() );
    unsigned right_order( ptr_RIGHT_RESIDUAL -> get_order() );

    if ( ( left_vars != eqn_order ) || ( right_vars != eqn_order ) ||
         ( left_order + right_order != eqn_order ) )
    {
     std::cout << "order " << eqn_order << std::endl;
     std::cout << "left nvars " << left_vars << std::endl;
     std::cout << "right nvars " << right_vars << std::endl;
     std::cout << "left order " << left_order << std::endl;
     std::cout << "right order " << right_order << std::endl;
     std::string problem;
     problem  = "The ODE_BVP equation and boundary conditions are not well\n";
     problem += "posed. The number of variables for each boundary condition\n";
     problem += "has to be the same as the order of the equation. Also the \n";
     problem += "order of both boundary conditions has to sum to the order\n";
     problem += "of the equation.\n";
     throw Error( problem );
    }
  }

	template <typename T, typename X>
  ODE_BVP<T, X>::~ODE_BVP()
	{}

  /* ---- Methods ----- */

  template <typename T, typename X>
  Mesh1D<T, X>& ODE_BVP<T, X>::solution()
  {
    return SOLUTION;
  }

  template <typename T, typename X>
  double& ODE_BVP<T, X>::tolerance()
  {
    return TOL;
  }

  template <typename T, typename X>
  int& ODE_BVP<T, X>::max_iterations()
  {
    return MAX_ITER;
  }

  template <typename T, typename X>
  void ODE_BVP<T, X>::set_monitor_det( bool monitor )
  {
    MONITOR_DET = monitor;
  }

  template <typename T, typename X>
  void ODE_BVP<T, X>::solve_bvp()
  {
    unsigned order( ptr_EQUATION -> get_order() );
    unsigned n( SOLUTION.nnodes() );
    double max_residual( 1.0 );
    int counter( 0 );
    int det_sign( LAST_DET_SIGN );
    // Banded LHS matrix
    BandedMatrix<T> a( n * order, 2 * order - 1, 2 * order - 1 );
    // RHS vector
    Vector<T> b( n * order, 0.0 );
    // Loop until convergence
    do {
      ++counter;
      assemble_matrix_problem( a, b );
      //TODO actions before linear solve ?
      max_residual = b.norm_inf();
      b = a.solve( b );
      // Put the solution into the 1D mesh
      for ( std::size_t var = 0; var < order; ++var )
      {
        for ( std::size_t i = 0; i < n; ++i )
        {
          SOLUTION( i, var ) += b[ i * order + var ];
        }
      }

    } while( ( max_residual > TOL ) && ( counter < MAX_ITER ) );

    if ( counter >= MAX_ITER )
    {
      std::string problem( "solve_bvp method error: Too many iterations" );
      throw Error( problem );
    }

    if ( MONITOR_DET )
    {
      T det = a.det();
      if ( std::real( det ) < 0 ) {
        det_sign = -1;
      } else {
        det_sign = 1;
      }
      if ( det_sign * LAST_DET_SIGN < 0 )
      {
        LAST_DET_SIGN = det_sign;
        std::string problem;
        problem = "Determinant monitor has changed signs in ODE_BVP.\n";
        problem += "Bifurcation detected.\n";
        throw Error( problem );
      }
      LAST_DET_SIGN = det_sign;
    }
  }

  typedef std::complex<double> cmplx;
  template<>
  std::pair<unsigned, unsigned> ODE_BVP<double, cmplx>::adapt(
                                                       const double& adapt_tol )
  {
    std::string problem;
    problem = " The ODE_BVP.adapt method has been called for a problem in \n";
    problem += " the complex plane. This doesn't make sense as the path is \n";
    problem += " not geometrically defined. \n";
    throw Error( problem );
  }

  template<>
  std::pair<unsigned, unsigned> ODE_BVP<cmplx, cmplx>::adapt(
                                                       const double& adapt_tol )
  {
    std::string problem;
    problem = " The ODE_BVP.adapt method has been called for a problem in \n";
    problem += " the complex plane. This doesn't make sense as the path is \n";
    problem += " not geometrically defined. \n";
    throw Error( problem );
  }

  template <typename T, typename X>
  std::pair<unsigned, unsigned> ODE_BVP<T, X>::adapt( const double& adapt_tol )
  {
    unsigned order( ptr_EQUATION -> get_order() );
    unsigned n( SOLUTION.nnodes() );

    Vector<T> F_node( order, 0.0 );
    Vector<T> R_node( order, 0.0 );

    std::vector<bool> refine( n, false );
    std::vector<bool> unrefine( n, false );

    Vector<double> Res2( n, 0.0 );

    for ( std::size_t node = 1; node <= n - 2; node += 2 )
    {
      // set the current solution at this node
      for ( unsigned var = 0; var < order; ++var )
      {
         F_node[ var ] = SOLUTION( node, var );
      }
      // set the y-pos in the eqn
      ptr_EQUATION -> coord(0) = SOLUTION.coord( node );
      // Update the equation to the nodal position
      ptr_EQUATION -> update( F_node );
      // step size centred at the node
      X invh = 1. / ( SOLUTION.coord( node + 1 ) - SOLUTION.coord( node - 1 ) );
      // loop over all the variables
      Vector<T> temp( order, 0.0 );
      for ( unsigned var = 0; var < order; ++var )
      {
        temp[ var ] = ptr_EQUATION -> residual()[ var ] - (
                SOLUTION( node + 1, var ) - SOLUTION( node - 1, var ) ) * invh;
      }
      Res2[ node ] = temp.norm_inf();
      if ( Res2[ node ] > adapt_tol )
      {
        refine[ node ] = true;
      }
      else if ( Res2[ node ] < TOL / 10. )
      {
        unrefine[ node ] = true;
      }
    }

    std::size_t no_refined( 0 ), no_unrefined( 0 );
    for ( std::size_t i = 0; i < refine.size(); ++i )
    {
      if ( refine[ i ] == true ){ no_refined++; }
      if ( unrefine[ i ] == true ){ no_unrefined++; }
    }

    // make a new refined/unrefined mesh
    Vector<X> nodes( SOLUTION.nodes() );
    Vector<X> newX;
    newX.push_back( nodes[ 0 ] );
    for ( std::size_t i = 1; i < n - 1; ++i )
    {
      if ( refine[ i ] )
      {
        if ( !refine[ i - 1 ] )
        {
          // Have not already refined to the left
          // so refine left AND right with new nodes
          X left( nodes[ i - 1 ] );
          X right( nodes[ i + 1 ] );
          X here( nodes[ i ] );
          newX.push_back( ( left + here ) / 2. );
          newX.push_back( here );
          newX.push_back( ( right + here ) / 2. );
        }
        else
        {
          // already left refined, so just refine right
          X right( nodes[ i + 1 ] );
          X here( nodes[ i ] );
          newX.push_back( here );
          newX.push_back( ( right + here ) / 2. );
        }
      }
      else if ( !unrefine[ i ] )
      {
        newX.push_back( nodes[ i ] );
      }
    }
    newX.push_back( nodes[ n - 1 ] );

    SOLUTION.remesh( newX );

    LAST_DET_SIGN = 0;

    std::pair< unsigned, unsigned > feedback;
    feedback.first = no_refined;
    feedback.second = no_unrefined;
    return feedback;
  }

  template<>
  void ODE_BVP<double, cmplx>::adapt_until( const double& adapt_tol )
  {
    std::string problem;
    problem = " The ODE_BVP.adapt_until method has been called for a problem\n";
    problem += " in the complex plane. This doesn't make sense as the path\n";
    problem += " is not geometrically defined. \n";
    throw Error( problem );
  }

  template<>
  void ODE_BVP<cmplx, cmplx>::adapt_until( const double& adapt_tol )
  {
    std::string problem;
    problem = " The ODE_BVP.adapt_until method has been called for a problem\n";
    problem += " in the complex plane. This doesn't make sense as the path\n";
    problem += " is not geometrically defined. \n";
    throw Error( problem );
  }

  template <typename T, typename X>
  void ODE_BVP<T, X>::adapt_until( const double& adapt_tol )
  {
    std::pair<unsigned, unsigned> changes;
    do {
      changes = adapt( adapt_tol );
      solve_bvp();
      std::cout << " * Adapting: " << changes.first << " " << changes.second
                << std::endl;
    }while ( changes.first + changes.second != 0 );
  }

	template <typename T, typename X>
	void ODE_BVP<T,X>::init_arc( T* ptr_param, const double& ds,
															 const double& max_ds )
  {
    Vector<T> state( SOLUTION.vars_as_vector() );
    this -> Arclength<T>::init_arc( state, ptr_param, ds, max_ds );
  }

	template <typename T, typename X>
	double ODE_BVP<T, X>::arclength_solve( const double& step )
	{
		this -> ds() = step;
		std::size_t order( ptr_EQUATION -> get_order() );  // Order of the equation
		std::size_t n( SOLUTION.nnodes() );            		 // Number of nodes
		std::size_t Size( order * n );

		BandedMatrix<T> Jac( n * order, 2 * order - 1, 2 * order - 1 );                 // Jacobian matrix
		// Residuals over all nodes
		Vector<T> Res1( Size, 0.0 );
		Vector<T> Res2( Size, 0.0 );
		Vector<T> dRes_dp( Size, 0.0 );
		// RHS vectors for the linear solvers
		Vector<T> y( Size, 0.0 );
		Vector<T> z( Size, 0.0 );
		Vector<T> Jac_E;
		// Make backups in case we can't find a converged solution
		Vector<T> backup_state( SOLUTION.vars_as_vector() );
		T backup_parameter( *( this -> ptr_PARAM ) );
		// Generate a 1st order guess for the next state and parameter
		Vector<T> x( this -> LAST_X + this -> X_DERIV_S * this -> DS );
		*( this -> ptr_PARAM ) = this -> LAST_PARAM + this -> PARAM_DERIV_S * this -> DS;

		SOLUTION.set_vars_from_vector( x );            // Update the solution mesh
		bool step_succeeded( false );                  // Check for success
		std::size_t itn( 0 );                          // Iteration counter

		do
		{
			++itn;
			double E1 = this -> arclength_residual( x );            // Arclength residual
			assemble_matrix_problem( Jac, Res1 );                   // Assemble matrix problem
			if ( Res1.norm_inf() < TOL && itn > 1 )
			{
				step_succeeded = true;
				break;
			}
			y = Jac.solve( Res1 );
			// Derivatives wrt parameter
			const double delta( 1.e-8 );
			*( this -> ptr_PARAM ) += delta;
			assemble_matrix_problem( Jac, Res2 );
			double E2 = this -> arclength_residual( x );
			*( this -> ptr_PARAM ) -= delta;
			dRes_dp = ( Res2 - Res1 ) / delta;
			double dE_dp = ( E2 - E1 ) / delta;
			z = Jac.solve( dRes_dp );
			Jac_E = this -> Jac_arclength_residual( x );
			T delta_p = -( E1 + Jac_E.dot( y ) ) / ( dE_dp + Jac_E.dot( z ) );
			Vector<T> delta_x = y + z * delta_p;
			// Update the state variables and the parameter with the corrections
			x += delta_x;
			*( this -> ptr_PARAM ) += delta_p;
			SOLUTION.set_vars_from_vector( x );
			// Check for convergence
			if ( delta_x.norm_inf() < TOL )
			{
				step_succeeded = true;
				break;
			}
			if ( itn > MAX_ITER )
			{
				step_succeeded = false;
				break;
			}
		 }while( true );
		 // If this isn't a successful step
		 if ( !step_succeeded )
		 {
				// Restore things using the backups we made earlier
				SOLUTION.set_vars_from_vector( backup_state );            // Restore state
				*( this -> ptr_PARAM ) = backup_parameter;                // Restore parameter
				this -> DS /= this -> ARCSTEP_MULTIPLIER;                 // Reduce step length
		 }
		 else
		 {
				this -> update( SOLUTION.vars_as_vector() );
				if ( itn > 8 || std::abs( this -> DS ) > this -> MAX_DS )
				{
					// Converging too slowly, so decrease DS
					this -> DS /= this -> ARCSTEP_MULTIPLIER;
				}
				if ( itn < 4 ) // Converging too quickly, so increase DS
				{
					if ( std::abs( this -> DS * this -> ARCSTEP_MULTIPLIER ) < this -> MAX_DS )
					{
						this -> DS *= this -> ARCSTEP_MULTIPLIER;
					}
				}
		 }
		 return this -> DS;
	}

  /* ----- Private methods ----- */

  template <typename T, typename X>
  void ODE_BVP<T, X>::assemble_matrix_problem( BandedMatrix<T>& a,
                                               Vector<T>& b )
  {
    a.fill( 0.0 );
    unsigned order( ptr_EQUATION -> get_order() );
    unsigned n( SOLUTION.nnodes() );
    std::size_t row( 0 );

    Matrix<T> Jac_midpt( order, order, 0.0 );
    Vector<T> F_midpt( order, 0.0 );
    Vector<T> R_midpt( order, 0.0 );
    Vector<T> state_dy( order, 0.0 );
    Matrix<T> h0( order, order, 0.0 );
    // update the BC residuals
    ptr_LEFT_RESIDUAL -> update( SOLUTION.get_nodes_vars( 0 ) );
    // add the LHS BCs to the matrix problem
    for ( std::size_t i = 0; i < ptr_LEFT_RESIDUAL -> get_order(); ++i )
    {
      // loop thru variables at LHS of the domain
      for ( std::size_t var = 0; var < order; ++var )
      {
        a( row, var ) = ptr_LEFT_RESIDUAL -> jacobian()( i, var );
      }
      b[ row ] = - ptr_LEFT_RESIDUAL -> residual()[ i ];
      ++row;
    }
    // inner nodes
    for ( std::size_t node = 0; node <= n - 2; ++node )
    {
      const std::size_t lnode( node );
      const std::size_t rnode( node + 1 );
      // inverse of step size
      const X invh = 1. / ( SOLUTION.coord( rnode ) - SOLUTION.coord( lnode ) );
      for ( std::size_t var = 0; var < order; ++var )
      {
        F_midpt[ var ] = ( SOLUTION( lnode, var ) +
                           SOLUTION( rnode, var ) ) / 2.;
        state_dy[ var ] = ( SOLUTION( rnode, var ) -
                            SOLUTION( lnode, var ) ) * invh;
      }
      X y_midpt = ( SOLUTION.coord( lnode ) + SOLUTION.coord( rnode ) ) / 2.;
      ptr_EQUATION -> coord( 0 ) = y_midpt;
      ptr_EQUATION -> update( F_midpt );
      ptr_EQUATION -> get_jacobian_of_matrix0_mult_vector( F_midpt,
                                                           state_dy, h0 );
      for ( std::size_t var = 0; var < order; ++var )
      {
        for ( std::size_t i = 0; i < order; ++i )
        {
          // Jac of matrix * F_y * g terms
          int col_l( order * lnode + i );
          int col_r( order * rnode + i );
          a( row, col_l ) += h0( var, i )/2.;
          a( row, col_r ) += h0( var, i )/2.;
          // Matrix * g_y terms
          a( row, col_l ) -= ptr_EQUATION -> matrix0()( var, i ) * invh;
          a( row, col_r ) += ptr_EQUATION -> matrix0()( var, i ) * invh;
          // Jacobian of RHS terms
          a( row, col_l ) -= ptr_EQUATION -> jacobian()( var, i ) / 2.;
          a( row, col_r ) -= ptr_EQUATION -> jacobian()( var, i ) / 2.;
        }
        b[ row ] = ptr_EQUATION -> residual()[ var ]
                 - state_dy.dot( ptr_EQUATION -> matrix0()[ var ] );
        row++;
      }
    }
    // update the BC residuals
    ptr_RIGHT_RESIDUAL -> update( SOLUTION.get_nodes_vars( n - 1 ) );
    // add the RHS BCs to the matrix problem
    for ( std::size_t i = 0; i < ptr_RIGHT_RESIDUAL -> get_order(); ++i )
    {
      // loop thru variables at RHS of the domain
      for ( std::size_t var = 0; var < order; ++var )
      {
        int col( order * ( n - 1 ) + var );
        a( row, col ) = ptr_RIGHT_RESIDUAL -> jacobian()( i, var );
      }
      b[ row ] = - ptr_RIGHT_RESIDUAL -> residual()[ i ];
      ++row;
    }

    if ( row != n * order )
     {
       std::string problem;
       problem += "ODE_BVP error: incorrect number of boundary conditions.";
       throw Error( problem );
     }
  }

	template <typename T, typename X>
	void ODE_BVP<T, X>::solve( Vector<T>& state )
	{
		SOLUTION.set_vars_from_vector( state );
		solve_bvp();
		state = SOLUTION.vars_as_vector();
	}



}  // End of namespace Luna

#endif
