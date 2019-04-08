/// \file ODE_BVP.h
/// A templated object for boundary-value problems as systems of first-order
/// ordinary differential equations.
/// TODO more of description (mention double template), method of solution

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

  /// A Vector class for use with double and std::complex<double>
  class ODE_BVP //TODO : public Arclength<T>
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
    unsigned N( SOLUTION.nnodes() );

    Vector<T> F_node( order, 0.0 );
    Vector<T> R_node( order, 0.0 );

    std::vector<bool> refine( N, false );
    std::vector<bool> unrefine( N, false );

    Vector<double> Res2( N, 0.0 );

    for ( std::size_t node = 1; node <= N - 2; node += 2 )
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
    for ( std::size_t i = 1; i < N - 1; ++i )
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
    newX.push_back( nodes[ N - 1 ] );

    SOLUTION.remesh( newX );

    LAST_DET_SIGN = 0;

    std::pair< unsigned, unsigned > feedback;
    feedback.first = no_refined;
    feedback.second = no_unrefined;
    return feedback;
  }

  template <typename T, typename X>
  void ODE_BVP<T, X>::adapt_until( const double& adapt_tol )
  {

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



}  // End of namespace Luna

#endif
