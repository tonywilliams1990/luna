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
//TODO remove this #include "SparseMatrix.h"
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
    Mesh1D<T, X> SOLUTION;            // Solution mesh
    int MAX_ITER;			                // Maximum number of iterations
    double TOL;						            // Tolerance for convergence
    Equation<T, X> *ptr_EQUATION; 	  // Pointer to the ODE equation
    Residual<T> *ptr_LEFT_RESIDUAL;   // Pointer to the left residual
    Residual<T> *ptr_RIGHT_RESIDUAL;  // Pointer to the right residual
    int LAST_DET_SIGN;                // Sign of the determinant of the Jacobian
    bool MONITOR_DET;                 // Monitor the determinant if true

    // Solve the system for an initial guess by Newton iteration. This method is
    //  inherited from Arclength and points to solve_bvp
    void solve( Vector<T>& state );

    // Assemble the Jacobian matrix and the residual vector
    void assemble_matrix_problem( BandedMatrix<T>& A, Vector<T>& B );

  public:

    /// Constructor
    ODE_BVP( Equation<T, X>* ptr_equation, const Vector<X>& nodes,
             Residual<T>* ptr_left_residual, Residual<T>* ptr_right_residual );

    //TODO copy constructor

    /// Destructor
    virtual ~ODE_BVP();

    /* ----- Methods ----- */

    /// Return a pointer to the current solution mesh
    Mesh1D<T, X>& solution();

    //TODO tolerance, max_iter, solution etc

    /// Solve the BVP
    void solve_bvp();

  }; // End of class ODE_BVP

  template <typename T, typename X>
  ODE_BVP<T, X>::ODE_BVP( Equation<T, X>* ptr_equation, const Vector<X>& nodes,
           Residual<T>* ptr_left_residual, Residual<T>* ptr_right_residual ) :
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
    

  }

  /* ----- Private methods ----- */

}  // End of namespace Luna

#endif
