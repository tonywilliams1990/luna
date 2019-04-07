/// \file Equation_1matrix.h
/// A templated class for equations that can be inherited from
/// to allow instantiation of ODE/PDE objects using the resulting class.
/// An equation class is simply a square residual class with one
/// extra coordinate for the ODE_BVP solver.

#ifndef EQUATION_1MATRIX_H
#define EQUATION_1MATRIX_H

#include "Residual_with_coords.h"
#include "Vector.h"
#include "Matrix.h"

namespace Luna
{

  // An equation object base class used in the BVP/IVP classes.
  // An equation object is a square residual object with an independent variable
  // data member and access methods.
  template < typename T, typename X = double >
  class Equation_1matrix : public Residual_with_coords<T, X>
  {
  public:
   /// Constructor
   /// \param order The order of the system
   explicit Equation_1matrix( const unsigned &order );

   /// Destructor
   virtual ~Equation_1matrix();

   /// Update the Equation object for the current set of state variables
   /// \param state The state vector at which to set the equation object
   void update( const Vector<T> &state );

   /// Return a handle to the matrix
   const Matrix<T>& matrix0() const;

   /// Return the product of the Jacobian-of-the-matrix and a vector 'vec'
   /// when the equation has a given 'state'. The user should overload this
   /// if concerned about performance of the solver. If not overloaded, the
   /// default is to finite difference the Jacobian-of-the-matrix.
   /// \param state The current state variables
   /// \param vec The vector which is multiplied by the Jacobian of the matrix
   /// \param h The resulting 2D matrix
   virtual void get_jacobian_of_matrix0_mult_vector( const Vector<T> &state,
      const Vector<T> &vec, Matrix<T> &h ) const;


 protected:

   /// Define the matrix in terms of the current state vector.
   /// \param x The current state vector.
   /// \param m The matrix.
   virtual void matrix0( const Vector<T> &x, Matrix<T> &m ) const
   {
     std::string problem;
     problem = "The equation::matrix0 method has not been implemented!\n";
     problem += "You have to implement this method to define the equation.\n";
     throw Error( problem );
   }

 private:
   // Matrix0 evaluated for the last state vector
   Matrix<T> MATRIX0_AT_LAST_STATE;

 }; // End of class Equation

 template <typename T, typename X>
 Equation_1matrix<T, X>::Equation_1matrix( const unsigned& order ) :
       Residual_with_coords<T, X>( order, 1 )
 {
   unsigned sys_order( this -> ORDER_OF_SYSTEM );
   MATRIX0_AT_LAST_STATE = Matrix<T>( sys_order, sys_order, 0.0 );
 }

 template <typename T, typename X>
 Equation_1matrix<T, X>::~Equation_1matrix()
 {}

 template <typename T, typename X>
 void Equation_1matrix<T, X>::update( const Vector<T> &state )
 {
     Residual_with_coords<T, X>::update( state );
     matrix0( state, MATRIX0_AT_LAST_STATE );
 }

 template <typename T, typename X>
 inline const Matrix<T>& Equation_1matrix<T, X>::matrix0() const
 {
   return MATRIX0_AT_LAST_STATE;
 }

 template <typename T, typename X>
 void Equation_1matrix<T, X>::get_jacobian_of_matrix0_mult_vector(
  const Vector<T> &state, const Vector<T> &vec, Matrix<T> &h ) const
 {
   // We dont need state in the default implementation as its already been set
   // by the update method. You do need it for the user to overload this method
   // with an explicit analytical version however.

   Vector<T> copy_of_state( this -> LAST_STATE );
   Matrix<T> copy_of_matrix( MATRIX0_AT_LAST_STATE );
   std::vector< Matrix<T> > jacmatrix;
   // update the Jacobian of the mass matrix
   for ( std::size_t i = 0; i < this -> ORDER_OF_SYSTEM; ++i )
   {
     copy_of_state[ i ] += this -> DELTA;
     matrix0( copy_of_state, copy_of_matrix );
     copy_of_state[ i ] -= this -> DELTA;
     copy_of_matrix -= MATRIX0_AT_LAST_STATE;
     copy_of_matrix *= ( 1. / this -> DELTA );
     // the 3D object that represents the Jacobian of the mass matrix
     jacmatrix.push_back( copy_of_matrix );
   }
   // evaluate the jacabian of mass contribution
   for ( unsigned i = 0; i < this -> ORDER_OF_SYSTEM; ++i )
   {
     for ( unsigned j = 0; j < this -> ORDER_OF_SYSTEM; ++j )
     {
       h( i, j ) = jacmatrix[ j ][ i ].dot( vec );
     }
   }
 }

} // End of namespace Luna

#endif
