/// \file Equation.h
/// A templated class for equations that can be inherited from
/// to allow instantiation of ODE/PDE objects using the resulting class.
/// An equation class is simply a square residual class with one
/// extra coordinate for the ODE_BVP solver.


#ifndef EQUATION_H
#define EQUATION_H

#include "Residual_with_coords.h"
#include "Vector.h"
#include "Matrix.h"

namespace Luna
{

  // An equation object base class used in the BVP/IVP classes.
  // An equation object is a square residual object with an independent variable
  // data member and access methods.
  template < class T, class X = double >
  class Equation : public Residual_with_coords<T, X>
  {
  public:
    /// Constructor for equation class.
    /// \param order The order of the system
    Equation( const unsigned& order );

    /// An empty destructor, virtual since we have virtual methods.
    virtual ~Equation();

  }; // End of class Equation

  template <typename T, typename X>
  Equation<T, X>::Equation( const unsigned& order ) :
      Residual_with_coords<T, X>( order, 1 )
  {
  }

  template <typename T, typename X>
  Equation<T, X>::~Equation()
  {
  }

} // End of namespace Luna

#endif
