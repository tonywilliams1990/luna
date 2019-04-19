/// \file Basis.h
/// A base class to be inherited by classes specifying basis functions to be
/// used in spectral methods.


#ifndef BASIS_H
#define BASIS_H

#include <iostream>
#include <string>

#include "Matrix.h"
#include "Vector.h"
#include "Error.h"

namespace Luna
{
  /// A templated base class to be inherited by spectral basis objects
	template <class T>

	class Basis
	{
		private:
			std::string IDENTIFIER;		 // String specifying the type of basis function

    public:

      /// Constructor
			Basis();

      /// Destructor ( virtual since we have virtual functions )
      virtual ~Basis();

			/* ----- Operator overloading ----- */

			/// Evaluation operator for the nth basis function at point x
			/// \param x The collocation point where the basis function is evaluated
			/// \param n The degree of the basis function
			/// \return The value of the nth degree basis function at the point x
			virtual T operator()( const T& x, const std::size_t& n );

			/// Evaluation operator for the nth basis function at a Vector of points
			/// \param x A Vector of points where the basis function is evaluated
			/// \param n The degree of the basis function
			/// \return A Vector of values of the nth degree basis function
			virtual Vector<T> operator()( const Vector<T>& x, const std::size_t& n );

			/// Evaluation operator for the dth derivative of the nth basis
			/// function at point x
			/// \param x The collocation point where the basis function is evaluated
			/// \param n The degree of the basis function
			/// \param d The derivative to return
			/// \return The dth derivative of the nth degree basis function at the
			/// point x
			virtual T operator()( const T& x, const std::size_t& n,
														const std::size_t& d );

			/// Evaluation operator for the dth derivative of the nth basis
			/// function at a Vector of points
			/// \param x A Vector of points where the basis function is evaluated
			/// \param n The degree of the basis function
			/// \param d The derivative to return
			/// \return A Vector of the dth derivative of the nth degree basis
			/// function at the Vector of points
			virtual Vector<T> operator()( const Vector<T>& x, const std::size_t& n,
																		const std::size_t& d );

			/// Copy assignment
	    /// \param original Basis function to be copied
	    /// \return The new Basis function
	    Basis<T>& operator=( const Basis<T>& original );

			/* ----- Methods ----- */

			/// Set the basis function identification string
			/// \param identifier The identifier string specifying a type of basis
			void set_identifier( std::string identifier )
			{
				IDENTIFIER = identifier;
			}

			/// Get the basis function identification string
			std::string get_identifier()
			{
				return IDENTIFIER;
			}

  }; // End of class Basis

	template <typename T>
	Basis<T>::Basis()
	{
		IDENTIFIER = "Unspecified";
	}

	template <typename T>
	Basis<T>::~Basis()
	{}

	template <typename T>
	T Basis<T>::operator()( const T& x, const std::size_t& n )
	{
		std::string problem;
		problem += " The Basis function operator () has not been overloaded in \n";
		problem += " the derived class. \n";
		throw Error( problem );
	}

	template <typename T>
	Vector<T> Basis<T>::operator()( const Vector<T>& x, const std::size_t& n )
	{
		std::string problem;
		problem += " The Basis function Vector operator () has not been\n";
		problem += " overloaded in the derived class. \n";
		throw Error( problem );
	}

	template <typename T>
	T Basis<T>::operator()( const T& x, const std::size_t& n,
																	const std::size_t& d )
	{
		std::string problem;
		problem += " The Basis function derivative operator () has not been \n";
		problem += " overloaded in the derived class. \n";
		throw Error( problem );
	}

	template <typename T>
	Vector<T> Basis<T>::operator()( const Vector<T>& x, const std::size_t& n,
																	const std::size_t& d )
	{
		std::string problem;
		problem += " The Basis function Vector derivative operator () has not \n";
		problem += " been overloaded in the derived class. \n";
		throw Error( problem );
	}

	template <typename T>
	Basis<T>& Basis<T>::operator=( const Basis<T>& original )
	{
		if ( this == &original ){ return *this; }
    IDENTIFIER = original.IDENTIFIER;
    return *this;
	}

}  // End of namespace Luna

#endif
