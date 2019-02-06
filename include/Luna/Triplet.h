/// \file Triplet.h
/// A class for defining Triplets of the form (i,j,value) where i and j are
/// indices. Triplets are useful when filling a SparseMatrix.

#ifndef TRIPLET_H
#define TRIPLET_H

#include <iostream>
#include <string>
#include <cmath>

#include "Vector.h"
#include "Error.h"

namespace Luna
{
  /// A templated class to define triplets
	template <class T>

	class Triplet
	{
		private:
			std::size_t ROW;									// Row index
			std::size_t COL;									// Column index
			T VAL;														// Value

    public:

			/// Constructor for an unspecified Triplet
			Triplet() : ROW( 0 ), COL( 0 ), VAL( 0 ) {}

			/// Constructor for a Triplet with specified values
			/// \param row Row index
			/// \param col Column index
			/// \param val Value
			Triplet( const std::size_t& row, const std::size_t& col, const T& val )
				: ROW( row ), COL( col ), VAL( val ) {}

			/// Copy constructor
		  /// \param source The source Triplet to be copied
		  Triplet( const Triplet<T>& source );

			/// Destructor
			~Triplet(){}

			/* ----- Operator overloading ----- */

			//TODO output operator

			/// Access operator (write)
	    /// \param row Row index
	    /// \param col Column index
			/// \param val Value
	    void operator() ( const std::size_t& row, const std::size_t& col,
			 									const T& val )
			{
				ROW = row;
				COL = col;
				VAL = val;
			}

			//TODO equals operator



			/* ----- Methods ----- */

			/// The row index of the Triplet
			/// \return A reference to the row index
			std::size_t& row();

			/// The column index of the Triplet
			/// \return A reference to the column index
			std::size_t& col();

			/// The value stored in the Triplet
			/// \return A reference to the stored value
			T& val();


  }; // End of class Triplet

	template <typename T>
  inline Triplet<T>::Triplet( const Triplet<T>& source )
  {
    *this = source;
  }

	/* ----- Operator overloading ----- */


	/* ----- Methods ----- */

	template <typename T>
	inline std::size_t& Triplet<T>::row()
	{
		return ROW;
	}

	template <typename T>
	inline std::size_t& Triplet<T>::col()
	{
		return COL;
	}

	template <typename T>
	inline T& Triplet<T>::val()
	{
		return VAL;
	}

}  // End of namespace Luna

#endif
