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

			/// Output operator
			template <class Type>
	    friend std::ostream& operator<<( std::ostream& os,
	                                     const Triplet<Type>& trip );

			/// Access operator (write)
	    /// \param row Row index
	    /// \param col Column index
			/// \param val Value
	    void operator() ( const std::size_t& row, const std::size_t& col,
			 									const T& val );

			/// Copy assignment
	    /// \param original Triplet to be copied
	    /// \return The new Triplet
	    Triplet<T>& operator=( const Triplet<T>& original );

			/// Compare indices for column major sorting
			/// \param t_1 The first Triplet
			/// \param t_2 The second Triplet
			template <class Type>
			friend bool operator<( const Triplet<Type>& t_1,
														 const Triplet<Type>& t_2 );

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

			/// The row index of the Triplet (read only)
			/// \return The row index
			std::size_t get_row() const;

			/// The column index of the Triplet (read only)
			/// \return The column index
			std::size_t get_col() const;

			/// The value stored in the Triplet (read only)
			/// \return The stored value
			T get_val() const ;

  }; // End of class Triplet

	template <typename T>
  inline Triplet<T>::Triplet( const Triplet<T>& source )
  {
    *this = source;
  }

	/* ----- Operator overloading ----- */

	template <class Type>
  inline std::ostream& operator<<( std::ostream& os, const Triplet<Type>& trip )
  {
		os << "( " << trip.ROW << ", " << trip.COL << ", " << trip.VAL << " )";
    return os;
  }

	template <typename T>
	inline void Triplet<T>::operator() ( const std::size_t& row,
																			 const std::size_t& col, const T& val )
	{
		ROW = row;
		COL = col;
		VAL = val;
	}

	template <typename T>
  inline Triplet<T>& Triplet<T>::operator=( const Triplet<T>& original )
  {
    if ( this == &original ){ return *this; }
    ROW = original.ROW;
		COL = original.COL;
		VAL = original.VAL;
    return *this;
  }

	template <class Type>
	inline bool operator<( const Triplet<Type>& t_1, const Triplet<Type>& t_2 )
	{
		bool compare;
		if ( t_1.COL < t_2.COL )
		{
			compare = true;
		} else if ( t_1.COL == t_2.COL ) {
			if ( t_1.ROW < t_2.ROW )
			{
				compare = true;
			} else {
				compare = false;
			}
		} else {
			compare = false;
		}
		return compare;
	}


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

	template <typename T>
	inline std::size_t Triplet<T>::get_row() const
	{
		return ROW;
	}

	template <typename T>
	inline std::size_t Triplet<T>::get_col() const
	{
		return COL;
	}

	template <typename T>
	inline T Triplet<T>::get_val() const
	{
		return VAL;
	}

}  // End of namespace Luna

#endif
