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
  /// A templated class to define Triplets
	template <class T>

	class Triplet
	{
		private:
			std::size_t i;										// Row index
			std::size_t j;										// Column index
			T val;														// Value

    public:

			/// Constructor for an unspecified Triplet
			Triplet(){}

			//TODO specfied constructor

			//TODO copy constructor

			/// Destructor
			~Triplet(){}

			/* ----- Operator overloading ----- */



			/* ----- Methods ----- */



  }; // End of class Triplet


}  // End of namespace Luna

#endif
