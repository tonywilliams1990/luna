/// \file  Sparse_test.cpp
/// \ingroup Examples
/// TODO description / rename file

#include "Luna/Core"
#include "Luna/Sparse"

using namespace std;
using namespace Luna;

int main()
{
  cout << "--------------------- Sparse Matrices --------------------" << endl;

  Vector<double> val( { 3., 4., 7., 1., 5., 2., 9., 6., 5. } );
  Vector<std::size_t> row_ind( { 0, 1, 2, 0, 2, 0, 2, 4, 4 } );
  Vector<std::size_t> col_ptr( { 0, 1, 3, 5, 8, 9 } );

  SparseMatrix<double> sparse( 5, 5, val, row_ind, col_ptr );


  cout << " * val       = " << sparse.val() << endl;
  cout << " * row_index = " << sparse.row_index() << endl;
  cout << " * col_start = " << sparse.col_start() << endl;
  cout << " * nonzeros  = " << sparse.nonzero() << endl;

  sparse.insert( 1, 3, 2.0 );
  sparse.insert( 1, 0, 5.0 );

  cout << " * val       = " << sparse.val() << endl;
  cout << " * row_index = " << sparse.row_index() << endl;
  cout << " * col_start = " << sparse.col_start() << endl;
  cout << " * nonzeros  = " << sparse.nonzero() << endl;

  typedef Triplet<double> Td;
  Vector<Td> triplets { Td( 0, 3, 2.0 ), Td( 2, 3, 9.0 ), Td( 4, 3, 6.0 ),
                        Td( 4, 4, 5.0 ), Td( 0, 0, 3.0 ), Td( 1, 1, 4.0 ) };

  triplets( Td( 2, 1, 7.0 ) );
  triplets( Td( 0, 2, 1.0 ) );
  triplets( Td( 2, 2, 5.0 ) );

  SparseMatrix<double> sparse_trip( 5, 5, triplets );
  cout << " * val       = " << sparse_trip.val() << endl;
  cout << " * row_index = " << sparse_trip.row_index() << endl;
  cout << " * col_start = " << sparse_trip.col_start() << endl;
  cout << " * nonzeros  = " << sparse_trip.nonzero() << endl;

  Vector<double> x { 1, 2, 3, 4, 5 };
  cout << " * x = " << x << endl;

  Vector<double> b;
  b = sparse_trip.multiply( x );
  cout << " * A * x = " << b << endl;
  b = sparse_trip.transpose_multiply( x );
  cout << " * A^T * x = " << b << endl;

  SparseMatrix<double> transp;
  transp = sparse_trip.transpose();
  cout << " * val       = " << transp.val() << endl;
  cout << " * row_index = " << transp.row_index() << endl;
  cout << " * col_start = " << transp.col_start() << endl;
  cout << " * nonzeros  = " << transp.nonzero() << endl;

  b = transp.multiply( x );
  cout << " * A^T * x = " << b << endl;

  cout << "--- FINISHED ---" << endl;
}
