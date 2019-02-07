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
  triplets( Td( 3, 2, 1.0 ) );

  SparseMatrix<double> sparse_trip( 5, 5, triplets );
  cout << " * val       = " << sparse_trip.val() << endl;
  cout << " * row_index = " << sparse_trip.row_index() << endl;
  cout << " * col_start = " << sparse_trip.col_start() << endl;
  cout << " * nonzeros  = " << sparse_trip.nonzero() << endl;

  Vector<double> b { 1, 2, 3, 4, 5 };
  Vector<double> x( 5, 0.0 );

  int iter( 0 );
  double err( 0.0 );
  sparse_trip.solve_bcg( b, x, 1, 1e-8, 100, iter, err );

  cout << " * x = " << x << endl;
  cout << " * iter = " << iter << endl;
  cout << " * err = " << scientific << err << endl;

  //TODO need to test a complex sparse matrix

  typedef std::complex<double> cmplx;
  typedef Triplet<cmplx> Tc;

  Vector<Tc> trips;

  trips( Tc( 0, 0, cmplx( 1.0, 1.0 ) ) );
  trips( Tc( 0, 1, cmplx( -1.0, 0.0 ) ) );
  trips( Tc( 1, 0, cmplx( 1.0, -1.0 ) ) );
  trips( Tc( 1, 1, cmplx( 1.0, 1.0 ) ) );

  SparseMatrix<cmplx> sparse_cmplx( 2, 2, trips );
  Vector<cmplx> c{ 1.0, 1.0 };
  Vector<cmplx> y( 2, 0.0 );
  cout << " * c = " << c << endl;
  sparse_cmplx.solve_bcg( c, y, 1, 1e-8, 100, iter, err );
  cout << " * y = " << y << endl;
  cout << " * B * y = " << sparse_cmplx.multiply( y ) << endl;

  //cout << trips << endl;

  cout << "--- FINISHED ---" << endl;
}
