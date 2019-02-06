/// \file  Sparse_test.cpp
/// \ingroup Examples
/// TODO

#include "Luna/Core"
#include "Luna/Sparse"

using namespace std;
using namespace Luna;

int main()
{
  cout << "--------------------- Sparse Matrices --------------------" << endl;

  Vector<double> val( 9 );
  Vector<std::size_t> row_ind( 9 );
  Vector<std::size_t> col_ptr( 6 );

  val[ 0 ] = 3.; val[ 1 ] = 4.; val[ 2 ] = 7.; val[ 3 ] = 1.; val[ 4 ] = 5.;
  val[ 5 ] = 2.; val[ 6 ] = 9.; val[ 7 ] = 6.; val[ 8 ] = 5.;

  row_ind[ 0 ] = 0; row_ind[ 1 ] = 1; row_ind[ 2 ] = 2; row_ind[ 3 ] = 0;
  row_ind[ 4 ] = 2; row_ind[ 5 ] = 0; row_ind[ 6 ] = 2; row_ind[ 7 ] = 4;
  row_ind[ 8 ] = 4;

  col_ptr[ 0 ] = 0; col_ptr[ 1 ] = 1; col_ptr[ 2 ] = 3; col_ptr[ 3 ] = 5;
  col_ptr[ 4 ] = 8; col_ptr[ 5 ] = 9;



  SparseMatrix<double> sparse( 5, 5, val, row_ind, col_ptr );


  cout << " * val       = " << sparse.val() << endl;
  cout << " * row_index = " << sparse.row_index() << endl;
  cout << " * col_start = " << sparse.col_start() << endl;
  cout << " * nonzeros  = " << sparse.nonzero() << endl;

  //Vector<std::size_t> col_ind;
  //col_ind = sparse.col_index();
  //cout << " * col_index = " << col_ind << endl;
  //cout << " * col_start = " << sparse.col_start_from_index( col_ind ) << endl;

  sparse.insert( 1, 3, 2.0 );
  sparse.insert( 1, 0, 5.0 );

  cout << " * val       = " << sparse.val() << endl;
  cout << " * row_index = " << sparse.row_index() << endl;
  cout << " * col_start = " << sparse.col_start() << endl;
  cout << " * nonzeros  = " << sparse.nonzero() << endl;

  sparse.scale( 2.0 );
  cout << " * val       = " << sparse.val() << endl;

  cout << "--- FINISHED ---" << endl;
}
