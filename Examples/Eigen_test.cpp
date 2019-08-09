/// \file Eigen_test.cpp
/// \ingroup Examples
/// TODO description / rename file

#include "Luna/Core"
#include "Luna/Eigenvalue"
#include "Luna/Sparse"

using namespace std;
using namespace Luna;

int main()
{
  cout << "----------------------- Eigensystems ----------------------" << endl;

  LinearEigensystem<double> eig_sys;

  cout << " * Compute the eigenvalues and eigenvectors of a real" << endl
       << " * symmetric matrix." << endl;

  std::size_t n( 4 );
  Matrix<double> A( n, n, 0.0 );
  // Lehmer matrix
  for ( std::size_t i = 0; i < n; i++ )
  {
    for ( std::size_t j = 0; j < n; j++ )
    {
      if ( j >= i ){
        A( i, j ) = 1. * ( i + 1 ) / ( j + 1 );
      } else {
        A( i, j ) = 1. * ( j + 1 ) / ( i + 1 );
      }
    }
  }

  cout << " * A = " << A << endl;

  eig_sys.compute_real_symmetric( A );

  cout << " * evals = " << eig_sys.eigenvalues().real() << endl;
  cout << " * evecs = " << eig_sys.eigenvector_matrix() << endl;

  Vector<double> evec_0 = eig_sys.eigenvectors()[ 0 ].real();
  double eval_0 = eig_sys.eigenvalues().real()[ 0 ];


  cout << " * A * evec_0 - eval_0 * evec_0 = ["
       << A * evec_0 - eval_0 * evec_0 << "]^T" << endl;

  cout << " * Computer the eigenvalues and eigenvectors of a real " << endl
       << " * symmetric tridiagonal matrix with main diagonal " << endl;
  Vector<double> diag( 4, 1.0 );
  diag[ 0 ] = diag[ 2 ] = 2.0;
  cout << " * diag = " << diag << endl;
  cout << " * and off diagonals " << endl;
  Vector<double> offdiag( 3, -1.0 );
  offdiag[ 1 ] = - 2.0; offdiag[ 2 ] = - 3.0;
  cout << " * offdiag = " << offdiag << endl;
  eig_sys.compute_real_symmetric_tridiagonal( diag, offdiag );
  cout << " * evals = " << eig_sys.eigenvalues().real() << endl;
  cout << " * evecs = " << eig_sys.eigenvector_matrix() << endl;

  cout << " * Compute the eigenvalues of a real nonsymmetric matrix. " << endl;
  Matrix<double> B( 3, 3, 1.0 );
  B( 0, 1 ) = B( 1, 2 ) = B( 2, 0 ) = 2.0;
  B( 0, 2 ) = B( 1, 0 ) = B( 2, 1 ) = 3.0;
  cout << " * B = " << B << endl;

  eig_sys.compute_real( B, false );
  cout << " * evals = " << eig_sys.eigenvalues() << endl;


  cout << "--- FINISHED ---" << endl;
}
