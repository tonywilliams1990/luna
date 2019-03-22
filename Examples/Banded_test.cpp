/// \file  Banded_test.cpp
/// \ingroup Examples
/// TODO description / rename file

#include "Luna/Core"
#include "Luna/Sparse"

using namespace std;
using namespace Luna;

int main()
{
  cout << "--------------------- Banded Matrices --------------------" << endl;

  std::size_t n( 4 );
  Matrix<double> dense( n, n, 0.0 );
  dense.fill_band( -2, 1.0 );
  dense.fill_band( -1, 2.0 );
  dense.fill_band(  0, 3.0 );
  dense.fill_band(  1, 4.0 );

  cout << " dense = " << dense << endl;

  // Create a BandedMatrix
  std::size_t bands_below( 2 );
  std::size_t bands_above( 1 );
  BandedMatrix<double> band( dense, bands_below, bands_above );

  cout << " band = " << band.compact() << endl;
  Vector<double> b( n, 1.0 );
  cout << " b = " << b << endl;
  Vector<double> x;
  x = band.solve( b );
  cout << " x = " << x << endl;
  cout << " det = " << band.det() << endl;

  cout << " band.size() = " << band.size() << endl;
  cout << " band.size_below() = " << band.size_below() << endl;
  cout << " band.size_above() = " << band.size_above() << endl;

  /*typedef std::complex<double> cmplx;
  Matrix<cmplx> dense_cmplx( n, n, 0.0 );
  dense_cmplx.fill_band( -2, 1.0 );
  dense_cmplx.fill_band( -1, 2.0 );
  dense_cmplx.fill_band(  0, cmplx( 3.0, 1.0 ) );
  dense_cmplx.fill_band(  1, 4.0 );
  cout << " dense_cmplx = " << dense_cmplx << endl;

  BandedMatrix<cmplx> band_cmplx( dense_cmplx, bands_below, bands_above );
  cout << " band_cmplx = " << band_cmplx.compact() << endl;
  Vector<cmplx> b_cmplx( n, 1.0 );
  cout << " b_cmplx = " << b_cmplx << endl;
  Vector<cmplx> x_cmplx;
  x_cmplx = band_cmplx.solve( b_cmplx );
  cout << " x_cmplx = " << x_cmplx << endl;
  cout << " det_cmplx = " << band_cmplx.det() << endl;*/

  cout << "--- FINISHED ---" << endl;
}
