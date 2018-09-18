/// \file  Vector_algebra.cpp
/// \ingroup Examples
/// Some simple vector algebra to check the basic functionality

#include "Luna/Core"

using namespace std;
using namespace Luna;

int main()
{
  cout << "----- Vector algebra -----" << endl;

  Vector<double> a_vec( 10, 1.0 );
  Vector<double> b_vec( a_vec );                     // Copy constructor
  Vector< std::complex<double> > cmplx_vec( 10, 1.0 );

  for (std::size_t i = 0; i < a_vec.size(); ++i )
  {
    a_vec[ i ] = i;
    b_vec[ i ] = i * i;
    cmplx_vec[ i ] = std::complex<double> ( i, 2 * i );
  }

  Vector<double> c_vec;
  c_vec = - b_vec;                                      // Copy assignment
  a_vec = b_vec - c_vec;
  a_vec = a_vec / 3;
  c_vec /= 3;

  c_vec.push_back( 100 );
  b_vec.resize( 11 );
  a_vec.clear();
  a_vec = c_vec.abs();
  b_vec.swap( 0, 5 );

  b_vec += - 1;
  a_vec.reverse();
  c_vec.assign( 11, 1.0 );
  a_vec = cmplx_vec.real();
  a_vec.powspace( 0, 1, 11, 1.2 );

  cout << "a_vec = " << a_vec << endl;
  cout << "b_vec = " << b_vec << endl;
  cout << "c_vec = " << c_vec << endl;
  cout << "cmplx_vec = " << cmplx_vec << endl;
  cout << "cmplx_vec.conjugate() = " << cmplx_vec.conjugate() << endl;

  cout << "b_vec.sum() = " << b_vec.sum( ) << endl;
  cout << "b_vec.dot(c_vec) = " << b_vec.dot(c_vec) << endl;

  a_vec.output( "./a_vec.dat", 15 );
  c_vec.read( "./a_vec.dat", true );
  cout << "c_vec = " << c_vec << endl;
  cout << "b_vec.max() = " << b_vec.max() << endl;
  cout << "b_vec.min() = " << b_vec.min() << endl;
  cout << "b_vec.max_index() = " << b_vec.max_index() << endl;
  cout << "b_vec.min_index() = " << b_vec.min_index() << endl;
  cout << "b_vec.square() = " << b_vec.square() << endl;
  cout << "b_vec.power() = " << b_vec.power( 3 ) << endl;

  cout << "b_vec.norm_1() = " << b_vec.norm_1() << endl;
  cout << "b_vec.norm_2() = " << b_vec.norm_2() << endl;
  cout << "b_vec.norm_p( p = 3.14 ) = " << b_vec.norm_p( 3.14 ) << endl;
  cout << "b_vec.norm_inf() = " << b_vec.norm_inf() << endl;

  cout << "FINISHED" << endl;

}
