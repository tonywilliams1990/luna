/// \file  Vector_algebra.cpp
/// \ingroup Examples
/// Some simple vector algebra to check the basic functionality

#include "Luna/Core"

using namespace std;
using namespace Luna;

int main()
{
  cout << "----- Vector algebra -----" << endl;

  Vector<double> a( 4, 0.0 );                       // Specified constructor
  Vector<double> b( a );                            // Copy constructor
  Vector<double> c;                                 // Empty constructor

  // Fill the Vectors
  a[ 0 ] = 1.0; a[ 1 ] = 3.0; a[ 2 ] = 1.0; a[ 3 ] = - 1.0;
  b[ 0 ] = 2.0; b[ 1 ] = 1.0; b[ 2 ] = 3.0; b[ 3 ] = 4.0;
  cout << " * a = " << a << endl;
  cout << " * b = " << b << endl;

  // Do some basic algebra
  c = a + b;
  cout << " * c = a + b = " << c << endl;
  c = b - a;
  cout << " * c = b - a = " << c << endl;
  c = 3 * a;
  cout << " * c = 3 * a = " << c << endl;

  // Add some new elements and swap them around
  a.push_back( 2.0 );
  b.resize( 5 );
  b[ 4 ] = - 3.0;
  b.swap( 3, 4 );
  cout << " * a = " << a << endl;
  cout << " * b = " << b << endl;

  // Magnitude of the Vectors
  cout << " * |a| = " << a.norm_2() << endl;
  cout << " * |b| = " << b.norm_p( 2 ) << endl;

  // Remove elements
  a.pop_back();
  a.pop_back();
  b.resize( 3 );
  cout << " * a = " << a << endl;
  cout << " * b = " << b << endl;

  // Angle between the two Vectors
  double angle;
  angle = a.dot( b ) / ( a.norm_2() * b.norm_2() );
  angle = std::acos( angle );
  angle = angle * 180 / M_PI;
  cout << " * angle between a and b = " << angle << " degrees "<< endl;

  // Create Vectors of spaced elements
  a.linspace( 0, 1, 6 );
  b.powspace( 0, 1, 6, 1.5 );
  cout << " * Linearly spaced elements = " << a << endl;
  cout << " * Power law spaced elements (p=1.5) = " << b << endl;

  cout << "--- FINISHED ---" << endl;

}
