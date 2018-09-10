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

  for (std::size_t i = 0; i < a_vec.size(); ++i )
  {
    a_vec[ i ] = i;
    b_vec[ i ] = i * i;
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

  b_vec += 1;

  cout << "a_vec = " << a_vec << endl;
  cout << "b_vec = " << b_vec << endl;
  cout << "c_vec = " << c_vec << endl;

  //throw Error( "This is an error!" );


  cout << "FINISHED" << endl;

}
