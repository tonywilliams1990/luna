/// \file Utility.h
/// A file specifying some utility functions that are sometimes useful
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <unistd.h>


namespace Luna
{
  namespace Utility
  {
    /// Return an integer value as a string
    /// \param val Integer value
    /// \return The integer as a string
    std::string stringify( const int &val );

    /// Return a double value as a string
    /// \param val Double value
    /// \param p Precision to be displayed in the string
    /// \param str Floatfield format flag i.e. (fixed or scientific)
    /// \return The double as a string
    std::string stringify( const double &val, int p, std::string str="" );

    /// Check if a file exists
    /// \param name The name of the file
    /// \return true if a file exists with that name otherwise false
    bool file_exists( const std::string& name );

  } // End of namespace Utility

  std::string Utility::stringify( const int &val )
  {
    std::stringstream temp;
    temp << val;
    return temp.str();
  }

  std::string Utility::stringify( const double &val, int p, std::string str )
  {
    std::stringstream temp;
    temp.precision( p );
    if ( str == "fixed" ) {
      temp << std::fixed << val;
    }
    else if ( str == "scientific" ) {
      temp << std::scientific << val;
    }
    else {
      temp << val;
    }
    return temp.str();
  }

  bool Utility::file_exists( const std::string& name )
  {
    std::ifstream f( name.c_str() );
    return f.good();
  }

} // End of namespace Luna
