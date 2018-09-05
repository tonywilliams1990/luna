/*
  Error.h - Defines a class for simple error handling.
*/

#ifndef ERROR_H
#define ERROR_H

#include <string>
#include <iostream>
#include <stdexcept>

namespace Luna
{

  /// A generic runtime exception
  class Error : public std::runtime_error
  {
  public:

    /// Constructor
    /// \param problem Error string to be produced
    Error( const std::string &problem ) : std::runtime_error( problem )
    {
      error_header();
      std::cout << problem << std::endl;
    }

  private:

    void error_header()
    {
      std::cout << "------------------------------------------" << std::endl;
      std::cout << " Error: A Luna routine has had a problem! " << std::endl;
      std::cout << "------------------------------------------" << std::endl;
    }

  };	// End of class Error
}  // End of namespace Luna

#endif
