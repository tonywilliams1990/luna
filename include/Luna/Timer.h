/// \file Timer.h
/// Defines a timer class for timing methods

#ifndef TIMER_H
#define TIMER_H

#include <iostream>
#include <chrono>
#include <ctime>
#include <cmath>

namespace Luna
{
  /// A simple timer class for timing methods
  class Timer
  {

  private:
    std::chrono::time_point<std::chrono::system_clock> START_TIME;
    std::chrono::time_point<std::chrono::system_clock> END_TIME;
    bool                                               RUNNING = false;

  public:

    /* ----- Constructors and Destructor ----- */

		/// Constructor
		Timer()
		{}

		/// Destructor
		~Timer()
		{}

		/* ----- Methods ----- */

    /// Start the timer
    void start();

    /// Stop the timer
    void stop();

    /// Reset the timer
		void reset();

    /// Output the time to the screen
		void print() const;

    /// Output the time to the screen with a specified message
		/// \param message The message to be printed before the time
		void print( const std::string& message ) const;

    /// Return the time in milliseconds
		/// \return The time in ms
    double get_time();

    /// Return the time in seconds
    /// \return The time in s
    double get_time_in_seconds();

  }; // End of Timer class

  inline void Timer::start()
  {
      START_TIME = std::chrono::system_clock::now();
      RUNNING = true;
  }

  inline void Timer::stop()
  {
      END_TIME = std::chrono::system_clock::now();
      RUNNING = false;
  }

  inline void Timer::reset()
	{
		RUNNING = true;
		START_TIME = std::chrono::system_clock::now();
	}

  inline void Timer::print() const
	{
		std::cout.precision(4);
		Timer temp( *this );
		const double elapsed_time_in_ms( temp.get_time() );
		if ( elapsed_time_in_ms > 1000 )
    {
    	std::cout << "  * TOTAL CPU time taken = " << elapsed_time_in_ms / 1000.
								<< " s\n";
    }
    else
    {
    	std::cout << "  * TOTAL CPU time taken = " << elapsed_time_in_ms
								<< " ms\n";
    }
	}

  inline void Timer::print( const std::string& message ) const
	{
		std::cout.precision(4);
		Timer temp( *this );
		const double elapsed_time_in_ms( temp.get_time() );
		if ( elapsed_time_in_ms > 1000 )
    {
    	std::cout << "  * " << message << " = " << elapsed_time_in_ms / 1000.
								<< " s\n";
    }
    else
    {
    	std::cout << "  * " << message << " = " << elapsed_time_in_ms << " ms\n";
    }
	}

  inline double Timer::get_time()
  {
      std::chrono::time_point<std::chrono::system_clock> endTime;

      if( RUNNING )
      {
          endTime = std::chrono::system_clock::now();
      }
      else
      {
          endTime = END_TIME;
      }

      return std::chrono::duration_cast<std::chrono::milliseconds>(endTime -
             START_TIME).count();
  }

  inline double Timer::get_time_in_seconds()
  {
      return get_time() / 1000.0;
  }

} // End of namespace Luna

#endif
