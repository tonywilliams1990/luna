/// \file Timer.h
/// Defines a timer class for timing methods

#ifndef TIMER_H
#define TIMER_H

#include <ctime>
#include <iostream>
#include <string>

namespace Luna
{
	/// A simple timer class for timing methods
	class Timer
	{

	private:
		clock_t START_TIME;						// Start time
		clock_t PAUSED_TIME;					// Time that the clock was paused
		bool STARTED;									// Boolean stating if the clock has been started
		bool PAUSED;									// Boolean stating if the clock has been paused

	public:

		/* ----- Constructors and Destructor ----- */

		/// Constructor
		Timer() : START_TIME( 0 ), PAUSED_TIME( 0 ), STARTED( false ), PAUSED( false )
		{}

		/// Destructor
		~Timer()
		{}

		/* ----- Methods ----- */

		/// Check if the timer has started
		/// \return true if the time has been started and false otherwise
		bool is_started();

		/// Check if the timer is stopped
		/// \return true if the timer is stopped and false otherwise
		bool is_stopped();

		/// Check if the timer is paused
		/// \return true if the timer is paused and false otherwise
		bool is_paused();

		/// Check if the timer is active
		/// \return true if the timer is still active and false otherwise
		bool is_active();

		/// Pause the timer
		void pause();

		/// Resume the timer
		void resume();

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
		double get_time() const;

		/// Return the number of clock ticks
		/// \return The number of clock ticks
		clock_t get_ticks();

	}; // End of Timer class

	/* ----- Inline definitions ----- */

	inline bool Timer::is_started()
	{
		return STARTED;
	}

	inline bool Timer::is_stopped()
	{
		return !STARTED;
	}

	inline bool Timer::is_paused()
	{
		return PAUSED;
	}

	inline bool Timer::is_active()
	{
		return !PAUSED & STARTED;
	}

	inline void Timer::pause()
	{
		if( PAUSED || !STARTED )
		{
			return;
		}
		PAUSED = true;
		PAUSED_TIME = clock();
	}

	inline void Timer::resume()
	{
		if( !PAUSED )
		{
			return;
		}
		PAUSED = false;
		START_TIME += clock() - PAUSED_TIME;
	}

	inline void Timer::start()
	{
		if( STARTED )
		{
			return;
		}
		STARTED = true;
		PAUSED = false;
		START_TIME = clock();
	}

	inline void Timer::stop()
	{
		STARTED = false;
	}

	inline void Timer::reset()
	{
		PAUSED = false;
		START_TIME = clock();
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

	inline double Timer::get_time() const
	{
		Timer temp( *this );
		return 1.e3 * temp.get_ticks() / CLOCKS_PER_SEC;
	}

	inline clock_t Timer::get_ticks()
	{
		if( !STARTED )
		{
			return 0;
		}

		if( PAUSED )
		{
			return PAUSED_TIME - START_TIME;
		}

		return clock() - START_TIME;
	}


}	// End of namespace Luna

#endif
