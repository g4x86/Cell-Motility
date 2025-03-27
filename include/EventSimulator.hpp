#ifndef EVENTSIMULATOR_HPP_
#define EVENTSIMULATOR_HPP_

#include <cstddef>

namespace simulation
{

/// EventSimulator is a framework of simulating time-dependent
/// continuous or discrete events based on Monte Carlo Simulation.

class EventSimulator
{
  protected:

	/// The status of Monte Carlo simulator:
	///  0: simulation is terminated when maximum loop step is reached.
	///  1: simulation is terminated when maximum time duration is reached.
	///  2: simulation is terminated when maximum time slice is reached.
	int status;

	double time_moment, time_slice, max_time_duration;

	double record_time_moment_interval, last_record_time_moment;

	size_t loop_step, max_loop_step;

	size_t record_loop_step_interval, last_record_loop_step;

  protected:

	/// This function records moment-wise simulation information
	/// every given moment interval and records step-wise simulation
	/// information every given step interval.
	void record();

	/// This function records moment-wise simulation information.
	virtual void time_record();

	/// This function records step-wise simulation information.
	virtual void step_record();

	/// This function does all necessary initializing work before
	/// the simulation loop starts.
	virtual void initialize();

	/// This function does all necessary cleanup work after the
	/// simulation loop finishes, including making a final call
	/// to the time_record function and setting the last record
	/// variables.
	virtual void finalize();

	/// This function does actual simulation work and returns either
	/// the time step of continuous processes or the waiting period
	/// of stochastic processes. This value is set as the time slice
	/// of current simulation step.
	///
	/// This function must be defined in derived simulator class.
	virtual double update() = 0;

	/// The simulator constructor for continuous-time event.
	EventSimulator(double max_duration, double time_step, double record_time_interval, size_t record_step_interval);

	/// The simulator constructor for discrete-time event.
	EventSimulator(double max_duration, size_t max_step, double record_time_interval, size_t record_step_interval);

	virtual ~EventSimulator() throw();

  public:

	/// This the interface function between simulator and external
	/// program.
	void run();

	/// This function returns the exit status of simulator.
	int get_status() const throw();
};

}

#endif /*EVENTSIMULATOR_HPP_*/
