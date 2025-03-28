#ifndef EVENT_HPP_
#define EVENT_HPP_

namespace simulation
{

/// Event defines a framework of continuous or discrete recurrent
/// event in terms of event duration, event rate and event action.
/// This class performs two tasks: the first task is to calculate
/// event rate based on its internal mechanism, and the second task
/// is to provide an interface to define the work that an event
/// actually does when it occurs.

class Event
{
  protected:

	/// The time duration of event occurrence.
	double duration;

	/// The recurrence rate of event.
	double rate;

	/// Event should also have a variable representing its states.
	/// This state variable is modified by the action function.
	/// The change of this state variable can be either continuous
	/// or discrete. The actual definition of state variable is up
	/// to the derived class which implements the virtual action
	/// function.

  protected:

	/// This function returns the event rate calculated using
	/// internal mechanism and the actual implementation can be
	/// redefined.
	virtual double compute_rate();

  public:

	Event(double r, double dur = 0);

	virtual ~Event() throw();

	/// This get function returns the cached event rate.
	double get_rate();

	/// This set function sets event rate to a given value.
	void set_rate(double r);

	/// This update function updates event rate by using the
	/// internal mechanism provided by compute_rate.
	void update_rate();

	/// This function fulfills the actual task of event.
	virtual void action() = 0;
};

}

#endif /*EVENT_HPP_*/
