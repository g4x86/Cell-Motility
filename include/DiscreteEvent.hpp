#ifndef DISCRETEEVENT_HPP_
#define DISCRETEEVENT_HPP_

#include <list>
#include <Event.hpp>

namespace simulation
{

class DiscreteEvent;
typedef std::list<DiscreteEvent*> DiscreteEvents;
typedef DiscreteEvents::iterator DiscreteEvent_iterator;
typedef std::list<DiscreteEvent_iterator> DiscreteEvent_iterators;
typedef DiscreteEvent_iterators::iterator DiscreteEvent_iterator_iterator;
typedef std::list<DiscreteEvent_iterator_iterator> DiscreteEvent_iterator_iterators;
typedef DiscreteEvent_iterator_iterators::iterator DiscreteEvent_iterator_iterator_iterator;

class DiscreteEvent : public Event
{
  private:

	/// The waiting period between event recurrence.
	double period;

	/// The boolean flag of the change of an event.
	/// false:  an event is not changed.
	/// true:   an event is changed.
	bool state;

	/// Three lists of events with which this effect is associated.
	/// The modified event list and the destroyed event list reflect
	/// two types of events affected by the action of this event.
	/// Modified events is the events whose rates and periods are
	/// changed by the action of this event, and destroyed events
	/// are the events destroyed by the action of this event. The
	/// actuator event list is a quick reference to those events
	/// which modify or destroy this event.

	/// The list of the events to be modified by the action of
	/// this event which may also include this event itself to
	/// indicate self modification.
	DiscreteEvent_iterators modified_events;

	/// The list of the events to be destroyed by the action of
	/// this event which may also include this event itself to
	/// indicate self destruction.
	DiscreteEvent_iterators destroyed_events;

	/// The list of the events whose actions affect this event.
	DiscreteEvent_iterators actuator_events;

	/// The list of modified events and the list of destroyed
	/// events should not contain the same events, which means
	/// an event should either be modified or be destroyed, but
	/// not both. Therefore all events in the actuator list are
	/// guaranteed to be unique.

	/// The pointer to corresponding event in the sorted event
	/// list in the simulator. This pointer provides a fast
	/// access to this event itself in the simulator.
	DiscreteEvent_iterator_iterator sorted_event;

  private:

	/// This function calculates waiting period by sampling the
	/// exponential probability density function of the waiting
	/// period determined by event rate, based on Gillespie's
	/// algorithm. The algorithm is fixed and cannot be redefined.
	double compute_period();

  protected:

	/// This function returns the event rate calculated using
	/// internal mechanism and the actual implementation can be
	/// redefined.
	virtual double compute_rate();

  public:

	virtual ~DiscreteEvent() throw();

	DiscreteEvent(double r = 0, double dur = 0);

	/// These functions provides access to event period.

	double get_period();

	/// This update function updates event period by using the
	/// internal mechanism provided by compute_period.
	void update_period();

	/// This function calls both update_rate and update_period
	/// in order to update both rate and period.
	void update();

	/// These functions provides access to event state.

	bool get_state();

	void set_state(bool s);

	/// These functions defines various comparison operations based
	/// on event period.

	DiscreteEvent& operator=(const DiscreteEvent& de);

	bool operator<(const DiscreteEvent& de) const;

	bool operator<=(const DiscreteEvent& de) const;

	bool operator>(const DiscreteEvent& de) const;

	bool operator>=(const DiscreteEvent& de) const;

	bool operator==(const DiscreteEvent& de) const;

	bool operator!=(const DiscreteEvent& de) const;

	/// This pure virtual function performs the actual tasks of an event and
	/// must be defined in derived class. The tasks of an event may include
	/// creating new events, modifying and/or destroying exisiting events.
	/// There are some guideline rules to follow as listed below:
	/// 1) To create new events, this function should use add_event function
	/// provided by the simulator to create new events and add them into the
	/// event list of the simulator, and should use add_modification and
	/// add_destruction functions provided by the simulator to create specific
	/// types of interactions based on the interaction relationships between
	/// current event and the new events it creates. To perform any extra tasks
	/// before or after new events are added to the event list of the simulator,
	/// please redefine the virtual functions pre_add_event and post_add_event
	/// provided by the simulator.
	/// 2) To modify existing events before or after their recurrence rates are
	/// changed, please redefine the virtual functions pre_modify_event and
	/// post_remove_event provided by the simulator to perform the tasks, and
	/// the modify function provided by the simulator takes care of the rest.
	/// So do not modify event rate explicitly inside this function. And since
	/// it is also possible to perform any actual tasks before event rates are
	/// changed inside this function, it is recommended to use those virtual
	/// functions provided by the simulator to perform these tasks in order to
	/// maintain clear consistency.
	/// 3) To destroy existing events, nothing needs to be done inside this
	/// function because everything is handled by the destroy function provided
	/// by the simulator. However if there are any extra work needs to be done
	/// before or after current event destroys exisiting events, then redefine
	/// the virtual functions pre_remove_event and post_remove_event provided by
	/// the simulator to perform the tasks. For the same consistency reason,
	/// please use those virtual functions provided by the simulator to perform
	/// these tasks.
	virtual void action() = 0;

	/// DiscreteEventSimulator have direct access to all data
	/// members of this class to improve efficiency.
	friend class DiscreteEventSimulator;
};

}

#endif /*DISCRETEEVENT_HPP_*/
