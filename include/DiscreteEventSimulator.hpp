#ifndef STOCHASTICPROCESSMULATOR_HPP_
#define STOCHASTICPROCESSMULATOR_HPP_

#include <EventSimulator.hpp>
#include <DiscreteEvent.hpp>

namespace simulation
{

/// DiscreteEventSimulator simulates the time-dependent behavior of
/// a dynamic system comprised by multiple interactive events using
/// Gillespie's Frist Method. The composition of this system is fully
/// dynamic due to the creation of new events, and the modification
/// and destruction of existing events.
///
/// This class simulates multiple dynmaic events of the system by
/// sampling these events according to their rates and updating the
/// selected event at each step. If the selected event interacts with
/// other events in the system, all affected events will be updated
/// accordingly. A sorted list of all events based on their waiting
/// time is maintained to speed up the sampling process.

class DiscreteEventSimulator : public EventSimulator
{
  private:

	/// The flag to indicate that the action of an event modifies
	/// itself.
	bool self_modify_flag;

	/// The flag to indicate that the action of an event destroys
	/// itself.
	bool self_destroy_flag;

	/// The pool of multiple dynamic events.
	DiscreteEvents events;

	/// The list of the pointers to sorted events.
	DiscreteEvent_iterators sorted_events;

  private:

	/// Set the state of all events.
	void set_event_state(bool s);

	/// Modification and destruction functions updates the affected
	/// events based on the lists of modified and destroyed events
	/// associated with each event, so the behaviors of these two
	/// functions are fixed.

	/// Modify existing events affected by an event action based on
	/// the list of modified events of this event.
	void modify(DiscreteEvent_iterator event_ptr);

	/// Destroy existing events affected by an event action based on
	/// the list of modified events of this event.
	void destroy(DiscreteEvent_iterator event_ptr);

	/// Trigger the action of an event specified by the simulator.
	void trigger(DiscreteEvent_iterator event_ptr);

	/// Update the events affected by the action of current event by
	/// destroying exisiting events, modifying exisiting events and
	/// creating new events.
	void synchronize(DiscreteEvent_iterator event_ptr);

	/// This function executes an event specified by the simulator,
	/// updates the recurrence rate and the waiting period of all
	/// affected events, and return the waiting period of current
	/// executed event before its action is executed.
	double exec(DiscreteEvent_iterator event_ptr);

	/// Bubble sorting algorithm.
	void bubble_sort(DiscreteEvent_iterators& events);

	/// This function finds the position to insert an event in order
	/// to maintain an event list sorted from the minimum period to
	/// the maximum one.
	DiscreteEvent_iterator_iterator_iterator locate_event_position(DiscreteEvent_iterator_iterator event_ptr_ptr, DiscreteEvent_iterator_iterator_iterator begin_event_ptr_ptr_ptr, DiscreteEvent_iterator_iterator_iterator end_event_ptr_ptr_ptr);

	/// This function moves the given event, rate, period and state
	/// to a new position in their lists, as indicated by the sorted
	/// period list.
	DiscreteEvent_iterator_iterator relocate_event(DiscreteEvent_iterator_iterator event_insert_ptr_ptr, DiscreteEvent_iterator_iterator event_ptr_ptr);

	/// This function executes the event with the minimum period and
	/// keep the event list sorted from the minimum waiting period
	/// to the maximum one.
	/// Prerequisites:
	/// The initial event list must not be empty and must be sorted
	/// before the first call to this function.
	double update();

  protected:

	/// Creation and connection functions use the interaction
	/// relationships among events to create new events and update
	/// the connections of existing events. Since such relationships
	/// are only available when derived classes are defined, these
	/// two functions are set as virtual function because they may
	/// be overriden in derived classes.

	/// Create new events triggered by an event action. This
	/// function should be redefined based on event interactions.
	virtual void create(DiscreteEvent_iterator event_ptr);

	/// Update the interaction relationships of current event with
	/// the events its action affects by using those interaction
	/// functions provided to update three association lists of
	/// current event. This function should be redefined based on
	/// event interactions.
	virtual void connect(DiscreteEvent_iterator event_ptr);

	/// Define a series of functions to implement the creation,
	/// modification and destroy of event and event interaction.
	/// Some functions will be implemented by derived class and
	/// all the functions are meant for internal use.

	/// Create new events.

	/// Extra work needs to do before adding a new event from
	/// current event.
	virtual void pre_add_event(DiscreteEvent_iterator event_ptr, DiscreteEvent* added_event);

	/// Add a new event into the event list from current event.
	DiscreteEvent_iterator add_event(DiscreteEvent_iterator event_ptr, DiscreteEvent* added_event);

	/// Extra work needs to do after adding a new event from
	/// current event.
	virtual void post_add_event(DiscreteEvent_iterator event_ptr, DiscreteEvent_iterator added_event_ptr);

	/// Extra work needs to do before self adding a new event.
	virtual void pre_add_event(DiscreteEvent* event);

	/// Self add a new event into the event list.
	DiscreteEvent_iterator add_event(DiscreteEvent* event);

	/// Extra work needs to do after self adding a new event.
	virtual void post_add_event(DiscreteEvent_iterator event_ptr);

	/// Modify existing events.

	/// Extra work needs to do before updating an event.
	virtual void pre_update_event(DiscreteEvent_iterator event_ptr, DiscreteEvent_iterator modified_event_ptr);

	/// Update the event affected by current event.
	void update_event(DiscreteEvent_iterator event_ptr, DiscreteEvent_iterator modified_event_ptr);

	/// Update the events affected by current event.
	void update_event(DiscreteEvent_iterator event_ptr, DiscreteEvent_iterators& modified_event_ptrs);

	/// Extra work needs to do after updating an event.
	virtual void post_update_event(DiscreteEvent_iterator event_ptr, DiscreteEvent_iterator modified_event_ptr);

	/// Extra work needs to do before self updating.
	virtual void pre_update_event(DiscreteEvent_iterator event_ptr);

	/// Self update current event.
	void update_event(DiscreteEvent_iterator event_ptr);

	/// Extra work needs to do after self updating.
	virtual void post_update_event(DiscreteEvent_iterator event_ptr);

	/// Destroy existing events.

	/// Extra work needs to do before removing an event.
	virtual void pre_remove_event(DiscreteEvent_iterator event_ptr, DiscreteEvent_iterator destroyed_event_ptr);

	/// Remove an event destroyed by current event.
	void remove_event(DiscreteEvent_iterator event_ptr, DiscreteEvent_iterator destroyed_event_ptr);

	/// Remove some events destroyed by current event.
	void remove_event(DiscreteEvent_iterator event_ptr, DiscreteEvent_iterators& destroyed_event_ptrs);

	/// Extra work needs to do after removing an event.
	virtual void post_remove_event(DiscreteEvent_iterator event_ptr);

	/// Extra work needs to do before self removing.
	virtual void pre_remove_event(DiscreteEvent_iterator event_ptr);

	/// Self remove current event from the event list.
	void remove_event(DiscreteEvent_iterator event_ptr);

	/// Extra work needs to do after self removing.
	virtual void post_remove_event();

	/// Event interactions.

	/// Add the connection of event 1 modifies event 2.
	void add_modification(DiscreteEvent_iterator event_ptr1, DiscreteEvent_iterator event_ptr2);

	/// Remove the connection of event 1 modifies event 2.
	void remove_modification(DiscreteEvent_iterator event_ptr1, DiscreteEvent_iterator event_ptr2);

	/// Empty the modified events of an event.
	void empty_modified_events(DiscreteEvent_iterator event_ptr);

	/// Add the connection of event 1 destroys event 2.
	void add_destruction(DiscreteEvent_iterator event_ptr1, DiscreteEvent_iterator event_ptr2);

	/// Remove the connection of event 1 destroys event 2.
	void remove_destruction(DiscreteEvent_iterator event_ptr1, DiscreteEvent_iterator event_ptr2);

	/// Empty the destroyed events of an event.
	void empty_destroyed_events(DiscreteEvent_iterator event_ptr);

	/// Remove the connection of an event with all other events.
	void remove_connection(DiscreteEvent_iterator event_ptr);

	/// These functions are inheried from EventSimulator such that
	/// derived class can make use of them.

	virtual void step_record();

	virtual void time_record();

	/// Initially sort the event list.
	/// Postrequisites:
	/// This function must be called by the overridden initialize
	/// function in derived class after the event list is initially
	/// filled such that the event list is sorted before the update
	/// function is reached, as required by the update function.
	virtual void initialize();

	virtual void finalize();

  public:

	DiscreteEventSimulator(double max_duration, size_t max_step, double record_time_interval, size_t record_step_interval);

	virtual ~DiscreteEventSimulator() throw();

	/// Get the number of events in the event list.
	size_t size() const;
};

}

#endif /*STOCHASTICPROCESSMULATOR_HPP_*/
