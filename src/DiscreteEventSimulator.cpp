#include <constants.hpp>
#include <DiscreteEventSimulator.hpp>

namespace simulation
{

DiscreteEventSimulator::DiscreteEventSimulator(double max_duration, size_t max_step, double record_time_interval, size_t record_step_interval) : EventSimulator(max_duration, max_step, record_time_interval, record_step_interval)
{
	self_modify_flag = false;
	self_destroy_flag = false;
	events.clear();
	sorted_events.clear();
}

DiscreteEventSimulator::~DiscreteEventSimulator() throw()
{
	for(DiscreteEvent_iterator event_ptr = events.begin(); event_ptr != events.end(); ++event_ptr)
	{
		if((*event_ptr) != 0) delete (*event_ptr);
	}
}

void DiscreteEventSimulator::set_event_state(bool s)
{
	for(DiscreteEvent_iterator event_ptr = events.begin(); event_ptr != events.end(); ++event_ptr)
	{
		(*event_ptr)->state = s;
	}
}

size_t DiscreteEventSimulator::size() const
{
	return events.size();
}

void DiscreteEventSimulator::bubble_sort(DiscreteEvent_iterators& events)
{
	if(events.size() > 1)
	{
		DiscreteEvent_iterator_iterator event_ptr_ptr = ++events.begin();
		while(event_ptr_ptr != events.end())
		{
			DiscreteEvent_iterator_iterator event_begin_search_ptr_ptr = event_ptr_ptr;
			--event_begin_search_ptr_ptr;
			DiscreteEvent_iterator_iterator event_insert_ptr_ptr = event_begin_search_ptr_ptr;
			while(event_insert_ptr_ptr != events.end() && ***event_ptr_ptr < ***event_insert_ptr_ptr) --event_insert_ptr_ptr;
			if(event_insert_ptr_ptr != event_begin_search_ptr_ptr)
			{
				++event_insert_ptr_ptr;
				events.insert(event_insert_ptr_ptr, *event_ptr_ptr);
				// Update the event at the new position.
				DiscreteEvent_iterator_iterator event_new_ptr_ptr = event_insert_ptr_ptr;
				--event_new_ptr_ptr;
				(**event_new_ptr_ptr)->sorted_event = event_new_ptr_ptr;
				// Remove the event at the old position.
				DiscreteEvent_iterator_iterator event_old_ptr_ptr = event_ptr_ptr++;
				events.erase(event_old_ptr_ptr);
			}
			else ++event_ptr_ptr;
		}
	}
}

void DiscreteEventSimulator::step_record()
{
	EventSimulator::step_record();
}

void DiscreteEventSimulator::time_record()
{
	EventSimulator::time_record();
}

void DiscreteEventSimulator::initialize()
{
	EventSimulator::initialize();
	// Initially sort the sorted event list.
	bubble_sort(sorted_events);
	// Reset the state of all events.
	set_event_state(false);
}

void DiscreteEventSimulator::finalize()
{
	EventSimulator::finalize();
}

DiscreteEvent_iterator_iterator_iterator DiscreteEventSimulator::locate_event_position(DiscreteEvent_iterator_iterator event_ptr_ptr, DiscreteEvent_iterator_iterator_iterator begin_event_ptr_ptr_ptr, DiscreteEvent_iterator_iterator_iterator end_event_ptr_ptr_ptr)
{
	// Recursive dichotomy searching algorithm.
	// The searching range [begin_event, end_event] must maintain
	// monotonically increasing sorted and the iterators within
	// this range must not contain the event to be inserted later.
	DiscreteEvent_iterator_iterator_iterator target_event_ptr_ptr_ptr;
	DiscreteEvent_iterator_iterator_iterator::difference_type dist = distance(begin_event_ptr_ptr_ptr, end_event_ptr_ptr_ptr);
	if(dist == 0)
	{
		target_event_ptr_ptr_ptr = begin_event_ptr_ptr_ptr;
		if(***event_ptr_ptr > ****begin_event_ptr_ptr_ptr) ++target_event_ptr_ptr_ptr;
	}
	else if(dist == 1)
	{
		if(***event_ptr_ptr <= ****begin_event_ptr_ptr_ptr) target_event_ptr_ptr_ptr = begin_event_ptr_ptr_ptr;
		else if(***event_ptr_ptr <= ****end_event_ptr_ptr_ptr) target_event_ptr_ptr_ptr = end_event_ptr_ptr_ptr;
		else
		{
			target_event_ptr_ptr_ptr = end_event_ptr_ptr_ptr;
			++target_event_ptr_ptr_ptr;
		}
	}
	else
	{
		if(***event_ptr_ptr <= ****begin_event_ptr_ptr_ptr) target_event_ptr_ptr_ptr = begin_event_ptr_ptr_ptr;
		else if(***event_ptr_ptr > ****end_event_ptr_ptr_ptr)
		{
			target_event_ptr_ptr_ptr = end_event_ptr_ptr_ptr;
			++target_event_ptr_ptr_ptr;
		}
		else
		{
			DiscreteEvent_iterator_iterator_iterator half_event_ptr_ptr_ptr = begin_event_ptr_ptr_ptr;
			advance(half_event_ptr_ptr_ptr, dist / 2);
			if(***event_ptr_ptr <= ****half_event_ptr_ptr_ptr)
			{
				target_event_ptr_ptr_ptr = locate_event_position(event_ptr_ptr, begin_event_ptr_ptr_ptr, half_event_ptr_ptr_ptr);
			}
			else
			{
				target_event_ptr_ptr_ptr = locate_event_position(event_ptr_ptr, half_event_ptr_ptr_ptr, end_event_ptr_ptr_ptr);
			}
		}
	}
	return target_event_ptr_ptr_ptr;
}

DiscreteEvent_iterator_iterator DiscreteEventSimulator::relocate_event(DiscreteEvent_iterator_iterator event_insert_ptr_ptr, DiscreteEvent_iterator_iterator event_ptr_ptr)
{
	// Move an event to its proper position in the list to maintain
	// a sorted event list.
	// Insert this event at the new position.
	sorted_events.insert(event_insert_ptr_ptr, *event_ptr_ptr);
	// Update the event at the new position.
	DiscreteEvent_iterator_iterator event_new_ptr_ptr = event_insert_ptr_ptr;
	--event_new_ptr_ptr;
	(**event_new_ptr_ptr)->sorted_event = event_new_ptr_ptr;
	// Remove the event at the old position.
	sorted_events.erase(event_ptr_ptr);
	// Return the new position.
	return event_new_ptr_ptr;
}

double DiscreteEventSimulator::update()
{
	double minimum_waiting_period;
	size_t n_event = sorted_events.size();
	if(n_event < 1) minimum_waiting_period = DBL_INF_POSITIVE;
	else
	{
		// Execute the event with the minimum waiting period.
		DiscreteEvent_iterator exec_event_ptr = *sorted_events.begin();
		minimum_waiting_period = exec(exec_event_ptr);
		if(n_event < 2) (*exec_event_ptr)->state = false;
		else
		{
			// Filter out the changed events.
			DiscreteEvent_iterator_iterators changed_event_ptr_ptrs;
			DiscreteEvent_iterator_iterators unchanged_event_ptr_ptrs;
			for(DiscreteEvent_iterator_iterator event_ptr_ptr = sorted_events.begin(); event_ptr_ptr != sorted_events.end(); ++event_ptr_ptr)
			{
				if((**event_ptr_ptr)->get_state()) changed_event_ptr_ptrs.push_back(event_ptr_ptr);
				else unchanged_event_ptr_ptrs.push_back(event_ptr_ptr);
			}
			if(unchanged_event_ptr_ptrs.size() > 0)
			{
				// If there are unchanged events available, insert each changed
				// event into the list of unchanged events and keep the list
				// sorted.
				for(DiscreteEvent_iterator_iterator_iterator event_ptr_ptr_ptr = changed_event_ptr_ptrs.begin(); event_ptr_ptr_ptr != changed_event_ptr_ptrs.end(); ++event_ptr_ptr_ptr)
				{
					// Find the position in the unchanged event list at which current
					// changed event should be inserted.
					DiscreteEvent_iterator_iterator_iterator event_insert_ptr_ptr_ptr = locate_event_position(*event_ptr_ptr_ptr, unchanged_event_ptr_ptrs.begin(), --unchanged_event_ptr_ptrs.end());
					// Convert this position in the unchanged event list to that in
					// the sorted event.
					DiscreteEvent_iterator_iterator event_insert_ptr_ptr;
					if(event_insert_ptr_ptr_ptr != unchanged_event_ptr_ptrs.end()) event_insert_ptr_ptr = *event_insert_ptr_ptr_ptr;
					else
					{
						DiscreteEvent_iterator_iterator_iterator last_event_ptr_ptr_ptr = event_insert_ptr_ptr_ptr;
						event_insert_ptr_ptr = *(--last_event_ptr_ptr_ptr);
						++event_insert_ptr_ptr;
					}
					// Relocate current changed event to the new position in the sorted
					// event list.
					DiscreteEvent_iterator_iterator event_relocated_ptr_ptr = relocate_event(event_insert_ptr_ptr, *event_ptr_ptr_ptr);
					// Reset the state of this event.
					(**event_relocated_ptr_ptr)->state = false;
					// After relocation, the content of event_ptr_ptr_ptr becomes invalid,
					// but event_ptr_ptr_ptr itself is still valid, so it is safe to use
					// a 'for' loop to go through the entire list of changed events.
					// Add this relocated event to the unchanged event list.
					unchanged_event_ptr_ptrs.insert(event_insert_ptr_ptr_ptr, event_relocated_ptr_ptr);
				}
			}
			else
			{
				// If all events are changed, then use bubble sorting algorithm
				// to sort all events at one time.
				bubble_sort(sorted_events);
				// Reset the state of all events to false.
				set_event_state(false);
			}
		}
	}
	// Return the minimum waiting period.
	return minimum_waiting_period;
}

// Create new events.

void DiscreteEventSimulator::pre_add_event(DiscreteEvent_iterator event_ptr, DiscreteEvent* added_event) {}

DiscreteEvent_iterator DiscreteEventSimulator::add_event(DiscreteEvent_iterator event_ptr, DiscreteEvent* added_event)
{
	// Pre-processing.
	pre_add_event(event_ptr, added_event);
	// Processing.
	events.push_back(added_event);
	// Synchronize with the sorted event list.
	DiscreteEvent_iterator added_event_ptr = --events.end();
	sorted_events.push_back(added_event_ptr);
	(*added_event_ptr)->sorted_event = --(sorted_events.end());
	// Post-processing.
	post_add_event(event_ptr, added_event_ptr);
	return added_event_ptr;
}

void DiscreteEventSimulator::post_add_event(DiscreteEvent_iterator event_ptr, DiscreteEvent_iterator added_event_ptr) {}

void DiscreteEventSimulator::pre_add_event(DiscreteEvent* event) {}

DiscreteEvent_iterator DiscreteEventSimulator::add_event(DiscreteEvent* event)
{
	// Pre-processing.
	pre_add_event(event);
	// Processing.
	events.push_back(event);
	// Synchronize with the sorted event list.
	DiscreteEvent_iterator event_ptr = --events.end();
	sorted_events.push_back(event_ptr);
	(*event_ptr)->sorted_event = --(sorted_events.end());
	// Post-processing.
	post_add_event(event_ptr);
	return event_ptr;
}

void DiscreteEventSimulator::post_add_event(DiscreteEvent_iterator event_ptr) {}

void DiscreteEventSimulator::create(DiscreteEvent_iterator event_ptr) {}

// Modify existing events.

void DiscreteEventSimulator::pre_update_event(DiscreteEvent_iterator event_ptr, DiscreteEvent_iterator modified_event_ptr) {}

void DiscreteEventSimulator::update_event(DiscreteEvent_iterator event_ptr, DiscreteEvent_iterator modified_event_ptr)
{
	// Pre-processing.
	pre_update_event(event_ptr, modified_event_ptr);
	// Processing.
	// Update event rate.
	(*modified_event_ptr)->update();
	// Post-processing.
	post_update_event(event_ptr, modified_event_ptr);
}

void DiscreteEventSimulator::update_event(DiscreteEvent_iterator event_ptr, DiscreteEvent_iterators& modified_event_ptrs)
{
	for(DiscreteEvent_iterator_iterator modified_event_ptr_ptr = modified_event_ptrs.begin(); modified_event_ptr_ptr != modified_event_ptrs.end(); ++modified_event_ptr_ptr)
	{
		update_event(event_ptr, *modified_event_ptr_ptr);
	}
}

void DiscreteEventSimulator::post_update_event(DiscreteEvent_iterator event_ptr, DiscreteEvent_iterator modified_event_ptr) {}

void DiscreteEventSimulator::pre_update_event(DiscreteEvent_iterator event_ptr) {}

void DiscreteEventSimulator::update_event(DiscreteEvent_iterator event_ptr)
{
	// Pre-processing.
	pre_update_event(event_ptr);
	// Processing.
	// Update event rate.
	(*event_ptr)->update();
	// Post-processing.
	post_update_event(event_ptr);
}

void DiscreteEventSimulator::post_update_event(DiscreteEvent_iterator event_ptr) {}

void DiscreteEventSimulator::modify(DiscreteEvent_iterator event_ptr)
{
	// Modify all other affected events.
	self_modify_flag = false;
	DiscreteEvent_iterators filtered_events;
	DiscreteEvent_iterators& modified_events = (*event_ptr)->modified_events;
	for(DiscreteEvent_iterator_iterator event_ptr_ptr = modified_events.begin(); event_ptr_ptr != modified_events.end(); ++event_ptr_ptr)
	{
		if(*event_ptr_ptr != event_ptr) filtered_events.push_back(*event_ptr_ptr);
		else self_modify_flag = true;
	}
	update_event(event_ptr, filtered_events);
	// Modify current event when no self destruction.
	if(!self_destroy_flag && self_modify_flag) update_event(event_ptr);
}

// Destroy existing events.

void DiscreteEventSimulator::pre_remove_event(DiscreteEvent_iterator event_ptr, DiscreteEvent_iterator destroyed_event_ptr) {}

void DiscreteEventSimulator::remove_event(DiscreteEvent_iterator event_ptr, DiscreteEvent_iterator destroyed_event_ptr)
{
	if(destroyed_event_ptr != events.end())
	{
		// Pre-processing.
		pre_remove_event(event_ptr, destroyed_event_ptr);
		// Processing.
		// Disconnect this event from the events it is associated with.
		remove_connection(destroyed_event_ptr);
		// Synchronize with the sorted event list.
		sorted_events.erase((*destroyed_event_ptr)->sorted_event);
		// Remove the memory copy of this event.		
		if((*destroyed_event_ptr) != 0) delete (*destroyed_event_ptr);
		// Finally remove this event.
		events.erase(destroyed_event_ptr);
		// Post-processing.
		post_remove_event(event_ptr);
	}
}

void DiscreteEventSimulator::remove_event(DiscreteEvent_iterator event_ptr, DiscreteEvent_iterators& destroyed_event_ptrs)
{
	for(DiscreteEvent_iterator_iterator destroyed_event_ptr_ptr = destroyed_event_ptrs.begin(); destroyed_event_ptr_ptr != destroyed_event_ptrs.end(); ++destroyed_event_ptr_ptr)
	{
		remove_event(event_ptr, *destroyed_event_ptr_ptr);
	}
}

void DiscreteEventSimulator::post_remove_event(DiscreteEvent_iterator event_ptr) {}

void DiscreteEventSimulator::pre_remove_event(DiscreteEvent_iterator event_ptr) {}

void DiscreteEventSimulator::remove_event(DiscreteEvent_iterator event_ptr)
{
	if(event_ptr != events.end())
	{
		// Pre-processing.
		pre_remove_event(event_ptr);
		// Processing.
		// Disconnect this event from the events it is associated with.
		remove_connection(event_ptr);
		// Synchronize with the sorted event list.
		sorted_events.erase((*event_ptr)->sorted_event);
		// Remove the memory copy of this event.		
		if((*event_ptr) != 0) delete (*event_ptr);
		// Finally remove this event.
		events.erase(event_ptr);
		// Post-processing.
		post_remove_event();
	}
}

void DiscreteEventSimulator::post_remove_event() {}

void DiscreteEventSimulator::destroy(DiscreteEvent_iterator event_ptr)
{
	// Destroy all other affected events except for this event itself.
	self_destroy_flag = false;
	DiscreteEvent_iterators filtered_events;
	DiscreteEvent_iterators& destroyed_events = (*event_ptr)->destroyed_events;
	for(DiscreteEvent_iterator_iterator event_ptr_ptr = destroyed_events.begin(); event_ptr_ptr != destroyed_events.end(); ++event_ptr_ptr)
	{
		if(*event_ptr_ptr != event_ptr) filtered_events.push_back(*event_ptr_ptr);
		else self_destroy_flag = true;
	}
	remove_event(event_ptr, filtered_events);
}

// Event interactions

void DiscreteEventSimulator::add_modification(DiscreteEvent_iterator event_ptr1, DiscreteEvent_iterator event_ptr2)
{
	// Add event2 to the modified event list of event1.
	(*event_ptr1)->modified_events.push_back(event_ptr2);
	// Add event1 to the actuator event list of event2.
	(*event_ptr2)->actuator_events.push_back(event_ptr1);
}

void DiscreteEventSimulator::remove_modification(DiscreteEvent_iterator event_ptr1, DiscreteEvent_iterator event_ptr2)
{
	// Remove event2 from the modified event list of event1.
	(*event_ptr1)->modified_events.remove(event_ptr2);
	// Remove event1 from the actuator event list of event2.
	(*event_ptr2)->actuator_events.remove(event_ptr1);
}

void DiscreteEventSimulator::empty_modified_events(DiscreteEvent_iterator event_ptr)
{
	// Remove this event from the actuator event list of its modified events.
	DiscreteEvent_iterators& modified_events = (*event_ptr)->modified_events;
	for(DiscreteEvent_iterator_iterator modified_ptr_ptr = modified_events.begin(); modified_ptr_ptr != modified_events.end(); ++modified_ptr_ptr)
	{
		if(*modified_ptr_ptr != event_ptr) (**modified_ptr_ptr)->actuator_events.remove(event_ptr);
	}
	// Empty the modified event list.
	modified_events.clear();
}

void DiscreteEventSimulator::add_destruction(DiscreteEvent_iterator event_ptr1, DiscreteEvent_iterator event_ptr2)
{
	// Add event2 to the destroyed event list of event1.
	(*event_ptr1)->destroyed_events.push_back(event_ptr2);
	// Add event1 to the actuator event list of event2.
	(*event_ptr2)->actuator_events.push_back(event_ptr1);
}

void DiscreteEventSimulator::remove_destruction(DiscreteEvent_iterator event_ptr1, DiscreteEvent_iterator event_ptr2)
{
	// Remove event2 from the destroyed event list of event1.
	(*event_ptr1)->destroyed_events.remove(event_ptr2);
	// Remove event1 from the actuator event list of event2.
	(*event_ptr2)->actuator_events.remove(event_ptr1);
}

void DiscreteEventSimulator::empty_destroyed_events(DiscreteEvent_iterator event_ptr)
{
	// Remove this event from the actuator event list of its destroyed events.
	DiscreteEvent_iterators& destroyed_events = (*event_ptr)->destroyed_events;
	for(DiscreteEvent_iterator_iterator destroyed_ptr_ptr = destroyed_events.begin(); destroyed_ptr_ptr != destroyed_events.end(); ++destroyed_ptr_ptr)
	{
		if(*destroyed_ptr_ptr != event_ptr) (**destroyed_ptr_ptr)->actuator_events.remove(event_ptr);
	}
	// Empty the destroyed event list.
	destroyed_events.clear();
}

void DiscreteEventSimulator::remove_connection(DiscreteEvent_iterator event_ptr)
{
	// Remove this event from the actuator event list of its modified events.
	DiscreteEvent_iterators& modified_events = (*event_ptr)->modified_events;
	for(DiscreteEvent_iterator_iterator modified_ptr_ptr = modified_events.begin(); modified_ptr_ptr != modified_events.end(); ++modified_ptr_ptr)
	{
		if(*modified_ptr_ptr != event_ptr) (**modified_ptr_ptr)->actuator_events.remove(event_ptr);
	}
	// Empty the modified event list.
	modified_events.clear();
	// Remove this event from the actuator event list of its destroyed events.
	DiscreteEvent_iterators& destroyed_events = (*event_ptr)->destroyed_events;
	for(DiscreteEvent_iterator_iterator destroyed_ptr_ptr = destroyed_events.begin(); destroyed_ptr_ptr != destroyed_events.end(); ++destroyed_ptr_ptr)
	{
		if(*destroyed_ptr_ptr != event_ptr) (**destroyed_ptr_ptr)->actuator_events.remove(event_ptr);
	}
	// Empty the destroyed event list.
	destroyed_events.clear();
	// Remove this event from the effector event lists of its actuator events.
	DiscreteEvent_iterators& actuator_events = (*event_ptr)->actuator_events;
	for(DiscreteEvent_iterator_iterator actuator_ptr_ptr = actuator_events.begin(); actuator_ptr_ptr != actuator_events.end(); ++actuator_ptr_ptr)
	{
		if(*actuator_ptr_ptr != event_ptr)
		{
			// Remove this event from the modified event list of its actuator event.
			(**actuator_ptr_ptr)->modified_events.remove(event_ptr);
			// Remove this event from the destroyed event list of its actuator event.
			(**actuator_ptr_ptr)->destroyed_events.remove(event_ptr);
		}
	}
	// Empty the actuator event list.
	actuator_events.clear();
}

// Trigger, connect and synchronize event.

void DiscreteEventSimulator::trigger(DiscreteEvent_iterator event_ptr)
{
	(*event_ptr)->action();
}

void DiscreteEventSimulator::connect(DiscreteEvent_iterator event_ptr) {}

void DiscreteEventSimulator::synchronize(DiscreteEvent_iterator event_ptr)
{
	// Theoretically the functions of destroy, modify and create
	// are independent of each other, such that the order to call
	// them does not affect the final output.
	// Practically destroy should be called first such that the
	// events to be destroyed are not subject to modify function
	// to improve the efficiency of this algorithm. This is
	// particularly important for a situation where an event will
	// destroy itself after its action is executed.
	destroy(event_ptr);
	modify(event_ptr);
	create(event_ptr);
}

double DiscreteEventSimulator::exec(DiscreteEvent_iterator event_ptr)
{
	// Get the waiting period of the event to be executed.
	double waiting_period = (*event_ptr)->period;
	// Execute the action of this event.
	trigger(event_ptr);
	// Update the connection of this event to other events.
	connect(event_ptr);
	// Process the impact of the action of this event on both itself
	// and on other events through event interactions.
	synchronize(event_ptr);
	// Self destruction if necessary.
	if(self_destroy_flag)
	{
		remove_event(event_ptr);
		self_destroy_flag = false;
	}
	return waiting_period;
}

}
