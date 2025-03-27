#include <cfloat>
#include <cmath>
#include <EventSimulator.hpp>
#include <constants.hpp>

namespace simulation
{

EventSimulator::EventSimulator(double max_duration, double time_step, double record_time_interval, size_t record_step_interval)
{
	status = 0;
	max_time_duration = max_duration;
	time_slice = time_step;
	double max_loop_step_double = max_time_duration / time_slice;
	double max_loop_step_int = std::floor(max_loop_step_double);
	max_loop_step = static_cast<size_t>(max_loop_step_int);
	if(max_loop_step_double > max_loop_step_int) ++max_loop_step;
	record_time_moment_interval = record_time_interval;
	record_loop_step_interval = record_step_interval;
	time_moment = 0;
	last_record_time_moment = time_moment;
	loop_step = 0;
	last_record_loop_step = loop_step;
}

EventSimulator::EventSimulator(double max_duration, size_t max_step, double record_time_interval, size_t record_step_interval)
{
	status = 0;
	max_time_duration = max_duration;
	max_loop_step = max_step;
	record_time_moment_interval = record_time_interval;
	record_loop_step_interval = record_step_interval;
	time_slice = 0;
	time_moment = 0;
	last_record_time_moment = time_moment;
	loop_step = 0;
	last_record_loop_step = loop_step;
}

EventSimulator::~EventSimulator() throw() {}

void EventSimulator::record()
{
	if(time_moment - last_record_time_moment >= record_time_moment_interval)
	{
		time_record();
		last_record_time_moment = time_moment;
	}
	if(loop_step - last_record_loop_step >= record_loop_step_interval)
	{
		step_record();
		last_record_loop_step = loop_step;
	}
}

void EventSimulator::time_record() {}

void EventSimulator::step_record() {}

void EventSimulator::initialize() {}

void EventSimulator::finalize()
{
	time_record();
	last_record_time_moment = time_moment;
	last_record_loop_step = loop_step;
}

void EventSimulator::run()
{
	initialize();
	for(loop_step = 1; loop_step <= max_loop_step; ++loop_step)
	{
		time_slice = update();
		if(time_slice >= DBL_INF_POSITIVE)
		{
			status = 2;
			break;
		}
		time_moment += time_slice;
		if(time_moment >= max_time_duration)
		{
			status = 1;
			break;
		}
		record();
	}
	finalize();
}

int EventSimulator::get_status() const throw()
{
	return status;
}

}
