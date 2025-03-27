#include <cstdlib>
#include <cmath>
#include <DiscreteEvent.hpp>
#include <constants.hpp>

namespace simulation
{

DiscreteEvent::DiscreteEvent(double r, double dur) : Event(r, dur)
{
	period = DBL_INF_POSITIVE;
	state = false;
	modified_events.clear();
	destroyed_events.clear();
	actuator_events.clear();
}

DiscreteEvent::~DiscreteEvent() throw() {}

double DiscreteEvent::compute_period()
{
	/// Probability density function: p(t) = r * exp(-r*t)
	/// Probability function: P(t) = 1 - exp(-r*t)
	/// Sampling function: t = -ln(u) / r, where u is a random
	/// number uniformly distributed between 0 and 1.
	double u = static_cast<double>(::random()) / RAND_MAX;
	return (-std::log(u) / rate);
}

double DiscreteEvent::get_period()
{
	return period;
}

void DiscreteEvent::update_period()
{
	period = compute_period();
}

double DiscreteEvent::compute_rate()
{
	return Event::compute_rate();
}

void DiscreteEvent::update()
{
	update_rate();
	update_period();
	state = true;
}

bool DiscreteEvent::get_state()
{
	return state;
}

void DiscreteEvent::set_state(bool s)
{
	state = s;
}

DiscreteEvent& DiscreteEvent::operator=(const DiscreteEvent& de)
{
	if(&de != this)
	{
		rate = de.rate;
		duration = de.duration;
		period = de.period;
		state = de.state;
	}
	return *this;
}

bool DiscreteEvent::operator<(const DiscreteEvent& de) const
{
	return period < de.period;
}

bool DiscreteEvent::operator<=(const DiscreteEvent& de) const
{
	return period <= de.period;
}

bool DiscreteEvent::operator>(const DiscreteEvent& de) const
{
	return period > de.period;
}

bool DiscreteEvent::operator>=(const DiscreteEvent& de) const
{
	return period >= de.period;
}

bool DiscreteEvent::operator==(const DiscreteEvent& de) const
{
	return period == de.period;
}

bool DiscreteEvent::operator!=(const DiscreteEvent& de) const
{
	return period != de.period;
}

}
