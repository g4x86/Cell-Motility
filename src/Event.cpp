#include <Event.hpp>

namespace simulation
{

Event::Event(double r, double dur)
{
	rate = r;
	duration = dur;
}

Event::~Event() throw() {}

double Event::compute_rate()
{
	return rate;
}

double Event::get_rate()
{
	return rate;
}

void Event::set_rate(double r)
{
	rate = r;
}

void Event::update_rate()
{
	rate = compute_rate();
}

}
