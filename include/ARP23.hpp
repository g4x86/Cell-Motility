#ifndef ARP23_HPP_
#define ARP23_HPP_

#include <typedefs.hpp>
#include <Molecule.hpp>

namespace motility
{

class ARP23 : public Molecule
{
  public:

	ARP23(int s = 1, double diam = 0.01, const CartesianCoordinate& loc = CartesianCoordinate(), bool mob = true, const char* t = "ARP23");

	// If ARP23 is assumed to have only two states: active and
	// inactive. Then several virtual member functions in base
	// class are enough to get and set the state of ARP23. If
	// the meaning of 'state' changes later, we can overload
	// getState() and setState(), and override those virtual
	// member functions in base class.
};

inline ARP23Handle arp23_handle_null {};
inline ARP23ConstHandle arp23_const_handle_null {};

}

#endif /*ARP23_HPP_*/
