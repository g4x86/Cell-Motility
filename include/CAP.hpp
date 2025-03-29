#ifndef CAP_HPP_
#define CAP_HPP_

#include <typedefs.hpp>
#include <Molecule.hpp>

namespace motility
{

class CAP : public Molecule
{
  public:

	CAP(int s = 1, double diam = 0.01, const CartesianCoordinate& loc = CartesianCoordinate(), bool mob = true, const char* t = "CAP");

	// If CAP is assumed to have only two states: active and
	// inactive. Then several virtual member functions in base
	// class are enough to get and set the state of CAP. If
	// the meaning of 'state' changes later, we can overload
	// getState() and setState(), and override those virtual
	// member functions in base class.
};

inline CAPHandle cap_handle_null {};
inline CAPConstHandle cap_const_handle_null {};

}

#endif /*CAP_HPP_*/
