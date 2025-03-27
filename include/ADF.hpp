#ifndef ADF_HPP_
#define ADF_HPP_

#include <Molecule.hpp>

namespace motility
{

class ADF : public Molecule
{
  public:

	ADF(int s = 1, double diam = 0.0055, const CartesianCoordinate& loc = CartesianCoordinate(), bool mob = true, const char* t = "ADF");

	// If CAP is assumed to have only two states: active and
	// inactive. Then several virtual member functions in base
	// class are enough to get and set the state of CAP. If
	// the meaning of 'state' changes later, we can overload
	// getState() and setState(), and override those virtual
	// member functions in base class.
};

}

#endif /*ADF_HPP_*/
