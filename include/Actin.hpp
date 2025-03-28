#ifndef ACTIN_HPP_
#define ACTIN_HPP_

#include <typedefs.hpp>
#include <Molecule.hpp>

namespace motility
{

class Actin : public Molecule
{
	// 0: not bound,
	// 1: bound by ADF,
	int bstate;

  public:

	Actin(const char* s = "ATP", double diam = 0.0055, const CartesianCoordinate& loc = CartesianCoordinate(), bool mob = true, const char* bs = "none", const char* t = "Actin");

	const char* getState() const;

	void setState(const char* s);

	const char* getBoundState() const;

	void setBoundState(const char* s);
};

inline ActinHandle actin_sentinel {};
inline ActinConstHandle actin_const_sentinel {};

}

#endif /*ACTIN_HPP_*/
