#ifndef MOLECULE_HPP_
#define MOLECULE_HPP_

#include <cmath>
#include <iostream>
#include <string>
#include <Coordinate.hpp>

namespace motility
{

class Molecule
{
	std::string type;

	CartesianCoordinate location;

	// 0: inactive
	// The meaning of other values of 'state' is determined
	// by the derived classes by defining two virtual member
	// functions activate() and deactivate().
	int state;

	// Unit: micrometer.
	double diameter;

	bool mobility;

  protected:

	int getState() const;

	void setState(int s);

	void setMobility(bool mob);

  public:

	Molecule(int s = 0, double diam = 0.01, const CartesianCoordinate& loc = CartesianCoordinate(), bool mob = true, const char* t = "generic Molecule");

	virtual ~Molecule();

	std::string getType() const;

	void setType(const char* t);

	double getDiameter() const;

	const CartesianCoordinate& getLocation() const;

	void setLocation(const CartesianCoordinate& loc);

	void translocate(const CartesianCoordinate& dLoc);

	bool isMobile() const;

	void mobilize();

	void immobilize();

	virtual double getVolume() const;

	virtual bool isActive() const;

	// Because activate() and deactivate() operate on 'state'
	// and the meaning of 'state' is interpreted in the derived
	// classes, these two functions should NOT be used directly.
	// Instead, the derived classes will overload these functions
	// to provide meaningful operations while interpreting 'state'.

	virtual void deactivate();

	virtual void activate();

	friend std::ostream& operator<<(std::ostream& os, const Molecule& m);
};

}

#endif /*MOLECULE_HPP_*/
