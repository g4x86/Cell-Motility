#include <Molecule.hpp>

namespace motility
{

Molecule::Molecule(int s, double diam, const CartesianCoordinate& loc, bool mob, const char* t)
{
	type = t;
	state = s;
	diameter = diam;
	location = loc;
	mobility = mob;
}

Molecule::~Molecule() {}

int Molecule::getState() const
{
	return state;
}

void Molecule::setState(int s)
{
	state = s;
}

void Molecule::setMobility(bool mob)
{
	mobility = mob;
}

std::string Molecule::getType() const
{
	return type;
}

void Molecule::setType(const char* t)
{
	type = t;
}

double Molecule::getDiameter() const
{
	return diameter;
}

const CartesianCoordinate& Molecule::getLocation() const
{
	return location;
}

void Molecule::setLocation(const CartesianCoordinate& loc)
{
	location = loc;
}

void Molecule::translocate(const CartesianCoordinate& dLoc)
{
	setLocation(location + dLoc);
}

bool Molecule::isMobile() const
{
	return mobility;
}

void Molecule::mobilize()
{
	setMobility(true);
}

void Molecule::immobilize()
{
	setMobility(false);
}

double Molecule::getVolume() const
{
	return M_PI * diameter * diameter * diameter / 6;
}

bool Molecule::isActive() const
{
	return (getState() == 0) ? false : true;
}

void Molecule::deactivate()
{
	setState(0);
}

void Molecule::activate()
{
	setState(1);
}

std::ostream& operator<<(std::ostream& os, const Molecule& m)
{
	const char* str1;
	if(m.state) str1 = "Active";
	else str1 = "Inactive";
	const char* str2;
	if(m.mobility) str2 = "diffusive";
	else str2 = "not diffusive";
	return os << str1 << " " << m.getType() << " : " << m.diameter << "um at " << m.location << ", " << str2 << ".";
}

}
