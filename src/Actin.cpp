#include <Actin.hpp>

namespace motility
{

Actin::Actin(const char* s, double diam, const CartesianCoordinate& loc, bool mob, const char* bs, const char* t) : Molecule(0, diam, loc, mob, t)
{
	setState(s);
	setBoundState(bs);
}

const char* Actin::getState() const
{
	const char* str;
	switch(Molecule::getState())
	{
		case 0: { str = "inactive"; break; }
		case 1: { str = "ATP"; break; }
		case 2: { str = "ADP"; break; }
		case 3: { str = "ADPi"; break; }
		default: { str = ""; }
	}
	return str;
}

void Actin::setState(const char* s)
{
	std::string str(s);
	if(str == "ATP") Molecule::setState(1);
	else if(str == "ADP") Molecule::setState(2);
	else if(str == "ADPi") Molecule::setState(3);
	else if(str == "inactive") Molecule::setState(0);
	else Molecule::setState(0);
}

const char* Actin::getBoundState() const
{
	const char* str;
	switch(bstate)
	{
		case 0: { str = "none"; break; }
		case 1: { str = "ADP"; break; }
		default: { str = ""; }
	}
	return str;
}

void Actin::setBoundState(const char* s)
{
	std::string str(s);
	if(str == "none") bstate = 0;
	else if(str == "ADP") bstate = 1;
	else bstate = 0;
}

}
