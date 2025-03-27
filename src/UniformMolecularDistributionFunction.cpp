#include <UniformMolecularDistributionFunction.hpp>

namespace motility
{

UniformMolecularDistributionFunction::UniformMolecularDistributionFunction(double v)
{
	value = v;
}

UniformMolecularDistributionFunction::~UniformMolecularDistributionFunction() {}

double UniformMolecularDistributionFunction::operator()(const CartesianCoordinate& loc) const
{
	return value;
}

}
