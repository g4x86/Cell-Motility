#include <UniformMolecularDistribution.hpp>

namespace motility
{

UniformMolecularDistribution::UniformMolecularDistribution(const SpatialBoundary& bound, double dens) : MolecularDistribution(bound)
{
	density = dens;
}

double UniformMolecularDistribution::getDensity(const CartesianCoordinate& loc)
{
	double dens;
	if(boundary.isInside(loc)) dens = density;
	else dens = 0;
	return dens;
}

void UniformMolecularDistribution::setDensity(double dens)
{
	density = dens;
}

}
