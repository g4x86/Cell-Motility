#ifndef MOLECULARDISTRIBUTIONFUNCTION_HPP_
#define MOLECULARDISTRIBUTIONFUNCTION_HPP_

#include <functional>
#include <Coordinate.hpp>

namespace motility
{

/// MolecularDistribution class describes static molecular distribution.
///
/// NonUniformMolecularDistribution uses the instance of the derived class from
/// this class to initialize the molecular density on all grid points. 
struct MolecularDistributionFunction
{
	virtual ~MolecularDistributionFunction();

	virtual double operator()(const CartesianCoordinate& loc) const = 0;
};

}

#endif /*MOLECULARDISTRIBUTIONFUNCTION_HPP_*/
