#ifndef UNIFORMMOLECULARDISTRIBUTION_HPP_
#define UNIFORMMOLECULARDISTRIBUTION_HPP_

#include <MolecularDistribution.hpp>

namespace motility
{

/// UniformMolecularDistribution class describes uniform molecular distribution.
///
/// Uniform molecular distribution is described by a single qunatity, usually
/// molecular density, throught out the spatial region between boundaries.
class UniformMolecularDistribution : public MolecularDistribution 
{
  private:

	double density;

  public:

	/// UniformMolecularDistribution constructor function.
	///
	/// \param bound the spatial boundary.
	/// \param dens the molecular density.
	UniformMolecularDistribution(const SpatialBoundary& bound, double dens);

	double getDensity(const CartesianCoordinate& loc);

	void setDensity(double dens);
};

}

#endif /*UNIFORMMOLECULARDISTRIBUTION_HPP_*/
