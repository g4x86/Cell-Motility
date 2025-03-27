#ifndef MOLECULARDISTRIBUTION_HPP_
#define MOLECULARDISTRIBUTION_HPP_

#include <Coordinate.hpp>
#include <SpatialBoundary.hpp>

namespace motility
{

/// MolecularDistribution class defines the interface of molecular distribution.
///
/// This class describes the abstract base of actual molecular distributions
/// in space. The spatial distribution of molecules can be uniform, non-uniform,
/// discrete, continuous, or combination of them.
class MolecularDistribution
{
  protected:

	/// The boundary of molecular distribution.
	SpatialBoundary boundary;

  protected:

	/// Protected MolecularDistribution constructor function.
	///
	/// \param bound the spatial boundary.
	MolecularDistribution(const SpatialBoundary& bound);

	/// Protected MolecularDistribution destructor function.
	virtual ~MolecularDistribution();

  public:

	/// This function returns the molecular density at specified location.
	///
	/// \param loc the spatial location.
	/// \return molecular density.
	virtual double getDensity(const CartesianCoordinate& loc) = 0;
};

}

#endif /*MOLECULARDISTRIBUTION_HPP_*/
