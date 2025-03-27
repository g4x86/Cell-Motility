#ifndef UNIFORMMOLECULARDISTRIBUTIONFUNCTION_HPP_
#define UNIFORMMOLECULARDISTRIBUTIONFUNCTION_HPP_

#include <MolecularDistributionFunction.hpp>

namespace motility
{

class UniformMolecularDistributionFunction : public MolecularDistributionFunction
{
  private:

	double value;

  public:

	UniformMolecularDistributionFunction(double v);

	virtual ~UniformMolecularDistributionFunction();

	double operator()(const CartesianCoordinate& loc) const;
};

}

#endif /*UNIFORMMOLECULARDISTRIBUTIONFUNCTION_HPP_*/
