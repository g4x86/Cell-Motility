#ifndef SPATIALBOUNDARY_HPP_
#define SPATIALBOUNDARY_HPP_

#include <Coordinate.hpp>

namespace motility
{

class SpatialBoundary
{
	double x_min, x_max;
	double y_min, y_max;
	double z_min, z_max;

  public:

	SpatialBoundary(double xmin = 0, double xmax = 0, double ymin = 0, double ymax = 0, double zmin = 0, double zmax = 0);

	double getXmin() const;

	double getXmax() const;

	double getYmin() const;

	double getYmax() const;

	double getZmin() const;

	double getZmax() const;

	double getXRange() const;

	double getYRange() const;

	double getZRange() const;

	bool isInside(const CartesianCoordinate& loc) const;
};

}

#endif /*SPATIALBOUNDARY_HPP_*/
