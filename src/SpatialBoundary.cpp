#include <SpatialBoundary.hpp>

namespace motility
{

SpatialBoundary::SpatialBoundary(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax)
{
	x_min = xmin; x_max = xmax;
	y_min = ymin; y_max = ymax;
	z_min = zmin; z_max = zmax;
}

double SpatialBoundary::getXmin() const
{
	return x_min;
}

double SpatialBoundary::getXmax() const
{
	return x_max;
}

double SpatialBoundary::getYmin() const
{
	return y_min;
}

double SpatialBoundary::getYmax() const
{
	return y_max;
}

double SpatialBoundary::getZmin() const
{
	return z_min;
}

double SpatialBoundary::getZmax() const
{
	return z_max;
}

double SpatialBoundary::getXRange() const
{
	return (x_max - x_min);
}

double SpatialBoundary::getYRange() const
{
	return (y_max - y_min);
}

double SpatialBoundary::getZRange() const
{
	return (z_max - z_min);
}

bool SpatialBoundary::isInside(const CartesianCoordinate& loc) const
{
	bool r;
	if(x_min <= loc.x && loc.x <= x_max && y_min <= loc.y && loc.y <= y_max && z_min <= loc.z && loc.z <= z_max) r = true;
	else r = false;
	return r;
}

}
