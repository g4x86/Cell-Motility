#ifndef COORDINATE_HPP_
#define COORDINATE_HPP_

#include <cfloat>
#include <cmath>
#include <list>
#include <iostream>

namespace motility
{

struct SphericalCoordinate;

struct CylindricalCoordinate;

struct CartesianCoordinate
{
	double x, y, z;

	CartesianCoordinate(double xx = 0, double yy = 0, double zz = 0);

	CartesianCoordinate(const SphericalCoordinate& sc);

	CartesianCoordinate(const CylindricalCoordinate& cc);

	const SphericalCoordinate toSphericalCoordinate() const;

	const CylindricalCoordinate toCylindricalCoordinate() const;

	const CartesianCoordinate& operator+() const;

	const CartesianCoordinate operator-() const;

	const CartesianCoordinate operator+(const CartesianCoordinate& cc) const;

	const CartesianCoordinate operator-(const CartesianCoordinate& cc) const;

	CartesianCoordinate& operator+=(const CartesianCoordinate& cc);

	CartesianCoordinate& operator-=(const CartesianCoordinate& cc);

	const CartesianCoordinate operator*(double c) const;

	const CartesianCoordinate operator/(double c) const;

	CartesianCoordinate& operator*=(double c);

	CartesianCoordinate& operator/=(double c);

	CartesianCoordinate& operator=(const SphericalCoordinate& sc);

	CartesianCoordinate& operator=(const CylindricalCoordinate& cc);

	bool operator==(const CartesianCoordinate& cc) const;

	bool operator!=(const CartesianCoordinate& cc) const;

	friend double distance(const CartesianCoordinate& cc1, const CartesianCoordinate& cc2);

	friend double abs(const CartesianCoordinate& cc);

	friend void normalize(CartesianCoordinate& cc);

	friend double dotProd(CartesianCoordinate& cc1, CartesianCoordinate& cc2);

	friend std::ostream& operator<<(std::ostream& os, const CartesianCoordinate& cc);
};

struct SphericalCoordinate
{
	double rou, theta, phi;

	SphericalCoordinate(double rr = 0, double tt = 0, double pp = M_PI / 2);

	SphericalCoordinate(const CartesianCoordinate& cc);

	SphericalCoordinate(const CylindricalCoordinate& cc);

	const CartesianCoordinate toCartesianCoordinate() const;

	const CylindricalCoordinate toCylindricalCoordinate() const;

	SphericalCoordinate& operator=(const CartesianCoordinate& cc);

	SphericalCoordinate& operator=(const CylindricalCoordinate& cc);

	const SphericalCoordinate& operator+() const;

	const SphericalCoordinate operator-() const;

	bool operator==(const SphericalCoordinate& sc) const;

	bool operator!=(const SphericalCoordinate& sc) const;

	friend double abs(const SphericalCoordinate& sc);

	friend void normalize(SphericalCoordinate& sc);

	friend std::ostream& operator<<(std::ostream& os, const SphericalCoordinate& sc);
};

struct CylindricalCoordinate
{
	double rou, theta, z;

	CylindricalCoordinate(double rr = 0, double tt = 0, double zz = 0);

	CylindricalCoordinate(const CartesianCoordinate& cc);

	CylindricalCoordinate(const SphericalCoordinate& sc);

	const CartesianCoordinate toCartesianCoordinate() const;

	const SphericalCoordinate toSphericalCoordinate() const;

	CylindricalCoordinate& operator=(const CartesianCoordinate& cc);

	CylindricalCoordinate& operator=(const SphericalCoordinate& sc);

	const CylindricalCoordinate& operator+() const;

	const CylindricalCoordinate operator-() const;

	bool operator==(const CylindricalCoordinate& cc) const;

	bool operator!=(const CylindricalCoordinate& cc) const;

	friend double abs(const CylindricalCoordinate& cc);

	friend void normalize(CylindricalCoordinate& sc);

	friend std::ostream& operator<<(std::ostream& os, const CylindricalCoordinate& cc);
};

struct Orientation
{
	double theta, phi;

	Orientation(double tt = 0, double pp = M_PI / 2);

	const Orientation& operator+() const;

	const Orientation operator-() const;

	bool operator==(const Orientation& orient) const;

	bool operator!=(const Orientation& orient) const;

	friend std::ostream& operator<<(std::ostream& os, const Orientation& orient);
};

struct GridCoordinate
{
	size_t i, j, k;

	GridCoordinate(size_t ii = 0, size_t jj = 0, size_t kk = 0);
};

typedef std::list<CartesianCoordinate> CartesianCoordinates;
typedef CartesianCoordinates::iterator CartesianCoordinateHandle;
typedef CartesianCoordinates::const_iterator CartesianCoordinateConstHandle;

inline CartesianCoordinateHandle cartesian_coordinate_handle_null {};
inline CartesianCoordinateConstHandle cartesian_coordinate_const_handle_null {};

typedef std::list<SphericalCoordinate> SphericalCoordinates;
typedef SphericalCoordinates::iterator SphericalCoordinateHandle;
typedef SphericalCoordinates::const_iterator SphericalCoordinateConstHandle;

inline SphericalCoordinateHandle spherical_coordinate_handle_null {};
inline SphericalCoordinateConstHandle spherical_coordinate_const_handle_null {};

typedef std::list<CylindricalCoordinate> CylindricalCoordinates;
typedef CylindricalCoordinates::iterator CylindricalCoordinateHandle;
typedef CylindricalCoordinates::const_iterator CylindricalCoordinateConstHandle;

inline CylindricalCoordinateHandle cylindrical_coordinate_handle_null {};
inline CylindricalCoordinateConstHandle cylindrical_coordinate_const_handle_null {};

typedef std::list<Orientation> Orientations;
typedef Orientations::iterator OrientationHandle;
typedef Orientations::const_iterator OrientationConstHandle;

inline OrientationHandle orientation_handle_null {};
inline OrientationConstHandle orientation_const_handle_null {};

typedef std::list<GridCoordinate> GridCoordinates;
typedef GridCoordinates::iterator GridCoordinateHandle;
typedef GridCoordinates::const_iterator GridCoordinateConstHandle;

inline GridCoordinateHandle grid_coordinate_handle_null {};
inline GridCoordinateConstHandle grid_coordinate_const_handle_null {};

}

#endif /*COORDINATE_HPP_*/
