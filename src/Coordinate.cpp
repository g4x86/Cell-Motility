#include <cfloat>
#include <cmath>
#include <Coordinate.hpp>

namespace motility
{

// CartesianCoordinate definition

CartesianCoordinate::CartesianCoordinate(double xx, double yy, double zz)
{
	x = xx; y = yy; z = zz;
}

CartesianCoordinate::CartesianCoordinate(const SphericalCoordinate& sc)
{
	CartesianCoordinate cc = sc.toCartesianCoordinate();
	x = cc.x; y = cc.y; z = cc.z;
}

CartesianCoordinate::CartesianCoordinate(const CylindricalCoordinate& cc)
{
	CartesianCoordinate cc0 = cc.toCartesianCoordinate();
	x = cc0.x; y = cc0.y; z = cc0.z;
}

const SphericalCoordinate CartesianCoordinate::toSphericalCoordinate() const
{
	double rou = std::sqrt(x * x + y * y + z * z);
	double theta;
	double phi = std::acos(z / rou);
	if(std::fabs(phi) < DBL_EPSILON)
	{
		phi = 0;
		theta = 0;
	}
	else if(std::fabs(phi - M_PI) < DBL_EPSILON)
	{
		phi = M_PI;
		theta = 0;
	}
	else
	{
		if(std::fabs(x) < DBL_EPSILON)
		{
			if(std::fabs(y) < DBL_EPSILON) theta = 0;
			else if(y > 0) theta = 0.5 * M_PI;
			else theta = 1.5 * M_PI;
		}
		else
		{
			if(std::fabs(y) < DBL_EPSILON)
			{
				if(x > 0) theta = 0;
				else theta = M_PI;
			}
			else
			{
				theta = std::atan(y / x);
				if(std::fabs(theta) < DBL_EPSILON) theta = 0;
				if((x < 0 && y > 0) || (x < 0 && y < 0)) theta += M_PI;
				else if(x > 0 && y < 0) theta += 2 * M_PI;
			}
		}
	}
	return SphericalCoordinate(rou, theta, phi);
}

const CylindricalCoordinate CartesianCoordinate::toCylindricalCoordinate() const
{
	double rou = std::sqrt(x * x + y * y);
	double theta;
	if(std::fabs(x) < DBL_EPSILON)
	{
		if(std::fabs(y) < DBL_EPSILON) theta = 0;
		else if(y > 0) theta = 0.5 * M_PI;
		else theta = 1.5 * M_PI;
	}
	else
	{
		if(std::fabs(y) < DBL_EPSILON)
		{
			if(x > 0) theta = 0;
			else theta = M_PI;
		}
		else
		{
			theta = std::atan(y / x);
			if(std::fabs(theta) < DBL_EPSILON) theta = 0;
			if((x < 0 && y > 0) || (x < 0 && y < 0)) theta += M_PI;
			else if(x > 0 && y < 0) theta += 2 * M_PI;
		}
	}
	return CylindricalCoordinate(rou, theta, z);
}

const CartesianCoordinate& CartesianCoordinate::operator+() const
{
	return *this;
}

const CartesianCoordinate CartesianCoordinate::operator-() const
{
	return CartesianCoordinate(-x, -y, -z);
}

const CartesianCoordinate CartesianCoordinate::operator+(const CartesianCoordinate& cc) const
{
	double xx = x + cc.x;
	double yy = y + cc.y;
	double zz = z + cc.z;
	if(std::fabs(xx) < DBL_EPSILON) xx = 0;
	if(std::fabs(yy) < DBL_EPSILON) yy = 0;
	if(std::fabs(zz) < DBL_EPSILON) zz = 0;
	return CartesianCoordinate(xx, yy, zz);
}

const CartesianCoordinate CartesianCoordinate::operator-(const CartesianCoordinate& cc) const
{
	double xx = x - cc.x;
	double yy = y - cc.y;
	double zz = z - cc.z;
	if(std::fabs(xx) < DBL_EPSILON) xx = 0;
	if(std::fabs(yy) < DBL_EPSILON) yy = 0;
	if(std::fabs(zz) < DBL_EPSILON) zz = 0;
	return CartesianCoordinate(xx, yy, zz);
}

CartesianCoordinate& CartesianCoordinate::operator+=(const CartesianCoordinate& cc)
{
	if(this != &cc)
	{
		x += cc.x; y += cc.y; z += cc.z;
		if(std::fabs(x) < DBL_EPSILON) x = 0;
		if(std::fabs(y) < DBL_EPSILON) y = 0;
		if(std::fabs(z) < DBL_EPSILON) z = 0;
	}
	return *this;
}

CartesianCoordinate& CartesianCoordinate::operator-=(const CartesianCoordinate& cc)
{
	if(this != &cc)
	{
		x -= cc.x; y -= cc.y; z -= cc.z;
		if(std::fabs(x) < DBL_EPSILON) x = 0;
		if(std::fabs(y) < DBL_EPSILON) y = 0;
		if(std::fabs(z) < DBL_EPSILON) z = 0;
	}
	return *this;
}

const CartesianCoordinate CartesianCoordinate::operator*(double c) const
{
	double xx = x * c;
	double yy = y * c;
	double zz = z * c;
	if(std::fabs(xx) < DBL_EPSILON) xx = 0;
	if(std::fabs(yy) < DBL_EPSILON) yy = 0;
	if(std::fabs(zz) < DBL_EPSILON) zz = 0;
	return CartesianCoordinate(xx, yy, zz);
}

const CartesianCoordinate CartesianCoordinate::operator/(double c) const
{
	double xx = x / c;
	double yy = y / c;
	double zz = z / c;
	if(std::fabs(xx) < DBL_EPSILON) xx = 0;
	if(std::fabs(yy) < DBL_EPSILON) yy = 0;
	if(std::fabs(zz) < DBL_EPSILON) zz = 0;
	return CartesianCoordinate(xx, yy, zz);
}

CartesianCoordinate& CartesianCoordinate::operator*=(double c)
{
	x *= c; y *= c; z *= c;
	if(std::fabs(x) < DBL_EPSILON) x = 0;
	if(std::fabs(y) < DBL_EPSILON) y = 0;
	if(std::fabs(z) < DBL_EPSILON) z = 0;
	return *this;
}

CartesianCoordinate& CartesianCoordinate::operator/=(double c)
{
	x /= c; y /= c; z /= c;
	if(std::fabs(x) < DBL_EPSILON) x = 0;
	if(std::fabs(y) < DBL_EPSILON) y = 0;
	if(std::fabs(z) < DBL_EPSILON) z = 0;
	return *this;
}

CartesianCoordinate& CartesianCoordinate::operator=(const SphericalCoordinate& sc)
{
	*this = sc.toCartesianCoordinate();
	return *this;
}

CartesianCoordinate& CartesianCoordinate::operator=(const CylindricalCoordinate& cc)
{
	*this = cc.toCartesianCoordinate();
	return *this;
}

bool CartesianCoordinate::operator==(const CartesianCoordinate& cc) const
{
	return (std::fabs(x - cc.x) < DBL_EPSILON && std::fabs(y - cc.y) < DBL_EPSILON && std::fabs(z - cc.z) < DBL_EPSILON);
}

bool CartesianCoordinate::operator!=(const CartesianCoordinate& cc) const
{
	return (std::fabs(x - cc.x) >= DBL_EPSILON || std::fabs(y - cc.y) >= DBL_EPSILON || std::fabs(z - cc.z) >= DBL_EPSILON);
}

double distance(const CartesianCoordinate& cc1, const CartesianCoordinate& cc2)
{
	double dx = cc1.x - cc2.x;
	double dy = cc1.y - cc2.y;
	double dz = cc1.z - cc2.z;
	return std::sqrt(dx * dx + dy * dy + dz * dz);
}

double abs(const CartesianCoordinate& cc)
{
	return std::sqrt(cc.x * cc.x + cc.y * cc.y + cc.z * cc.z);
}

void normalize(CartesianCoordinate& cc)
{
	cc /= abs(cc);
}

double dotProd(CartesianCoordinate& cc1, CartesianCoordinate& cc2)
{
	double sum = cc1.x * cc2.x + cc1.y * cc2.y + cc1.z * cc2.z;
	if(std::fabs(sum) < DBL_EPSILON) sum = 0;
	return sum;
}

std::ostream& operator<<(std::ostream& os, const CartesianCoordinate& cc)
{
	return os << "(" << cc.x << ", " << cc.y << ", " << cc.z << ")";
}


// SphericalCoordinate definition

SphericalCoordinate::SphericalCoordinate(double rr, double tt, double pp)
{
	rou = rr; theta = tt; phi = pp;
}

SphericalCoordinate::SphericalCoordinate(const CartesianCoordinate& cc)
{
	SphericalCoordinate sc = cc.toSphericalCoordinate();
	rou = sc.rou; theta = sc.theta; phi = sc.phi;
}

SphericalCoordinate::SphericalCoordinate(const CylindricalCoordinate& cc)
{
	SphericalCoordinate sc = cc.toSphericalCoordinate();
	rou = sc.rou; theta = sc.theta; phi = sc.phi;
}

const CartesianCoordinate SphericalCoordinate::toCartesianCoordinate() const
{
	double x = rou * std::cos(theta) * std::sin(phi);
	double y = rou * std::sin(theta) * std::sin(phi);
	double z = rou * std::cos(phi);
	if(std::fabs(x) < DBL_EPSILON) x = 0;
	if(std::fabs(y) < DBL_EPSILON) y = 0;
	if(std::fabs(z) < DBL_EPSILON) z = 0;
	return CartesianCoordinate(x, y, z);
}

const CylindricalCoordinate SphericalCoordinate::toCylindricalCoordinate() const
{
	double rr = rou * std::sin(phi);
	double z = rou * std::cos(phi);
	if(std::fabs(rr) < DBL_EPSILON) rr = 0;
	if(std::fabs(z) < DBL_EPSILON) z = 0;
	return CylindricalCoordinate(rr, theta, z);
}

SphericalCoordinate& SphericalCoordinate::operator=(const CartesianCoordinate& cc)
{
	*this = cc.toSphericalCoordinate();
	return *this;
}

SphericalCoordinate& SphericalCoordinate::operator=(const CylindricalCoordinate& cc)
{
	*this = cc.toSphericalCoordinate();
	return *this;
}

const SphericalCoordinate& SphericalCoordinate::operator+() const
{
	return *this;
}

const SphericalCoordinate SphericalCoordinate::operator-() const
{
	double rr = rou;
	double tt = theta + M_PI;
	if(tt > 2 * M_PI) tt -= 2 * M_PI;
	double pp = M_PI - phi;
	return SphericalCoordinate(rr, tt, pp);
}

bool SphericalCoordinate::operator==(const SphericalCoordinate& sc) const
{
	return (std::fabs(rou - sc.rou) < DBL_EPSILON && std::fabs(theta - sc.theta) < DBL_EPSILON && std::fabs(phi - sc.phi) < DBL_EPSILON);
}

bool SphericalCoordinate::operator!=(const SphericalCoordinate& sc) const
{
	return (std::fabs(rou - sc.rou) >= DBL_EPSILON || std::fabs(theta - sc.theta) >= DBL_EPSILON || std::fabs(phi - sc.phi) >= DBL_EPSILON);
}

double abs(const SphericalCoordinate& sc)
{
	return sc.rou;
}

void normalize(SphericalCoordinate& sc)
{
	sc.rou = 1;
}

std::ostream& operator<<(std::ostream& os, const SphericalCoordinate& sc)
{
	return os << "(" << sc.rou << ", " << sc.theta << ", " << sc.phi << ")";
}


// CylindricalCoordinate definition

CylindricalCoordinate::CylindricalCoordinate(double rr, double tt, double zz)
{
	rou = rr; theta = tt; z = zz;
}

CylindricalCoordinate::CylindricalCoordinate(const CartesianCoordinate& cc)
{
	CylindricalCoordinate cc0 = cc.toCylindricalCoordinate();
	rou = cc0.rou; theta = cc0.theta; z = cc0.z;
}

CylindricalCoordinate::CylindricalCoordinate(const SphericalCoordinate& sc)
{
	CylindricalCoordinate cc = sc.toCylindricalCoordinate();
	rou = cc.rou; theta = cc.theta; z = cc.z;
}

const CartesianCoordinate CylindricalCoordinate::toCartesianCoordinate() const
{
	double x = rou * std::cos(theta);
	double y = rou * std::sin(theta);
	if(std::fabs(x) < DBL_EPSILON) x = 0;
	if(std::fabs(y) < DBL_EPSILON) y = 0;
	return CartesianCoordinate(x, y, z);
}

const SphericalCoordinate CylindricalCoordinate::toSphericalCoordinate() const
{
	double rr = std::sqrt(rou * rou + z * z);
	double tt;
	double phi = std::acos(z / rr);
	if(std::fabs(phi) < DBL_EPSILON)
	{
		phi = 0;
		tt = 0;
	}
	else if(std::fabs(phi - M_PI) < DBL_EPSILON)
	{
		phi = M_PI;
		tt = 0;
	}
	else tt = theta;
	return SphericalCoordinate(rr, tt, phi);
}

CylindricalCoordinate& CylindricalCoordinate::operator=(const CartesianCoordinate& cc)
{
	*this = cc.toCylindricalCoordinate();
	return *this;
}

CylindricalCoordinate& CylindricalCoordinate::operator=(const SphericalCoordinate& sc)
{
	*this = sc.toCylindricalCoordinate();
	return *this;
}

const CylindricalCoordinate& CylindricalCoordinate::operator+() const
{
	return *this;
}

const CylindricalCoordinate CylindricalCoordinate::operator-() const
{
	double rr = rou;
	double tt = theta + M_PI;
	if(tt > 2 * M_PI) tt -= 2 * M_PI;
	double zz = -z;
	return CylindricalCoordinate(rr, tt, zz);
}

bool CylindricalCoordinate::operator==(const CylindricalCoordinate& cc) const
{
	return (std::fabs(rou - cc.rou) < DBL_EPSILON && std::fabs(theta - cc.theta) < DBL_EPSILON && std::fabs(z - cc.z) < DBL_EPSILON);
}

bool CylindricalCoordinate::operator!=(const CylindricalCoordinate& cc) const
{
	return (std::fabs(rou - cc.rou) >= DBL_EPSILON || std::fabs(theta - cc.theta) >= DBL_EPSILON || std::fabs(z - cc.z) < DBL_EPSILON);
}

double abs(const CylindricalCoordinate& cc)
{
	return std::sqrt(cc.rou * cc.rou + cc.z * cc.z);
}

void normalize(CylindricalCoordinate& cc)
{
	double k = cc.z / cc.rou;
	cc.rou = 1 / std::sqrt(1 + k * k);
	cc.z = cc.rou * k;
}

std::ostream& operator<<(std::ostream& os, const CylindricalCoordinate& cc)
{
	return os << "(" << cc.rou << ", " << cc.theta << ", " << cc.z << ")";
}


// Orientation definition

Orientation::Orientation(double tt, double pp)
{
	theta = tt;
	if(theta < 0) theta += 2 * M_PI;
	else if(theta > 2 * M_PI) theta -= 2 * M_PI;
	phi = pp;
}

const Orientation& Orientation::operator+() const
{
	return *this;
}

const Orientation Orientation::operator-() const
{
	double tt = theta + M_PI;
	if(tt > 2 * M_PI) tt -= 2 * M_PI;
	double pp = M_PI - phi;
	return Orientation(tt, pp);
}

bool Orientation::operator==(const Orientation& orient) const
{
	return (std::fabs(theta - orient.theta) < DBL_EPSILON && std::fabs(phi - orient.phi) < DBL_EPSILON);
}

bool Orientation::operator!=(const Orientation& orient) const
{
	return (std::fabs(theta - orient.theta) >= DBL_EPSILON || std::fabs(phi - orient.phi) >= DBL_EPSILON);
}

std::ostream& operator<<(std::ostream& os, const Orientation& orient)
{
	return os << "(" << orient.theta << ", " << orient.phi << ")";
}


// Orientation definition

GridCoordinate::GridCoordinate(size_t ii, size_t jj, size_t kk)
{
	i = ii;
	j = jj;
	k = kk;
}

}
