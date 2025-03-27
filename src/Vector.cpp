#include <cfloat>
#include <cmath>
#include <Vector.hpp>

namespace motility
{

Vector::Vector(double xx, double yy, double zz)
{
	x = xx; y = yy; z = zz;
	mag = 0;
	orient.theta = 0;
	orient.phi = M_PI / 2;
	mag_flag = false;
	orient_flag = false;
}

Vector::Vector(double mg, const Orientation& ot)
{
	mag = mg; orient = ot;
	CartesianCoordinate cc(SphericalCoordinate(mag, orient.theta, orient.phi));
	x = cc.x; y = cc.y; z = cc.z;
	mag_flag = true;
	orient_flag = true;
}

Vector::Vector(const CartesianCoordinate& bp, const CartesianCoordinate& ep)
{
	x = ep.x - bp.x;
	y = ep.y - bp.y;
	z = ep.z - bp.z;
	mag = 0;
	orient.theta = 0;
	orient.phi = M_PI / 2;
	mag_flag = false;
	orient_flag = false;
}

void Vector::clearMag()
{
	mag = 0;
	mag_flag = false;
}

void Vector::clearOrient()
{
	orient.theta = 0;
	orient.phi = M_PI / 2;
	orient_flag = false;
}

double Vector::getX() const
{
	return x;
}

double Vector::getY() const
{
	return y;
}

double Vector::getZ() const
{
	return z;
}

double Vector::getMag()
{
	if(!mag_flag)
	{
		mag = std::sqrt(x * x + y * y + z * z);
		mag_flag = true;
	}
	return mag;
}

Orientation& Vector::getOrient()
{
	if(!orient_flag)
	{
		SphericalCoordinate sc(CartesianCoordinate(x, y, z));
		mag = sc.rou; orient.theta = sc.theta; orient.phi = sc.phi;
		mag_flag = true;
		orient_flag = true;
	}
	return orient;
}

void Vector::setX(double xx)
{
	x = xx;
	if(mag_flag)
	{
		mag = std::sqrt(x * x + y * y + z * z);
	}
	if(orient_flag)
	{
		SphericalCoordinate sc(CartesianCoordinate(x, y, z));
		orient.theta = sc.theta; orient.phi = sc.phi;
	}
}
	
void Vector::setY(double yy)
{
	y = yy;
	if(mag_flag)
	{
		mag = std::sqrt(x * x + y * y + z * z);
	}
	if(orient_flag)
	{
		SphericalCoordinate sc(CartesianCoordinate(x, y, z));
		orient.theta = sc.theta; orient.phi = sc.phi;
	}
}
	
void Vector::setZ(double zz)
{
	z = zz;
	if(mag_flag)
	{
		mag = std::sqrt(x * x + y * y + z * z);
	}
	if(orient_flag)
	{
		SphericalCoordinate sc(CartesianCoordinate(x, y, z));
		orient.theta = sc.theta; orient.phi = sc.phi;
	}
}

void Vector::setMag(double mg)
{
	// This function does not change vector orientation,
	// but only changes x, y and z proportionally. So it
	// only modifies 'mag_flag', not 'orient_flag'.
	double l = getMag();
	mag = mg;
	double lambda = mag / l;
	x *= lambda; y *= lambda; z *= lambda;
}
	
void Vector::setOrient(const Orientation& ot)
{
	// This function does not change vector length,
	// but only changes its orientation. And it
	// modifies both 'orient_flag' and 'mag_flag'.
	orient = ot;
	double l = getMag();
	CartesianCoordinate cc(SphericalCoordinate(l, orient.theta, orient.phi));
	x = cc.x; y = cc.y; z = cc.z;
	orient_flag = true;
}

void Vector::normalize()
{
	setMag(1);
}

Vector Vector::normal() const
{
	return Vector(1, orient);
}

bool Vector::operator==(const Vector& v) const
{
	return (std::fabs(x - v.x) < DBL_EPSILON && std::fabs(y - v.y) < DBL_EPSILON && std::fabs(z - v.z) < DBL_EPSILON);
}

bool Vector::operator!=(const Vector& v) const
{
	return (std::fabs(x - v.x) >= DBL_EPSILON || std::fabs(y - v.y) >= DBL_EPSILON || std::fabs(z - v.z) >= DBL_EPSILON);
}

const Vector& Vector::operator+() const
{
	return *this;
}

const Vector Vector::operator-() const
{
	return Vector(-x, -y, -z);
}

const Vector Vector::operator+(const Vector& v) const
{
	return Vector(x + v.x, y + v.y, z + v.z);
}

const Vector Vector::operator-(const Vector& v) const
{
	return Vector(x - v.x, y - v.y, z - v.z);
}

const Vector Vector::operator*(double r) const
{
	return Vector(x * r, y * r, z * r);
}

const Vector Vector::operator/(double r) const
{
	return Vector(x / r, y / r, z / r);
}
	
Vector& Vector::operator+=(const Vector& v)
{
	if(this != &v)
	{
		x += v.x; y += v.y; z += v.z;
		if(mag_flag)
		{
			mag = std::sqrt(x * x + y * y + z * z);
		}
		if(orient_flag)
		{
			SphericalCoordinate sc(CartesianCoordinate(x, y, z));
			orient.theta = sc.theta; orient.phi = sc.phi;
		}
	}
	return *this;
}

Vector& Vector::operator-=(const Vector& v)
{
	if(this != &v)
	{
		x -= v.x; y -= v.y; z -= v.z;
		if(mag_flag)
		{
			mag = std::sqrt(x * x + y * y + z * z);
		}
		if(orient_flag)
		{
			SphericalCoordinate sc(CartesianCoordinate(x, y, z));
			orient.theta = sc.theta; orient.phi = sc.phi;
		}
	}
	return *this;
}

Vector& Vector::operator*=(double r)
{
	x *= r; y *= r; z *= r;
	if(mag_flag)
	{
		mag = std::sqrt(x * x + y * y + z * z);
	}
	if(orient_flag)
	{
		SphericalCoordinate sc(CartesianCoordinate(x, y, z));
		orient.theta = sc.theta; orient.phi = sc.phi;
	}
	return *this;
}

Vector& Vector::operator/=(double r)
{
	x /= r; y /= r; z /= r;
	if(mag_flag)
	{
		mag = std::sqrt(x * x + y * y + z * z);
	}
	if(orient_flag)
	{
		SphericalCoordinate sc(CartesianCoordinate(x, y, z));
		orient.theta = sc.theta; orient.phi = sc.phi;
	}
	return *this;
}

std::ostream& operator<<(std::ostream& os, const Vector& v)
{
	return os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
}

double dotProd(const Vector& v1, const Vector& v2)
{
	double sum = v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
	if(std::fabs(sum) < DBL_EPSILON) sum = 0;
	return sum;
}

Vector crossProd(const Vector& v1, const Vector& v2)
{
	double x = v1.y * v2.z - v2.y * v1.z;
	double y = v1.z * v2.x - v2.z * v1.x;
	double z = v1.x * v2.y - v2.x * v1.y;
	if(std::fabs(x) < DBL_EPSILON) x = 0;
	if(std::fabs(y) < DBL_EPSILON) y = 0;
	if(std::fabs(z) < DBL_EPSILON) z = 0;
	return Vector(x, y, z);
}

double abs(Vector& v)
{
	return v.getMag();
}

void normalize(Vector& v)
{
	v.normalize();
}

CartesianCoordinate operator+(const CartesianCoordinate& p, const Vector& v)
{
	return CartesianCoordinate(p.x + v.getX(), p.y + v.getY(), p.z + v.getZ());
}

CartesianCoordinate operator+(const Vector& v, const CartesianCoordinate& p)
{
	return CartesianCoordinate(p.x + v.getX(), p.y + v.getY(), p.z + v.getZ());
}

}
