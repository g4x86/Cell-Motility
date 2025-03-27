#ifndef VECTOR_HPP_
#define VECTOR_HPP_

#include <cfloat>
#include <cmath>
#include <iostream>
#include <Coordinate.hpp>

namespace motility
{

class Vector
{
	// A vector is represented by a magnitude and an
	// orientation in space.
	//
	// Note: a vector does not have a fixed location,
	// but can be located anywhere in space.

	double x, y, z;

	double mag;

	Orientation orient;
	// the default direction is posive X axis

	bool mag_flag, orient_flag;
	// Only when mag and orient are needed, they are
	// calculated.

  public:

	// Some member functions and friend functions which
	// get the value of the private data members of this
	// class, or perform some calculation, will modify
	// some data members. Therefore other classes which
	// use this class usually cannot pass 'const' reference
	// between functions.
	Vector(double xx = 0, double yy = 0, double zz = 0);

	Vector(double mg, const Orientation& ot);

	Vector(const CartesianCoordinate& bp, const CartesianCoordinate& ep);

	void clearMag();

	void clearOrient();

	double getX() const;

	double getY() const;

	double getZ() const;

	double getMag();

	Orientation& getOrient();
	// This function may modify member data

	void setX(double xx);
	
	void setY(double yy);
	
	void setZ(double zz);

	void setMag(double mg);
	
	void setOrient(const Orientation& ot);

	void normalize();

	Vector normal() const;

	bool operator==(const Vector& v) const;

	bool operator!=(const Vector& v) const;

	const Vector& operator+() const;

	const Vector operator-() const;

	const Vector operator+(const Vector& v) const;

	const Vector operator-(const Vector& v) const;

	const Vector operator*(double r) const;

	const Vector operator/(double r) const;

	Vector& operator+=(const Vector& v);

	Vector& operator-=(const Vector& v);

	Vector& operator*=(double r);

	Vector& operator/=(double r);

	friend std::ostream& operator<<(std::ostream& os, const Vector& v);

	friend double dotProd(const Vector& v1, const Vector& v2);

	friend Vector crossProd(const Vector& v1, const Vector& v2);

	friend double abs(Vector& v);

	friend void normalize(Vector& v);

	friend CartesianCoordinate operator+(const CartesianCoordinate& p, const Vector& v);

	friend CartesianCoordinate operator+(const Vector& v, const CartesianCoordinate& p);
};

}

#endif /*VECTOR_HPP_*/
