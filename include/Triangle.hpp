#ifndef TRIANGLE_HPP_
#define TRIANGLE_HPP_

#include <iostream>
#include <Coordinate.hpp>
#include <Vector.hpp>

namespace motility
{

class Triangle
{
	CartesianCoordinate vertices[3];
	// The order of three vertices determins the orientation of triangle,
	// which follows the right-hand rule. The direction of triangle is
	// determined by its normal vector.

	CartesianCoordinate center;

	Vector normal;

	double area;

	void initialize();

  public:

	Triangle();

	Triangle(const CartesianCoordinate& v0, const CartesianCoordinate& v1, const CartesianCoordinate& v2);

	Triangle(const CartesianCoordinate* vs);

	const CartesianCoordinate* getVertices() const;

	const CartesianCoordinate& getCenter() const;

	Vector& getNormal();

	double getArea() const;

	void getSideLength(double* side_length) const;

	double getRegularity() const;

	const Triangle& operator+() const;

	const Triangle operator-() const;

	bool operator==(const Triangle& t) const;

	bool operator!=(const Triangle& t) const;

	friend std::ostream& operator<<(std::ostream& os, const Triangle& t);
};

}

#endif /*TRIANGLE_HPP_*/
