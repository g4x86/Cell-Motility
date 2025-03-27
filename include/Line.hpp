#ifndef LINE_HPP_
#define LINE_HPP_

#include <iostream>
#include <Coordinate.hpp>
#include <Vector.hpp>

namespace motility
{

class Line
{
	/// A line is represented by a vector connecting a beginning point
	/// and an ending point.
	CartesianCoordinate begin, end;

	Vector vector;

  public:

	/// The default line starts from the origin in the direction of
	/// positive X axis and with a length of zero.
	Line();

	Line(const CartesianCoordinate& bp, const CartesianCoordinate& ep);

	Line(const CartesianCoordinate* ps);

	Line(const CartesianCoordinate& bp, const Vector& v);

	const CartesianCoordinate& getBegin() const;

	const CartesianCoordinate& getEnd() const;

	/// The input argument can be either the parameter of
	/// parameterized line equation:
	///
	/// X = X(bp) + V(bp->ep) * t
	///
	/// or the distance from any location on this line to
	/// the beginning point.
	///
	/// The behavior of this function is controlled by its
	/// second input argument 'param':
	///
	/// When 'param' is 'true', the first argument is used
	/// as the parameter of line equation.
	/// When 'param' is 'false', the first argument is used
	/// as the distance from any location on this line to
	/// the beginning point.
	CartesianCoordinate getLocation(double r, bool param = true);

	/// Since some member functions of Vector can modify its
	/// member data, this function must return a non-const
	/// reference.
	Vector& getVector();

	double length();

	Orientation& getOrient();

	void setBegin(const CartesianCoordinate& bp);

	void setEnd(const CartesianCoordinate& ep);

	const Line& operator+() const;

	const Line operator-() const;

	const Line operator+(const Line& l) const;

	bool operator==(const Line& l) const;

	bool operator!=(const Line& l) const;

	friend std::ostream& operator<<(std::ostream& os, const Line& l);
};

}

#endif /*LINE_HPP_*/
