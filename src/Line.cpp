#include <Line.hpp>

namespace motility
{

Line::Line() {}

Line::Line(const CartesianCoordinate& bp, const CartesianCoordinate& ep)
{
	begin = bp;
	end = ep;
	vector = Vector(begin, end);
}

Line::Line(const CartesianCoordinate* ps)
{
	begin = ps[0];
	end = ps[1];
	vector = Vector(begin, end);
}

Line::Line(const CartesianCoordinate& bp, const Vector& v)
{
	begin = bp;
	vector = v;
	end = bp + v;
}

const CartesianCoordinate& Line::getBegin() const
{
	return begin;
}

const CartesianCoordinate& Line::getEnd() const
{
	return end;
}

CartesianCoordinate Line::getLocation(double r, bool param)
{
	double t;
	if(param) t = r;
	else t = r / vector.getMag();
	return CartesianCoordinate(begin + vector * t);
}

Vector& Line::getVector()
{
	return vector;
}

double Line::length()
{
	return vector.getMag();
}

Orientation& Line::getOrient()
{
	return vector.getOrient();
}

void Line::setBegin(const CartesianCoordinate& bp)
{
	begin = bp;
	vector = Vector(begin, end);
}

void Line::setEnd(const CartesianCoordinate& ep)
{
	end = ep;
	vector = Vector(begin, end);
}

const Line& Line::operator+() const
{
	return *this;
}

const Line Line::operator-() const
{
	return Line(end, begin);
}

const Line Line::operator+(const Line& l) const
{
	Line ll;
	if(begin == l.begin)
	{
		ll.begin = begin;
		ll.end = end;
		ll.vector = vector + l.vector;
	}
	return ll;
}

bool Line::operator==(const Line& l) const
{
	return (begin == l.begin && end == l.end);
}

bool Line::operator!=(const Line& l) const
{
	return (begin != l.begin || end != l.end);
}

std::ostream& operator<<(std::ostream& os, const Line& l)
{
	return os << "[" << l.begin << ", " << l.end << "]";
}

}
