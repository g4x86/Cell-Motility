#include <algorithm>
#include <Triangle.hpp>

namespace motility
{

Triangle::Triangle()
{
	area = 0;
}

Triangle::Triangle(const CartesianCoordinate& v0, const CartesianCoordinate& v1, const CartesianCoordinate& v2)
{
	vertices[0] = v0;
	vertices[1] = v1;
	vertices[2] = v2;
	initialize();
}

Triangle::Triangle(const CartesianCoordinate* vs)
{
	vertices[0] = vs[0];
	vertices[1] = vs[1];
	vertices[2] = vs[2];
	initialize();
}

void Triangle::initialize()
{
	center = (vertices[0] + vertices[1] + vertices[2]) / 3;
	normal = crossProd(Vector(vertices[0], vertices[1]), Vector(vertices[0], vertices[2]));
	area = abs(normal) / 2;
	normalize(normal);
}

const CartesianCoordinate* Triangle::getVertices() const
{
	return vertices;
}

const CartesianCoordinate& Triangle::getCenter() const
{
	return center;
}

Vector& Triangle::getNormal()
{
	// The caller of this function may call the member function
	// of Vector, such as getMag(), getOrient() and normalize(),
	// and the caller may call the friend functions of Vector,
	// such as abs() and normalize(). All these functions may
	// modify the member data of the vector. Therefore this
	// function is not marked as 'const'.
	return normal;
}

double Triangle::getArea() const
{
	return area;
}

void Triangle::getSideLength(double* side_length) const
{
	for(size_t i = 0; i < 3; ++i)
	{
		size_t j = i;
		if(++j > 2) j = 0;
		side_length[i] = distance(vertices[i], vertices[j]);
	}
}

double Triangle::getRegularity() const
{
	//
	// Definition
	//
	// The regularity of a triangle is defined in terms of the
	// relative difference of the length of three sides:
	//
	//                     max{a(i)} - min{a(i)}
	//           1 / (1 + -----------------------)
	//                          min{a(i)}
	//
	// Therefore the value of regularity ranges (0, 1].
	//
	double side_length[3];
	getSideLength(side_length);
	double min_side_length = std::min(std::min(side_length[0], side_length[1]), side_length[2]);
	double max_side_length = std::max(std::max(side_length[0], side_length[1]), side_length[2]);
	double regularity = 1 / (1 + (max_side_length - min_side_length) / min_side_length);
	return regularity;
}

const Triangle& Triangle::operator+() const
{
	return *this;	
}

const Triangle Triangle::operator-() const
{
	return Triangle(vertices[0], vertices[2], vertices[1]);
}

bool Triangle::operator==(const Triangle& t) const
{
	return (vertices[0] == t.vertices[0] && vertices[1] == t.vertices[1] && vertices[2] == t.vertices[2]);
}

bool Triangle::operator!=(const Triangle& t) const
{
	return (vertices[0] != t.vertices[0] || vertices[1] != t.vertices[1] || vertices[2] != t.vertices[2]);
}

std::ostream& operator<<(std::ostream& os, const Triangle& t)
{
	return os << '[' << t.vertices[0] << ", " << t.vertices[1] << ", " << t.vertices[2] << ']';
}

}
