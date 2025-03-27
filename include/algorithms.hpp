//
// Principle
//
// All algorithm functions do NOT take care of the memory management of
// the input and output argument variables, which means that all the
// arguments are passed through references as much as possible. Therefore
// the calling functions should be responsible for allocating necessary
// memories for the arguments passed to algorithm functions.
//
//
// Attention
//
// Based on the reference-passing mechanism, some arguments may be changed
// after the calling of algorithm functions. A thumb rule is to make the
// input arguments as constant as possible while keep the output arguments
// changeable.
//
#ifndef ALGORITHMS_HPP_
#define ALGORITHMS_HPP_

#include <cassert>
#include <cfloat>
#include <cstdlib>
#include <cctype>
#include <cmath>
#include <string>
#include <sstream>
#include <list>
#include <algorithm>
#include <Coordinate.hpp>
#include <Vector.hpp>
#include <Line.hpp>
#include <Triangle.hpp>

// A convenient macro definition of switch-case syntax
// for char* and std::string.
#define SWITCH(str) { std::string __str_tmp = std::string(str); while(1) {
#define SWITCH_END()  break; } }
#define CASE(str) if(std::string(str) == __str_tmp)
#define DEFAULT()

namespace motility
{

// These error-handling functions are used to deal with the failure
// of opening a file to read or write.
void handleErrorEvent(const char* msg, int code = EXIT_FAILURE);

void handleErrorEvent(std::string msg, int code = EXIT_FAILURE);

void handleErrorEvent(const std::stringstream& msg, int code = EXIT_FAILURE);

double cot(double f);

// This function generates a random number uniformly distributed
// between 0 and 1.
double getUniformProbability();

bool isFractionZero(double f);

bool isEqual(double f1, double f2);

bool splitFileName(const std::string& file_name, std::string& main_name, std::string& ext_name);

// This function converts the given histogram to corresponding
// probability density distribution and return 'true' if the
// histogram is not not empty, otherwise 'false' is returned.
bool convertHistogramToProbability(const double* hist_dist, double* prob_dist, size_t size);

// This function return the index of the selected item based on
// the given probability density distribution.
size_t selectFromProbability(const double* prob_dist, size_t size);

// This function selects an element based on a given histogram
// distribution and return the status of selection.
bool selectFromHistogram(const double* hist_dist, size_t size, size_t& index);

// This function calculates the directional distance from
// a gvien point to the plane of a given triangle.
double distance(const CartesianCoordinate& point, Triangle& triangle);

// This function determines if a given point is in the place
// of a given triangle.
bool isInPlane(const CartesianCoordinate& point, Triangle& triangle);

// This function calculates the projection spot of a given
// point onto the plane of a given triangle. The projection
// point may lie outside of the given triangle.
CartesianCoordinate projection(const CartesianCoordinate& point, Triangle& triangle);

// This function calculates the homogeneous Barycentric
// coordinates of a given point with respect to a given
// triangle.
CartesianCoordinate barycentric(const CartesianCoordinate& point, Triangle& triangle);

bool isInsideTriangle(const CartesianCoordinate& point, Triangle& triangle);

// This function determines if the projection spot of a given
// point is inside a given triangle.
bool isProjectionInsideTriangle(const CartesianCoordinate& point, Triangle& triangle);

// This function uses the Barycentric coordinates to determine
// if a given line intersects with a given triangle.
bool isIntersecting(Line& line, Triangle& triangle);

// This function determins whether a given point is inside a
// given triangular box.
bool isInsideTriangularBox(const CartesianCoordinate& point, Triangle& triangle, double height);

// This function calculates the intersection point between the
// the orientation defined by the line and the plane defined by
// the triangle.
CartesianCoordinate intersection(Line& line, Triangle& triangle);

double angle(Vector& v1, Vector& v2);

double angle(Triangle& t1, Triangle& t2);

/// This function mapps a point from original coordinate system
/// to rotated coordinate system, or versus.
///
/// \param r a point in original or rotated coordinate system.
/// \param n the vector which the original coordinate system rotates
/// around.
/// \param alpha the angle by which the original coordinate system
/// rotates around n. The directions of n and alpha follow the right-
/// hand rule.
/// \param dir "true" means mapping from original coordinate system
/// to the rotated one, and "false" means reverse mapping.
/// \return the mapped point in rotated or original coordinate system.
///
/// Note: the vector n remains the same in both original and rotated
/// coordinate system.
CartesianCoordinate rotate(CartesianCoordinate r, Vector n, double alpha, bool dir = true);

int strtoi(const std::string& str);

long strtol(const std::string& str);

unsigned long strtoul(const std::string& str);

float strtof(const std::string& str);

double strtod(const std::string& str);

std::string tolower(const std::string& str);

std::string toupper(const std::string& str);

bool strtob(const std::string& str);

int uitoi(unsigned int i);

long ultol(unsigned long i);

template<typename T>
bool isContained(const std::list<T>& l, const T& x)
{
	bool contained_flag = false;
	for(typename std::list<T>::const_iterator it = l.begin(); it != l.end(); ++it)
	{
		if(x == *it)
		{
			contained_flag = true;
			break;
		}
	}
	return contained_flag;
}

template<typename T>
bool unique_append(std::list<T>& l, const T& x)
{
	bool action_flag = true;
	for(typename std::list<T>::iterator it = l.begin(); it != l.end(); ++it)
	{
		if(x == *it)
		{
			action_flag = false;
			break;
		}
	}
	if(action_flag) l.push_back(x);
	return action_flag;
}

template<typename T>
std::list<T> merge(const std::list<T>& l1, const std::list<T>& l2)
{
	std::list<T> l3;
	for(typename std::list<T>::const_iterator it1 = l1.begin(); it1 != l1.end(); ++it1)
	{
		bool redundant_flag = false;
		for(typename std::list<T>::const_iterator it3 = l3.begin(); it3 != l3.end(); ++it3)
		{
			if(*it1 == *it3)
			{
				redundant_flag = true;
				break;
			}
		}
		if(!redundant_flag) l3.push_back(*it1);
	}
	for(typename std::list<T>::const_iterator it2 = l2.begin(); it2 != l2.end(); ++it2)
	{
		bool redundant_flag = false;
		for(typename std::list<T>::const_iterator it3 = l3.begin(); it3 != l3.end(); ++it3)
		{
			if(*it2 == *it3)
			{
				redundant_flag = true;
				break;
			}
		}
		if(!redundant_flag) l3.push_back(*it2);
	}
	return l3;
}

}

#endif /*ALGORITHMS_HPP_*/
