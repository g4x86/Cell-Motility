#include <algorithms.hpp>
#include <ParameterTable.hpp>

namespace motility
{

void handleErrorEvent(const char* msg, int code)
{
	std::cerr << "ERROR: " << msg << "!" << std::endl;
	std::exit(code);
}

void handleErrorEvent(std::string msg, int code)
{
	std::cerr << "ERROR: " << msg << "!" << std::endl;
	std::exit(code);
}

void handleErrorEvent(const std::stringstream& msg, int code)
{
	std::cerr << "ERROR: " << msg.str() << "!" << std::endl;
	std::exit(code);
}

double cot(double f)
{
	return (1 / std::tan(f));
}

double getUniformProbability()
{
	return static_cast<double>(::random()) / RAND_MAX;
}

bool isFractionZero(double f)
{
	bool r;
	if(std::fabs(f - static_cast<long>(f)) < DBL_EPSILON) r = true;
	else r = false;
	return r;
}

bool isEqual(double f1, double f2)
{
	bool r;
	if(std::fabs(f1 - f2) < DBL_EPSILON) r = true;
	else r = false;
	return r;
}

bool splitFileName(const std::string& file_name, std::string& main_name, std::string& ext_name)
{
	bool action;
	size_t len = file_name.length();
	size_t pos = file_name.rfind('.', len);
	if(pos < len - 1)
	{
		size_t n = len - pos - 1;
		main_name = file_name.substr(0, pos);
		ext_name = file_name.substr(pos + 1, n);
		action = true;
	}
	else action = false;
	return action;
}

bool convertHistogramToProbability(const double* hist_dist, double* prob_dist, size_t size)
{
	bool action;
	double sum = 0;
	size_t idx;
	for(idx = 0; idx < size; ++idx) sum += hist_dist[idx];
	if(sum > 0)
	{
		action = true;
		for(idx = 0; idx < size; ++idx) prob_dist[idx] = hist_dist[idx] / sum;
		for(idx = 0; idx < size; ++idx) if(isEqual(prob_dist[idx], 0)) prob_dist[idx] = 0;
	}
	else action = false;
	return action;
}

size_t selectFromProbability(const double* prob_dist, size_t size)
{
	double* cumul_prob_dist = 0;
	cumul_prob_dist = new double [size];
	size_t idx;
	for(idx = 0; idx < size; ++idx)
	{
		if(idx == 0) cumul_prob_dist[idx] = prob_dist[idx];
		else cumul_prob_dist[idx] = cumul_prob_dist[idx - 1] + prob_dist[idx];
	}
	size_t selected_idx = 0;
	double u = getUniformProbability();
	for(idx = 0; idx < size; ++idx)
	{
		if(u <= cumul_prob_dist[idx])
		{
			selected_idx = idx;
			break;
		}
		else if(u <= cumul_prob_dist[idx + 1])
		{
			selected_idx = idx + 1;
			break;
		}
		else {}
	}
	if(cumul_prob_dist != 0) delete [] cumul_prob_dist;
	return selected_idx;
}

bool selectFromHistogram(const double* hist_dist, size_t size, size_t& index)
{
	bool select_flag;
	double* prob_dist = 0;
	prob_dist = new double [size];
	if(convertHistogramToProbability(hist_dist, prob_dist, size))
	{
		index = selectFromProbability(prob_dist, size);
		select_flag = true;
	}
	else select_flag = false;
	if(prob_dist != 0) delete [] prob_dist;
	return select_flag;
}

double distance(const CartesianCoordinate& point, Triangle& triangle)
{
	return dotProd(triangle.getNormal(), Vector((triangle.getVertices())[0], point));
}

bool isInPlane(const CartesianCoordinate& point, Triangle& triangle)
{
	return isEqual(distance(point, triangle), 0);
}

CartesianCoordinate projection(const CartesianCoordinate& point, Triangle& triangle)
{
	CartesianCoordinate p;
	double dist = distance(point, triangle);
	if(isEqual(dist, 0)) p = point;
	else
	{
		Vector v = triangle.getNormal() * (-dist);
		Line l(point, v);
		p = l.getEnd();
	}
	return p;
}

CartesianCoordinate barycentric(const CartesianCoordinate& point, Triangle& triangle)
{
	const CartesianCoordinate* vertices = triangle.getVertices();
	Triangle t0(point, vertices[1], vertices[2]);
	Triangle t1(point, vertices[2], vertices[0]);
	Triangle t2(point, vertices[0], vertices[1]);
	double area = triangle.getArea();
	double x = t0.getArea() / area;
	double y = t1.getArea() / area;
	double z = t2.getArea() / area;
	return CartesianCoordinate(x, y, z);
}

bool isInsideTriangle(const CartesianCoordinate& point, Triangle& triangle)
{
	//
	// This function uses the Barycentric coordinates to determine if
	// a point in the plane of a given triangle is located within the
	// triangle. Given the Barycentric coordinates of a point
	// P (x, y, z), if P is located in the interior of the triangle,
	// then the summation of x, y and z is 1, otherwise it is greater
	// than 1.
	//
	CartesianCoordinate p = barycentric(point, triangle);
	return isEqual(p.x + p.y + p.z, 1);
}

bool isProjectionInsideTriangle(const CartesianCoordinate& point, Triangle& triangle)
{
	return isInsideTriangle(projection(point, triangle), triangle);
}

bool isIntersecting(Line& line, Triangle& triangle)
{
	bool intersecting;
	Vector n = triangle.getNormal();
	Vector l = line.getVector();
	const CartesianCoordinate* vertices = triangle.getVertices();
	Vector v = Vector(line.getBegin(), vertices[0]);
	double r = dotProd(n, v) / dotProd(n, l);
	if(isEqual(r, 0)) r = 0;
	if(isEqual(r, 1)) r = 1;
	if(r < 0 || r > 1) intersecting = false;
	else intersecting = isInsideTriangle(line.getLocation(r), triangle);
	return intersecting;
}

CartesianCoordinate intersection(Line& line, Triangle& triangle)
{
	Vector n = triangle.getNormal();
	Vector l = line.getVector();
	const CartesianCoordinate* vertices = triangle.getVertices();
	Vector v = Vector(line.getBegin(), vertices[0]);
	double r = dotProd(n, v) / dotProd(n, l);
	if(isEqual(r, 0)) r = 0;
	if(isEqual(r, 1)) r = 1;
	return line.getLocation(r);
}

bool isInsideTriangularBox(const CartesianCoordinate& point, Triangle& triangle, double height)
{
	// The triangular box is defined as a box with a given
	// triangular facet and a given height in the direction
	// of the normal vector of the facet.
	bool inside;
	CartesianCoordinate p = projection(point, triangle);
	if(isProjectionInsideTriangle(p, triangle))
	{
		double dist = distance(point, triangle);
		if(isEqual(dist, 0)) inside = true;
		else
		{
			if(isEqual(height, 0)) inside = false;
			else
			{
				double sides = dist * height;
				if(sides < 0)
				{
					// The location of the given point and the box region are
					// at the opposite sides of the given triangle.
					inside = false;
				}
				else
				{
					if(std::fabs(dist) <= std::fabs(height)) inside = true;
					else inside = false;
				}
			}
		}
	}
	else inside = false;
	return inside;
}

double angle(Vector& v1, Vector& v2)
{
	double alpha = std::acos(dotProd(v1, v2) / (abs(v1) * abs(v2)));
	if(isEqual(alpha, 0)) alpha = 0;
	else if(isEqual(alpha, M_PI)) alpha = M_PI;
	return alpha;
}

double angle(Triangle& t1, Triangle& t2)
{
	// This function calculates the angle between two triangles
	// sharing the same edge.
	const CartesianCoordinate* vertices1 = t1.getVertices();
	const CartesianCoordinate* vertices2 = t2.getVertices();
	bool match[3] = {false, false, false};
	size_t i, j;
	for(i = 0; i < 3; ++i)
	{
		for(j = 0; j < 3; ++j)
		{
			if(vertices1[j] == vertices2[i])
			{
				match[i] = true;
				break;
			}
		}
	}
	size_t match_count = 0;
	for(i = 0; i < 3; ++i) if(match[i] == true) ++match_count;
	assert(match_count == 2);
	for(i = 0; i < 3; ++i) if(match[i] == false) break;
	double dist = distance(vertices2[i], t1);
	assert(std::fabs(dist) > DBL_EPSILON);
	// These two triangular facets cannot be parallel to each other.
	double alpha = angle(t1.getNormal(), t2.getNormal());
	double beta;
	if(dist > 0) beta = M_PI - alpha;
	else beta = M_PI + alpha;
	if(beta > 2 * M_PI) beta -= 2 * M_PI;
	return beta;
}

CartesianCoordinate rotate(CartesianCoordinate r, Vector n, double alpha, bool dir)
{
	if(!dir) alpha = -alpha;
	double x = n.getX();
	double y = n.getY();
	double z = n.getZ();
	double sin_alpha = std::sin(alpha);
	double cos_alpha = std::cos(alpha);
	double R[3][3];
	R[0][0] = cos_alpha + (1 - cos_alpha) * x * x;
	R[0][1] = (1 - cos_alpha) * x * y + sin_alpha * z;
	R[0][2] = (1 - cos_alpha) * x * z - sin_alpha * y;
	R[1][0] = (1 - cos_alpha) * y * x - sin_alpha * z;
	R[1][1] = cos_alpha + (1 - cos_alpha) * y * y;
	R[1][2] = (1 - cos_alpha) * y * z + sin_alpha * x;
	R[2][0] = (1 - cos_alpha) * z * x + sin_alpha * y;
	R[2][1] = (1 - cos_alpha) * z * y - sin_alpha * x;
	R[2][2] = cos_alpha + (1 - cos_alpha) * z * z;
	double u = R[0][0] * r.x + R[0][1] * r.y + R[0][2] * r.z;
	double v = R[1][0] * r.x + R[1][1] * r.y + R[1][2] * r.z;
	double w = R[2][0] * r.x + R[2][1] * r.y + R[2][2] * r.z;
	if(std::fabs(u) < DBL_EPSILON) u = 0;
	if(std::fabs(v) < DBL_EPSILON) v = 0;
	if(std::fabs(w) < DBL_EPSILON) w = 0;
	CartesianCoordinate s(u, v, w);
	return s;
}

int strtoi(const std::string& str)
{
	return std::atoi(str.c_str());
}

long strtol(const std::string& str)
{
	return std::strtol(str.c_str(), 0, 0);
}

unsigned long strtoul(const std::string& str)
{
	return std::strtoul(str.c_str(), 0, 0);
}

float strtof(const std::string& str)
{
	return std::strtof(str.c_str(), 0);
}

double strtod(const std::string& str)
{
	return std::strtod(str.c_str(), 0);
}

std::string tolower(const std::string& str)
{
	std::string str2(str);
	std::transform(str.begin(), str.end(), str2.begin(), (int(*)(int))std::tolower);
	return str2;
}

std::string toupper(const std::string& str)
{
	std::string str2(str);
	std::transform(str.begin(), str.end(), str2.begin(), (int(*)(int))std::toupper);
	return str2;
}

bool strtob(const std::string& str)
{
	bool ans = false;
	std::string str2 = tolower(str);
	assert(str2 == "true" || str2 == "false" );
	if(str2 == "true") ans = true;
	else ans = false;
	return ans;
}

int uitoi(unsigned int i)
{
	int j = static_cast<int>(i);
	assert(j >= 0);
	return j;
}

long ultol(unsigned long i)
{
	long j = static_cast<long>(i);
	assert(j >= 0);
	return j;
}

}
