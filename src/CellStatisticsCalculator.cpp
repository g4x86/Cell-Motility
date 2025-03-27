#include <CellStatisticsCalculator.hpp>
#include <MotileCell.hpp>
#include <ParameterTable.hpp>
#include <algorithms.hpp>

namespace motility
{

CellStatisticsCalculator::CellStatisticsCalculator(MotileCell* cell, UniformMolecularDistribution* ecs_dist, double t, const std::string& data_dir)
{
	// Step 2:
	// Initialize current time moment.
	time_moment = t;
	// Step 3:
	// Open output files.
	ParameterTable::Table& param_table = ParameterTable::instance();
	std::string cell_radius_dist_filename = data_dir + param_table[std::string("cell_radius_dist_filename")];
	cell_radius_dist_file = new OutputFile(cell_radius_dist_filename);
	std::string cell_velocity_dist_filename = data_dir + param_table[std::string("cell_velocity_dist_filename")];
	cell_velocity_dist_file = new OutputFile(cell_velocity_dist_filename);
	std::string cell_growing_dist_filename = data_dir + param_table[std::string("cell_growing_dist_filename")];
	cell_growing_dist_file = new OutputFile(cell_growing_dist_filename);
	std::string cell_devangle_dist_filename = data_dir + param_table[std::string("cell_devangle_dist_filename")];
	cell_devangle_dist_file = new OutputFile(cell_devangle_dist_filename);
	std::string cell_outward_dist_filename = data_dir + param_table[std::string("cell_outward_dist_filename")];
	cell_outward_dist_file = new OutputFile(cell_outward_dist_filename);
	std::string cell_stats_filename = data_dir + param_table[std::string("cell_stats_filename")];
	cell_stats_file = new OutputFile(cell_stats_filename);
	delimeter = param_table[std::string("delimeter")];
	// Step 4:
	// Calculate the number of periphery intervals.
	double periphery_degree_interval = strtod(param_table[std::string("periphery_degree_interval")]);
	double theta_min = strtod(param_table[std::string("theta_min")]);
	double theta_max = strtod(param_table[std::string("theta_max")]);
	double intervals = (theta_max - theta_min) / periphery_degree_interval;
	assert(intervals >= 0);
	double intervals_int = std::floor(intervals);
	n_periphery_interval = static_cast<size_t>(intervals_int);
	if(intervals > intervals_int) ++n_periphery_interval;
	// Step 5:
	// Initialize the radius distribution of a motile cell.
	radii = 0;
	radii = new double [n_periphery_interval];
	size_t* n_filament = 0;
	n_filament = new size_t [n_periphery_interval];
	size_t* n_growing_filament = 0;
	n_growing_filament = new size_t [n_periphery_interval];
	double* growing_percent = 0;
	growing_percent = new double [n_periphery_interval];
	double* deviation_angles = 0;
	deviation_angles = new double [n_periphery_interval];
	size_t* n_outward_filament = 0;
	n_outward_filament = new size_t [n_periphery_interval];
	double* outward_percent = 0;
	outward_percent = new double [n_periphery_interval];
	for(size_t i = 0; i < n_periphery_interval; ++i)
	{
		radii[i] = 0;
		n_filament[i] = 0;
		n_growing_filament[i] = 0;
		growing_percent[i] = 0;
		deviation_angles[i] = 0;
		n_outward_filament[i] = 0;
		outward_percent[i] = 0;
	}
	// Step 6:
	// Calculate the initial position of the leading edge in each
	// angular slice along cell periphery.
	collectSpreadingProperty(cell, ecs_dist, n_filament, radii, n_growing_filament, growing_percent, deviation_angles, n_outward_filament, outward_percent);
	size_t n_empty_periphery_radius = 0;
	radius = 0;
	for(size_t i = 0; i < n_periphery_interval; ++i)
	{
		if(n_filament[i] > 0) radius += radii[i];
		else
		{
			radii[i] = nan("");
			++n_empty_periphery_radius;
		}
	}
	assert(n_empty_periphery_radius < n_periphery_interval);
	radius /= (n_periphery_interval - n_empty_periphery_radius);
	// Step 7:
	// Initialize the output files for spreading-velocity distribution
	// and various statistics.
	// The first column of velocity distribution file is time moment,
	// and the rest is velocity distribution.
	// Get the output file stream.
	std::ofstream& cell_radius_dist_ostream = cell_radius_dist_file->getStream();
	std::ofstream& cell_velocity_dist_ostream = cell_velocity_dist_file->getStream();
	std::ofstream& cell_growing_dist_ostream = cell_growing_dist_file->getStream();
	std::ofstream& cell_devangle_dist_ostream = cell_devangle_dist_file->getStream();
	std::ofstream& cell_outward_dist_ostream = cell_outward_dist_file->getStream();
	std::ofstream& cell_stats_ostream = cell_stats_file->getStream();
	cell_radius_dist_ostream << theta_min;
	cell_velocity_dist_ostream << theta_min;
	cell_growing_dist_ostream << theta_min;
	cell_devangle_dist_ostream << theta_min;
	cell_outward_dist_ostream << theta_min;
	for(size_t i = 0; i < n_periphery_interval; ++i)
	{
		double theta = theta_min + (i + 1) * periphery_degree_interval;
		if(theta > theta_max) theta = theta_max;
		cell_radius_dist_ostream << delimeter << theta;
		cell_velocity_dist_ostream << delimeter << theta;
		cell_growing_dist_ostream << delimeter << theta;
		cell_devangle_dist_ostream << delimeter << theta;
		cell_outward_dist_ostream << delimeter << theta;
	}
	cell_radius_dist_ostream << std::endl;
	cell_velocity_dist_ostream << std::endl;
	cell_growing_dist_ostream << std::endl;
	cell_devangle_dist_ostream << std::endl;
	cell_outward_dist_ostream << std::endl;
	cell_radius_dist_ostream << time_moment;
	cell_growing_dist_ostream << time_moment;
	cell_devangle_dist_ostream << time_moment;
	cell_outward_dist_ostream << time_moment;
	for(size_t i = 0; i < n_periphery_interval; ++i)
	{
		cell_radius_dist_ostream << delimeter << radii[i];
		cell_growing_dist_ostream << delimeter << growing_percent[i];
		cell_devangle_dist_ostream << delimeter << deviation_angles[i];
		cell_outward_dist_ostream << delimeter << outward_percent[i];
	}
	cell_radius_dist_ostream << std::endl;
	cell_growing_dist_ostream << std::endl;
	cell_devangle_dist_ostream << std::endl;
	cell_outward_dist_ostream << std::endl;
	cell_stats_ostream << "Time" << delimeter;
	cell_stats_ostream << "SpreadingArea" << delimeter;
	cell_stats_ostream << "SpreadingVelocity" << delimeter;
	cell_stats_ostream << "nTotalFilaments" << delimeter;
	cell_stats_ostream << "nGrowingFilaments" << delimeter;
	cell_stats_ostream << "MembraneArea" << std::endl;
	// Step 8:
	// Clear up memory.
	if(n_outward_filament != 0) delete [] n_outward_filament;
	if(outward_percent != 0) delete [] outward_percent;
	if(deviation_angles != 0) delete [] deviation_angles;
	if(growing_percent != 0) delete [] growing_percent;
	if(n_growing_filament != 0) delete [] n_growing_filament;
	if(n_filament != 0) delete [] n_filament;
}

CellStatisticsCalculator::~CellStatisticsCalculator()
{
	if(radii != 0) delete [] radii;
	if(cell_radius_dist_file != 0) delete cell_radius_dist_file;
	if(cell_velocity_dist_file != 0) delete cell_velocity_dist_file;
	if(cell_growing_dist_file != 0) delete cell_growing_dist_file;
	if(cell_devangle_dist_file != 0) delete cell_devangle_dist_file;
	if(cell_outward_dist_file != 0) delete cell_outward_dist_file;
	if(cell_stats_file != 0) delete cell_stats_file;
}

void CellStatisticsCalculator::collectSpreadingProperty(MotileCell* cell, UniformMolecularDistribution* ecs_dist, size_t* n_filament, double* radii, size_t* n_growing_filament, double* growing_percent, double* deviation_angles, size_t* n_outward_filament, double* outward_percent)
{
	ParameterTable::Table& param_table = ParameterTable::instance();
	double periphery_degree_interval = strtod(param_table[std::string("periphery_degree_interval")]);
	SurfaceTopology& membrane_surface = cell->getMembraneSurface();
	Vertices& vertices = membrane_surface.getVertices();
	for(VertexHandle vh = vertices.begin(); vh != vertices.end(); ++vh)
	{
		CartesianCoordinate vertex_location(vh->getLocation());
		if(!isEqual(ecs_dist->getDensity(vertex_location), 0))
		{
			//
			// For a fibroblast cell spreading on a plane of glass slide, the
			// best characterization of spreading process is the distribution
			// of spreading velocity along cell periphery over 360 degree on the
			// plane. This spreading velocity is a 2D velocity vector. Therefore
			// the 'rou' of the cylindrical coordinate of the leading edge of
			// cell membrane can be used to determine the new_radius location
			// of the leading edge. The position of leading edge is estimated
			// by the filament end with the maximum radial coordinate.
			//
			CylindricalCoordinate p(vertex_location);
			volatile size_t i = static_cast<size_t>(std::floor(p.theta * 180 / (M_PI * periphery_degree_interval)));
			assert(i < n_periphery_interval);
			// The leading edge is represented by the most outside cell periphery.
			if(p.rou > radii[i]) radii[i] = p.rou;
			n_filament[i] += 1;
			// Count the number of growing filaments.
			if(!vh->getFilament()->isCapped()) n_growing_filament[i] += 1;
			// Record the maximum deviation angle.
			double deviation_angle = membrane_surface.computeDeviationAngleOfFilamentGrowth(*(vh->getFilament()));
			if(deviation_angle > deviation_angles[i]) deviation_angles[i] = deviation_angle;
			// Count outward growing filaments.
			if(deviation_angle < M_PI / 2) n_outward_filament[i] += 1;
		}
	}
	for(size_t i = 0; i < n_periphery_interval; ++i)
	{
		if(n_filament[i] > 0)
		{
			// Calculate the percentage of the number of growing filaments.
			growing_percent[i] = static_cast<double>(n_growing_filament[i]) / static_cast<double>(n_filament[i]);
			// Calculate the percentage of the number of outward filaments.
			outward_percent[i] = static_cast<double>(n_outward_filament[i]) / static_cast<double>(n_filament[i]);
		}
		else
		{
			// Set growing percentage to nan.
			growing_percent[i] = nan("");
			// Set deviation angle to nan.
			deviation_angles[i] = nan("");
			// Set outward percentage to nan.
			outward_percent[i] = nan("");
			// Approximate the position of leading edge from the positions of the
			// neighboring peripheral leading edges.
			size_t prev_index, next_index;
			if(i == 0) prev_index = n_periphery_interval - 1;
			else prev_index = i -1;
			if(i == n_periphery_interval - 1) next_index = 0;
			else next_index = i + 1;
			if(n_filament[prev_index] > 0 && n_filament[next_index] > 0) radii[i] = (radii[prev_index] + radii[next_index]) / 2;
			else radii[i] = nan("");
		}
	}
}

void CellStatisticsCalculator::computeSpreadingVelocity(double* radii, double* new_radii, double dt, double* spreading_velocity)
{
	for(size_t i = 0; i < n_periphery_interval; ++i)
	{
		// Calculate the spreading velocity.
		double v = 0;
		if(std::isnan(radii[i]) == 0 && std::isnan(new_radii[i]) == 0)
		{
			if(std::fabs(dt) > DBL_EPSILON)
			{
				v = (new_radii[i] - radii[i]) / dt;
				if(std::fabs(v) < DBL_EPSILON) v = 0;
			}
			else v = nan("");
		}
		else v = 0;
		spreading_velocity[i] = v;
	}
}

void CellStatisticsCalculator::operator()(MotileCell* cell, UniformMolecularDistribution* ecs_dist, double t, bool file_saving_flag)
{
	// Get the output file stream.
	std::ofstream& cell_radius_dist_ostream = cell_radius_dist_file->getStream();
	std::ofstream& cell_velocity_dist_ostream = cell_velocity_dist_file->getStream();
	std::ofstream& cell_growing_dist_ostream = cell_growing_dist_file->getStream();
	std::ofstream& cell_devangle_dist_ostream = cell_devangle_dist_file->getStream();
	std::ofstream& cell_outward_dist_ostream = cell_outward_dist_file->getStream();
	std::ofstream& cell_stats_ostream = cell_stats_file->getStream();
	// Step 0:
	// Initialize the variables to calculate new radius distribution.
	double dt = t - time_moment;
	if(dt > DBL_EPSILON)
	{
		double* new_radii = 0;
		new_radii = new double [n_periphery_interval];
		double* spreading_velocity = 0;
		spreading_velocity = new double [n_periphery_interval];
		size_t* n_filament = 0;
		n_filament = new size_t [n_periphery_interval];
		size_t* n_growing_filament = 0;
		n_growing_filament = new size_t [n_periphery_interval];
		double* growing_percent = 0;
		growing_percent = new double [n_periphery_interval];
		double* deviation_angles = 0;
		deviation_angles = new double [n_periphery_interval];
		size_t* n_outward_filament = 0;
		n_outward_filament = new size_t [n_periphery_interval];
		double* outward_percent = 0;
		outward_percent = new double [n_periphery_interval];
		for(size_t i = 0; i < n_periphery_interval; ++i)
		{
			new_radii[i] = 0;
			spreading_velocity[i] = 0;
			n_filament[i] = 0;
			n_growing_filament[i] = 0;
			growing_percent[i] = 0;
			deviation_angles[i] = 0;
			n_outward_filament[i] = 0;
			outward_percent[i] = 0;
		}
		// Step 1:
		// Collect and calculate various spreading properties, e.g. radius,
		// percentage of growing filaments, deviation angle.
		collectSpreadingProperty(cell, ecs_dist, n_filament, new_radii, n_growing_filament, growing_percent, deviation_angles, n_outward_filament, outward_percent);
		computeSpreadingVelocity(radii, new_radii, dt, spreading_velocity);
		if(file_saving_flag)
		{
			cell_radius_dist_ostream << t;
			cell_velocity_dist_ostream << t;
			cell_growing_dist_ostream << t;
			cell_devangle_dist_ostream << t;
			cell_outward_dist_ostream << t;
			for(size_t i = 0; i < n_periphery_interval; ++i)
			{
				cell_radius_dist_ostream << delimeter << new_radii[i];
				cell_velocity_dist_ostream << delimeter << spreading_velocity[i];
				cell_growing_dist_ostream << delimeter << growing_percent[i];
				cell_devangle_dist_ostream << delimeter << (180 * deviation_angles[i] / M_PI);
				cell_outward_dist_ostream << delimeter << outward_percent[i];
			}
			cell_radius_dist_ostream << std::endl;
			cell_velocity_dist_ostream << std::endl;
			cell_growing_dist_ostream << std::endl;
			cell_devangle_dist_ostream << std::endl;
			cell_outward_dist_ostream << std::endl;
		}
		// Step 2:
		// Calculate, record and print out various cell statistics.
		double new_radius = 0;
		size_t n_empty_periphery_radius = 0;
		for(size_t i = 0; i < n_periphery_interval; ++i)
		{
			if(n_filament[i] > 0) new_radius += new_radii[i];
			else ++n_empty_periphery_radius;
		}
		assert(n_empty_periphery_radius < n_periphery_interval);
		new_radius /= (n_periphery_interval - n_empty_periphery_radius);
		double average_spreading_velocity = 0;
		size_t n_empty_spreading_velocity = 0;
		for(size_t i = 0; i < n_periphery_interval; ++i)
		{
			if(!std::isnan(spreading_velocity[i])) average_spreading_velocity += spreading_velocity[i];
			else ++n_empty_spreading_velocity;
		}
		assert(n_empty_spreading_velocity < n_periphery_interval);
		average_spreading_velocity /= (n_periphery_interval - n_empty_spreading_velocity);
		double average_spreading_area = M_PI * new_radius * new_radius;
		SurfaceTopology& membrane_surface = cell->getMembraneSurface();
		if(file_saving_flag)
		{
			cell_stats_ostream << t << delimeter;
			cell_stats_ostream << average_spreading_area << delimeter;
			cell_stats_ostream << average_spreading_velocity << delimeter;
			cell_stats_ostream << membrane_surface.getVertexSize() << delimeter;
			cell_stats_ostream << membrane_surface.getVolatileVertexSize() << delimeter;
			cell_stats_ostream << membrane_surface.getArea() << std::endl;
		}
		std::cout << "dt = " << dt << std::endl;
		std::cout << "total number of filaments on membrane = " << membrane_surface.getVertexSize() << std::endl;
		std::cout << "total number of growing filaments on membrane = " << membrane_surface.getVolatileVertexSize() << std::endl;
		std::cout << "membrane surface area = " << membrane_surface.getArea() << std::endl;
		std::cout << "total filament density on membrane = " << membrane_surface.getVertexSize() / membrane_surface.getArea() << std::endl;
		std::cout << "average lamellipodium diameter = " << 2 * new_radius << std::endl;
		std::cout << "lamellipodium spreading area = " << average_spreading_area << std::endl;
		std::cout << "average spreading velocity = " << average_spreading_velocity << " um per sec" << std::endl;
		// Step 3:
		// Save current results of radius distribution.
		for(size_t i = 0; i < n_periphery_interval; ++i) radii[i] = new_radii[i];
		radius = new_radius;
		// Step 4:
		// Clean up memory.
		if(n_outward_filament != 0) delete [] n_outward_filament;
		if(outward_percent != 0) delete [] outward_percent;
		if(deviation_angles != 0) delete [] deviation_angles;
		if(growing_percent != 0) delete [] growing_percent;
		if(n_growing_filament != 0) delete [] n_growing_filament;
		if(n_filament != 0) delete [] n_filament;
		if(spreading_velocity != 0) delete [] spreading_velocity;
		if(new_radii != 0) delete [] new_radii;
		// Step 6:
		// Save current time moment
		time_moment = t;
	}
}

}
