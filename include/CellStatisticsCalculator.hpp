#ifndef CELLSTATISTICSCALCULATOR_HPP_
#define CELLSTATISTICSCALCULATOR_HPP_

#include <string>
#include <UniformMolecularDistribution.hpp>
#include <OutputFile.hpp>

namespace motility
{

class MotileCell;

/// This class calculates various dynamic statistics of
/// a spreading cell and save the results into two files:
/// one containing various overall statistics of spreading
/// cell, and the other one containing the 2D distirubtion
/// of spreading velocity.
class CellStatisticsCalculator
{
  private:

	/// This variable records current time moment.
	double time_moment;

	/// This variable denotes the number of angular intervals along
	/// the periphery of the leading edge of a cell spreading on a
	/// 2D glass slide coated with fibronectin.
	volatile size_t n_periphery_interval;

	/// These two variables store cell radius and its angulr distribution.
	double radius;
	double* radii;

	/// The directory to store three output files.
	std::string data_dir;

	/// The radius distribution of the spreading cell.
	OutputFile* cell_radius_dist_file;

	/// The velocity distribution of the spreading cell.
	OutputFile* cell_velocity_dist_file;

	/// The distribution of the percentage of the number of growing
	/// filaments of the spreading cell.
	OutputFile* cell_growing_dist_file;

	/// The distribution of the deviation angle of filament growing
	/// direction from the raidal direction on glass slide.
	OutputFile* cell_devangle_dist_file;

	/// The distribution of the percentage of of the number of actin
	/// filaments growing outward.
	OutputFile* cell_outward_dist_file;

	/// The statistics of cell spreading.
	OutputFile* cell_stats_file;

	/// This variable specifies the delimiter used in output data files.
	std::string delimeter;

  private:

	void collectSpreadingProperty(MotileCell* cell, UniformMolecularDistribution* ecs_dist, size_t* n_filament, double* radii, size_t* n_growing_filament, double* growing_percent, double* deviation_angles, size_t* n_outward_filament, double* outward_percent);

	void computeSpreadingVelocity(double* radii, double* new_radii, double dt, double* spreading_velocity);

  public:

	/// The constructor initialize various variables needed to calculate cell statistics.
	CellStatisticsCalculator(MotileCell* cell, UniformMolecularDistribution* ecs_dist, double t, const std::string& data_dir);

	/// The destructor clean up used memory and close data files.
	virtual ~CellStatisticsCalculator();

	/// The function-like interface for main program to call.
	void operator()(MotileCell* cell, UniformMolecularDistribution* ecs_dist, double t, bool file_saving_flag);
};

}

#endif /*CELLSTATISTICSCALCULATOR_HPP_*/
