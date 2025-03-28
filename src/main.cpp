/*******************************************************************************
 *
 *    main.cpp
 *
 *    The program to simulate actin-baesd cell motility.
 *
 ******************************************************************************/

#include <cstdlib>
#include <iostream>
#include <string>
#include <constants.hpp>
#include <algorithms.hpp>
#include <InputFile.hpp>
#include <ParameterTable.hpp>
#include <ReactionTypeTable.hpp>
#include <MotileCell.hpp>
#include <SpatialBoundary.hpp>
#include <UniformMolecularDistribution.hpp>
#include <initializeParameterTable.hpp>
#include <initializeReactionTypeTable.hpp>
#include <CellStatisticsCalculator.hpp>
#include <MotileCell.hpp>

using namespace motility;
using namespace simulation;

void usage_help(char* program_name)
{
	std::cout << "Usage:" << std::endl;
	std::cout << "If MOTILITY_HOME is not set: " << std::endl;
	std::cout << program_name << " <paramter INI file> <reaction INI file>" << std::endl;
	std::cout << "If MOTILITY_HOME is set, then use the default input files: " << "${MOTILITY_HOME}/input/paramters.ini and ${MOTILITY_HOME}/input/reactions.ini" << std::endl;
	exit(-1);
}

int main(int argc, char* argv[])
{
	// Obtain the home directory of cell motility simulation from
	// environment settings.
	char* home_dir_char = std::getenv("MOTILITY_HOME");
	std::string home_dir;
	if(home_dir_char != 0) home_dir = home_dir_char;
	// Obtain paramter file names from input arguments.
	std::string parameter_filename, reaction_filename;
	if(argc < 2)
	{
		if(!home_dir.empty())
		{
			parameter_filename = home_dir + "/input/" + "parameters.ini";
			reaction_filename = home_dir + "/input/" + "reactions.ini";
		}
		else usage_help(argv[0]);
	}
	else if(argc < 3) usage_help(argv[0]);
	else if(argc < 4)
	{
		parameter_filename = argv[1];
		reaction_filename = argv[2];
	}
	else usage_help(argv[0]);
	// Read the input parameter file and reaction file.
	InputFile parameter_file(parameter_filename);
	initializeParameterTable(parameter_file.getStream());
	InputFile reaction_file(reaction_filename);
	initializeReactionTypeTable(reaction_file.getStream());
	// Obtain a reference to the global paramter table file.
	ParameterTable::Table& param_table = ParameterTable::instance();
	// Set up the directory and the names of cell geometry files.
	std::string data_dir;
	if(!home_dir.empty()) data_dir = home_dir + "/output/";
	else data_dir = "./";
	std::string nameBuf = param_table[std::string("cell_geom_filename")];
	std::string cell_geom_filename, cell_geom_filename_ext;
	bool splitting_flag = splitFileName(nameBuf, cell_geom_filename, cell_geom_filename_ext);
	if(!splitting_flag) handleErrorEvent("the geometry file does not have an extension name");
	// Initialize spatially distributed extracellular signaling molecules.
	double x_min = strtod(param_table[std::string("x_min")]);
	double x_max = strtod(param_table[std::string("x_max")]);
	double y_min = strtod(param_table[std::string("y_min")]);
	double y_max = strtod(param_table[std::string("y_max")]);
	double z_min = strtod(param_table[std::string("z_min")]);
	double leading_edge_thickness = strtod(param_table[std::string("leading_edge_thickness")]);
	SpatialBoundary fibronectin_boundary(x_min, x_max, y_min, y_max, z_min, z_min + leading_edge_thickness);
	double fibronectin_conc = strtod(param_table[std::string("fibronectin_conc")]);
	UniformMolecularDistribution fibronectin_dist(fibronectin_boundary, fibronectin_conc);
	double simulation_time = strtod(param_table[std::string("simulation_time")]);
	size_t simulation_step = strtoul(param_table[std::string("simulation_step")]);
	double record_time_interval = strtod(param_table[std::string("record_time_interval")]);
	size_t record_step_interval = strtoul(param_table[std::string("record_step_interval")]);
	// Initialize motile cell.
	MotileCell motile_cell(simulation_time, simulation_step, record_time_interval, record_step_interval, &fibronectin_dist, data_dir, cell_geom_filename, cell_geom_filename_ext);
	// Start simulating actin-based cell motility.
	motile_cell.run();
	return motile_cell.get_status();
}
