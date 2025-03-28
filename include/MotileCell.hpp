#ifndef MOTILECELL_HPP_
#define MOTILECELL_HPP_

#include <map>
#include <typedefs.hpp>
#include <DiscreteEventSimulator.hpp>
#include <BranchTree.hpp>
#include <FilamentBranch.hpp>
#include <UniformMolecularDistribution.hpp>
#include <SurfaceTopology.hpp>
#include <VertexEdgeFacet.hpp>

namespace motility
{

class CellStatisticsCalculator;

typedef simulation::DiscreteEvent_iterator FilamentReaction_iterator;
typedef simulation::DiscreteEvent_iterators FilamentReaction_iterators;
typedef simulation::DiscreteEvent_iterator_iterator FilamentReaction_iterator_iterator;

class MotileCell : public simulation::DiscreteEventSimulator
{
  private:

	// The pointer to extracellular signal distribution.
	UniformMolecularDistribution* ecs_dist;

	/// The actin cytoskeleton of motile cell.
	BranchTrees filament_network;

	/// The membrane surface of motile cell.
	SurfaceTopology membrane_surface;

	/// For simplicity, assume uniform distribution of actin,
	/// Arp23, CP and ADF.
	UniformMolecularDistribution* actin_dist;

	UniformMolecularDistribution* arp23_dist;

	UniformMolecularDistribution* cap_dist;

	UniformMolecularDistribution* adf_dist;

	/// The the name and the directory of geometry files.
	std::string data_dir, cell_geom_filename, cell_geom_filename_ext;

	/// The pointer to the calculator of cell statistics.
	CellStatisticsCalculator* cell_statistics_calculator;

  private:

	FilamentBranch makeNewFilament(double rou, double theta, double phi);

	void initializeFilamentNetwork();

	void initializeFilamentReaction();

	void associateNewFilamentReaction(FilamentReaction_iterator growing_reaction_ptr, FilamentReaction_iterator branching_reaction_ptr, FilamentReaction_iterator capping_reaction_ptr);

	void connect(FilamentReaction_iterator reaction_ptr);

	void create(FilamentReaction_iterator reaction_ptr);

	void pre_remove_event(FilamentReaction_iterator reaction_ptr, FilamentReaction_iterator destroyed_reaction_ptr);

	void pre_remove_event(FilamentReaction_iterator reaction_ptr);

	/// This function calculates the change of total energy of
	/// motile cell.
	///
	/// \param branch_handle the handle of an actin filament.
	/// \param type the type of a reaction.
	/// \return The change of the energy of total energy.
	double computeEnergyChange(FilamentBranchHandle branch_handle, const std::string& type);

	/// A set of overriden functions from DiscreteEventSimulator

	void step_record();

	void time_record();

	void initialize();

	void finalize();

  public:

	MotileCell(double max_duration, size_t max_step, double record_time_interval, size_t record_step_interval, UniformMolecularDistribution* ecsd_ptr, const std::string dir, const std::string geom_filename, const std::string geom_filename_ext);

	virtual ~MotileCell() throw();

	/// A set of utility functions to retrieve several internal
	/// data.

	UniformMolecularDistribution& getActinDist();

	UniformMolecularDistribution& getArp23Dist();

	UniformMolecularDistribution& getCapDist();

	UniformMolecularDistribution& getAdfDist();

	BranchTrees& getFilamentNetwork();

	SurfaceTopology& getMembraneSurface();

	/// A set of functions to calculate the rate of several
	/// filament reactions and to execute these reactions.
	/// These functions are called by FilamentReaction when
	/// corresponding reaction gets executed.

	/// This function calculates the rate of filament growing
	/// reaction.
	///
	/// \param branch_handle the handle of a filament.
	/// \return The rate of filament growing reaction.
	double computeFilamentGrowingRate(FilamentBranchHandle branch_handle);

	/// This function calculates the rate of filament branching
	/// reaction.
	///
	/// \param branch_handle the handle of a filament.
	/// \return The rate of filament branching reaction.
	double computeFilamentBranchingRate(FilamentBranchHandle branch_handle);

	/// This function calculates the rate of filament capping
	/// reaction.
	///
	/// \param branch_handle the handle of a filament.
	/// \return The rate of filament capping reaction.
	double computeFilamentCappingRate(FilamentBranchHandle branch_handle);

	/// This function updates the attachment status of capped
	/// filaments.
	///
	/// \param vertices a list of vertices.
	/// \param fibronectin_dist the spatial distribution of
	/// extracellular fibronectin.
	/// \return No value is returned.
	void updateCappedFilamentAttachmentToMembrane(VertexHandles& vertices);

	/// This function executes filament growing reaction.
	///
	/// \param branch_handle the handle of a filament.
	/// \param fibronectin_dist the spatial distribution of
	/// extracellular fibronectin.
	/// \return A list of vertices affected by filament growing.
	VertexHandles growFilament(FilamentBranchHandle branch_handle);

	/// This function executes filament branching reaction.
	///
	/// \param branch_handle the handle of a filament.
	/// \param fibronectin_dist the spatial distribution of
	/// extracellular fibronectin.
	/// \return A list of vertices affected by filament branching.
	VertexHandles branchFilament(FilamentBranchHandle branch_handle);

	/// This function executes filament capping reaction.
	///
	/// \param branch_handle the handle of a filament.
	/// \param fibronectin_dist the spatial distribution of
	/// extracellular fibronectin.
	/// \return A list of vertices affected by filament capping.
	VertexHandles capFilament(FilamentBranchHandle branch_handle);
};

}

#endif /*MOTILECELL_HPP_*/
