#include <iostream>
#include <string>
#include <cstdlib>
#include <unistd.h>
#include <MotileCell.hpp>
#include <CellStatisticsCalculator.hpp>
#include <constants.hpp>
#include <OutputFile.hpp>
#include <ParameterTable.hpp>
#include <ReactionTypeTable.hpp>
#include <FilamentReaction.hpp>
#include <algorithms.hpp>

namespace motility
{

MotileCell::MotileCell(double max_duration, size_t max_step, double record_time_interval, size_t record_step_interval, UniformMolecularDistribution* ecsd_ptr, const std::string dir, const std::string geom_filename, const std::string geom_filename_ext) : simulation::DiscreteEventSimulator(max_duration, max_step, record_time_interval, record_step_interval)
{
	ecs_dist = ecsd_ptr;
	data_dir = dir;
	cell_geom_filename = geom_filename;
	cell_geom_filename_ext = geom_filename_ext;
	ParameterTable::Table& param_table = ParameterTable::instance();
	double x_min = strtod(param_table[std::string("x_min")]);
	double x_max = strtod(param_table[std::string("x_max")]);
	double y_min = strtod(param_table[std::string("y_min")]);
	double y_max = strtod(param_table[std::string("y_max")]);
	double z_min = strtod(param_table[std::string("z_min")]);
	double z_max = strtod(param_table[std::string("z_max")]);
	SpatialBoundary intracellular_molecule_boundary(x_min, x_max, y_min, y_max, z_min, z_max);
	actin_dist = 0;
	double actin_conc = strtod(param_table[std::string("actin_conc")]);
	actin_dist = new UniformMolecularDistribution(intracellular_molecule_boundary, actin_conc);
	arp23_dist = 0;
	double arp23_conc = strtod(param_table[std::string("arp23_conc")]);
	arp23_dist = new UniformMolecularDistribution(intracellular_molecule_boundary, arp23_conc);
	cap_dist = 0;
	double cap_conc = strtod(param_table[std::string("cap_conc")]);
	cap_dist = new UniformMolecularDistribution(intracellular_molecule_boundary, cap_conc);
	adf_dist = 0;
	double adf_conc = strtod(param_table[std::string("adf_conc")]);
	adf_dist = new UniformMolecularDistribution(intracellular_molecule_boundary, adf_conc);
	cell_statistics_calculator = 0;
}

MotileCell::~MotileCell() throw()
{
	if(cell_statistics_calculator != 0) delete cell_statistics_calculator;
	if(actin_dist != 0) delete actin_dist;
	if(arp23_dist != 0) delete arp23_dist;
	if(cap_dist != 0) delete cap_dist;
	if(adf_dist != 0) delete adf_dist;
}

FilamentBranch MotileCell::makeNewFilament(double rou, double theta, double phi)
{
	ParameterTable::Table& param_table = ParameterTable::instance();
	double actin_diameter = strtod(param_table[std::string("actin_diameter")]);
	double arp23_diameter = strtod(param_table[std::string("arp23_diameter")]);
	Orientation orient(theta, phi);
	SphericalCoordinate local_position(rou, theta, phi);
	double init_cell_center_x = strtod(param_table[std::string("init_cell_center_x")]);
	double init_cell_center_y = strtod(param_table[std::string("init_cell_center_y")]);
	double init_cell_center_z = strtod(param_table[std::string("init_cell_center_z")]);
	CartesianCoordinate init_cell_center(init_cell_center_x, init_cell_center_y, init_cell_center_z);
	CartesianCoordinate actual_position = local_position.toCartesianCoordinate() + init_cell_center;
	ARP23 arp23(0, arp23_diameter, actual_position);
	Actin actin("ATP", actin_diameter);
	FilamentBranch branch(arp23, actin, orient);
	return branch;
}

void MotileCell::initializeFilamentNetwork()
{
	ParameterTable::Table& param_table = ParameterTable::instance();
	double init_cell_diameter = strtod(param_table[std::string("init_cell_diameter")]);
	size_t n_init_cell_mesh_horizontal = strtoul(param_table[std::string("n_init_cell_mesh_horizontal")]);
	size_t n_init_cell_mesh_vertical = strtoul(param_table[std::string("n_init_cell_mesh_vertical")]);
	double delta_theta = 2 * M_PI / n_init_cell_mesh_horizontal;
	double delta_phi = M_PI / n_init_cell_mesh_vertical;
	// Vary the value of theta and phi by delta_theta and delta_phi
	// respectively to construct a triangulated spherical surface.
	BranchTreeHandle new_tree_handle;
	FilamentBranchHandle new_branch_handle;
	double theta = 0, phi = 0;
	// The double loops are constructed to create all vertices needed
	// for the initial spherical geometry of fibroblast cell and then
	// connect the vertices into triangular facets.
	for(size_t layer_vertical = 0; layer_vertical <= n_init_cell_mesh_vertical; ++layer_vertical)
	{
		phi = delta_phi * layer_vertical;
		if(layer_vertical == 0)
		{
			// The vertex at the top of sphere
			filament_network.push_back(BranchTree(makeNewFilament(init_cell_diameter/2, theta, phi)));
			// a new BranchTree is created.
			new_tree_handle = --(filament_network.end());
			new_branch_handle = new_tree_handle->getLastBranchHandle();
			// Set the associated tree of the new filament branch.
			new_branch_handle->setTreeHandle(new_tree_handle);
			// Everytime when a new filament is added to cell surface as a
			// a vertex, this filament and corresponding vertex must be
			// bi-directionally linked to each other.
			//
			// 1) the vertex holds a FilamentBranchHandle-type pointer to the
			// filament. This operation is carried out within Vertex constructor.
			//
			// 2) the filament also holds a VertexHandle-type pointer to
			// the vertex. This operation is carried out by the member
			// function 'addVertex' of Membrane when this vertex is created
			// and added to cell surface.
			membrane_surface.addVertex(Vertex(new_branch_handle));
		}
		else if(layer_vertical == n_init_cell_mesh_vertical)
		{
			// The vertex at the bottom of sphere
			filament_network.push_back(BranchTree(makeNewFilament(init_cell_diameter/2, theta, phi)));
			new_tree_handle = --(filament_network.end());
			new_branch_handle = new_tree_handle->getLastBranchHandle();
			// Set the associated tree of the new filament branch.
			new_branch_handle->setTreeHandle(new_tree_handle);
			// Bi-directionally link the newly created filament and vertex.
			membrane_surface.addVertex(Vertex(new_branch_handle));
			// Now add corresponding triangular facets into cell surface
			VertexHandle vh1, vh2, vh3;
			vh2 = membrane_surface.getLastVertexHandle();
			for(size_t layer_horizontal = n_init_cell_mesh_horizontal; layer_horizontal > 0; --layer_horizontal)
			{
				vh1 = vh2;
				advance(vh1, -ultol(layer_horizontal));
				vh3 = vh1;
				if(layer_horizontal == 1) advance(vh3, 1 - ultol(n_init_cell_mesh_horizontal));
				else advance(vh3, 1);
				membrane_surface.addFacet(vh1, vh2, vh3);
			}
		}
		else
		{
			// The vertices at in the middle of sphere
			for(size_t layer_horizontal = 0; layer_horizontal < n_init_cell_mesh_horizontal; ++layer_horizontal)
			{
				theta = delta_theta * layer_horizontal;
				filament_network.push_back(BranchTree(makeNewFilament(init_cell_diameter/2, theta, phi)));
				new_tree_handle = --(filament_network.end());
				new_branch_handle = new_tree_handle->getLastBranchHandle();
				// Set the associated tree of the new filament branch.
				new_branch_handle->setTreeHandle(new_tree_handle);
				// Bi-directionally link the newly created filament and vertex.
				membrane_surface.addVertex(Vertex(new_branch_handle));
				// Now add corresponding triangular facets into cell surface
				VertexHandle vh1, vh2, vh3;
				if(layer_vertical == 1)
				{
					// The vertex at the first circle layer of sphere
					vh2 = membrane_surface.getLastVertexHandle();
					vh3 = membrane_surface.getFirstVertexHandle();
					if(layer_horizontal > 0)
					{
						vh1 = vh2;
						advance(vh1, -1);
						membrane_surface.addFacet(vh1, vh2, vh3);
					}
					if(layer_horizontal == n_init_cell_mesh_horizontal - 1)
					{
						vh1 = vh2;
						advance(vh2, 1 - ultol(n_init_cell_mesh_horizontal));
						membrane_surface.addFacet(vh1, vh2, vh3);
					}
				}
				else
				{
					// The vertex at all other vertical layers of sphere,
					// except for those at the top layer, the first layer
					// and the bottom layer.
					if(layer_horizontal == 0)
					{
						// The vertex at the first horizontal layer
						//
						// Type I triangular facet
						vh1 = membrane_surface.getLastVertexHandle();
						vh2 = vh1;
						advance(vh2, 1 - ultol(n_init_cell_mesh_horizontal));
						vh3 = vh1;
						advance(vh3, -ultol(n_init_cell_mesh_horizontal));
						membrane_surface.addFacet(vh1, vh2, vh3);
					}
					else
					{
						// The vertex at the middle horizontal layers
						//
						// Type II triangular facet
						vh2 = membrane_surface.getLastVertexHandle();
						vh1 = vh2;
						advance(vh1, -1);
						vh3 = vh2;
						advance(vh3, -ultol(n_init_cell_mesh_horizontal));
						membrane_surface.addFacet(vh1, vh2, vh3);
						//
						// Type I triangular facet
						if(layer_horizontal < n_init_cell_mesh_horizontal - 1)
						{
							vh1 = vh2;
							advance(vh2, 1 - ultol(n_init_cell_mesh_horizontal));
							vh3 = vh2;
							advance(vh3, -1);
							membrane_surface.addFacet(vh1, vh2, vh3);
						}
						else
						{
							vh1 = membrane_surface.getLastVertexHandle();
							vh3 = vh1;
							advance(vh3, -ultol(n_init_cell_mesh_horizontal));
							vh2 = vh3;
							advance(vh2, 1 - ultol(n_init_cell_mesh_horizontal));
							membrane_surface.addFacet(vh1, vh2, vh3);
						}
					}
					if(layer_horizontal == n_init_cell_mesh_horizontal - 1)
					{
						// The vertex at the last horizontal layer
						//
						// Type II triangular facet
						vh2 = vh1;
						advance(vh2, 1 - ultol(n_init_cell_mesh_horizontal));
						vh3 = vh2;
						advance(vh3, -ultol(n_init_cell_mesh_horizontal));
						membrane_surface.addFacet(vh1, vh2, vh3);
					}
				}
			}
		}
	}
}

void MotileCell::initializeFilamentReaction()
{
	for(BranchTreeHandle tree_handle = filament_network.begin(); tree_handle != filament_network.end(); ++tree_handle)
	{
		FilamentBranches& branches = tree_handle->getBranches();
		for(FilamentBranchHandle branch_handle = branches.begin(); branch_handle != branches.end(); ++branch_handle)
		{
			assert(branch_handle->isAttachedToMembrane());
			if(!branch_handle->isCapped())
			{
				FilamentReaction_iterator growing_reaction_ptr = add_event(new FilamentReaction("GROWING", branch_handle, this, ecs_dist));
				FilamentReaction_iterator branching_reaction_ptr = add_event(new FilamentReaction("BRANCHING", branch_handle, this, ecs_dist));
				FilamentReaction_iterator capping_reaction_ptr = add_event(new FilamentReaction("CAPPING", branch_handle, this, ecs_dist));
				branch_handle->addReaction(growing_reaction_ptr);
				branch_handle->addReaction(branching_reaction_ptr);
				branch_handle->addReaction(capping_reaction_ptr);
				associateNewFilamentReaction(growing_reaction_ptr, branching_reaction_ptr, capping_reaction_ptr);
			}
		}
	}
}

void MotileCell::associateNewFilamentReaction(FilamentReaction_iterator growing_reaction_ptr, FilamentReaction_iterator branching_reaction_ptr, FilamentReaction_iterator capping_reaction_ptr)
{
	// Link each reaction to its affected reactions.
	add_modification(growing_reaction_ptr, growing_reaction_ptr);
	add_modification(growing_reaction_ptr, branching_reaction_ptr);
	add_modification(growing_reaction_ptr, capping_reaction_ptr);
	add_modification(branching_reaction_ptr, growing_reaction_ptr);
	add_modification(branching_reaction_ptr, branching_reaction_ptr);
	add_modification(branching_reaction_ptr, capping_reaction_ptr);
	add_destruction(capping_reaction_ptr, growing_reaction_ptr);
	add_destruction(capping_reaction_ptr, branching_reaction_ptr);
	add_destruction(capping_reaction_ptr, capping_reaction_ptr);
}

void MotileCell::connect(FilamentReaction_iterator reaction_ptr)
{
	// Remove all existing modification and destruction connections
	// of a reaction.
	empty_modified_events(reaction_ptr);
	empty_destroyed_events(reaction_ptr);
	// Add the reactions associated with all affected filaments to
	// the modified and/or the destroyed reactions of this reaction.
	FilamentReaction* reaction = dynamic_cast<FilamentReaction*>(*reaction_ptr);
	FilamentBranchHandle filament_ptr = reaction->getFilament();
	const std::string& reaction_type = reaction->getType();
	FilamentBranchHandles& affected_filament_ptrs = reaction->getAffectedFilaments();
	for(FilamentBranchHandleHandle affected_filament_ptr_ptr = affected_filament_ptrs.begin(); affected_filament_ptr_ptr != affected_filament_ptrs.end(); ++affected_filament_ptr_ptr)
	{
		FilamentReaction_iterators& affected_filament_reaction_ptrs = (*affected_filament_ptr_ptr)->getReactions();
		for(FilamentReaction_iterator_iterator affected_filament_reaction_ptr_ptr = affected_filament_reaction_ptrs.begin(); affected_filament_reaction_ptr_ptr != affected_filament_reaction_ptrs.end(); ++affected_filament_reaction_ptr_ptr)
		{
			if(*affected_filament_ptr_ptr == filament_ptr)
			{
				if(reaction_type == "GROWING" || reaction_type == "BRANCHING") add_modification(reaction_ptr, *affected_filament_reaction_ptr_ptr);
				else if(reaction_type == "CAPPING") add_destruction(reaction_ptr, *affected_filament_reaction_ptr_ptr);
				else {}
			}
			else add_modification(reaction_ptr, *affected_filament_reaction_ptr_ptr);
		}
	}
}

void MotileCell::create(FilamentReaction_iterator reaction_ptr)
{
	FilamentReaction* reaction = dynamic_cast<FilamentReaction*>(*reaction_ptr);
	if(reaction->getType() == "BRANCHING")
	{
		FilamentBranchHandle filament_ptr = reaction->getFilament();
		BranchTreeHandle tree_ptr = filament_ptr->getTreeHandle();
		// Since this function is called right after a branching reaction
		// is executed, the last filmanent branch in branch tree is the
		// newly created filament. Also remember to associate these new
		// filament reactions with the reactions on the mother filament.
		FilamentBranchHandle new_branch_ptr = tree_ptr->getLastBranchHandle();
		FilamentReaction_iterator growing_reaction_ptr = add_event(new FilamentReaction("GROWING", new_branch_ptr, this, ecs_dist));
		FilamentReaction_iterator branching_reaction_ptr = add_event(new FilamentReaction("BRANCHING", new_branch_ptr, this, ecs_dist));
		FilamentReaction_iterator capping_reaction_ptr = add_event(new FilamentReaction("CAPPING", new_branch_ptr, this, ecs_dist));
		new_branch_ptr->addReaction(growing_reaction_ptr);
		new_branch_ptr->addReaction(branching_reaction_ptr);
		new_branch_ptr->addReaction(capping_reaction_ptr);
		associateNewFilamentReaction(growing_reaction_ptr, branching_reaction_ptr, capping_reaction_ptr);
	}
}
void MotileCell::pre_remove_event(FilamentReaction_iterator reaction_ptr, FilamentReaction_iterator destroyed_reaction_ptr)
{
	// Remove the reaction to be destroyed from the reaction list
	// of the filament to which this reaction belongs.
	FilamentReaction* destroyed_reaction = dynamic_cast<FilamentReaction*>(*destroyed_reaction_ptr);
	destroyed_reaction->getFilament()->getReactions().remove(destroyed_reaction_ptr);
}

void MotileCell::pre_remove_event(FilamentReaction_iterator reaction_ptr)
{
	// Remove the reaction to be destroyed from the reaction list
	// of the filament to which this reaction belongs.
	FilamentReaction* reaction = dynamic_cast<FilamentReaction*>(*reaction_ptr);
	reaction->getFilament()->getReactions().remove(reaction_ptr);
}

void MotileCell::step_record()
{
	simulation::DiscreteEventSimulator::step_record();
	std::cout << "At the step #" << loop_step << " (time " << time_moment << " sec) :" << std::endl;
	(*cell_statistics_calculator)(this, ecs_dist, time_moment, false);
	std::cout << std::endl;
}

void MotileCell::time_record()
{
	simulation::DiscreteEventSimulator::time_record();
	std::cout << "At the time " << time_moment << " sec (step #" << loop_step << ") :" << std::endl;
	// Record a snapshot of cell geometry.
	std::stringstream strs;
	strs << data_dir << cell_geom_filename << '-' << time_moment << 's' << '.' << cell_geom_filename_ext;
	std::string current_geom_filename = strs.str();
	std::cout << "Writing " << current_geom_filename << std::endl;
	OutputFile current_geom_file(current_geom_filename);
	membrane_surface.exportGeometry(current_geom_file.getStream());
	(*cell_statistics_calculator)(this, ecs_dist, time_moment, true);
	std::cout << std::endl;
}

void MotileCell::initialize()
{
	// Initialize filament network.
	initializeFilamentNetwork();
	// Initialize filament reaction.
	initializeFilamentReaction();
	// Sort the initial event pool.
	simulation::DiscreteEventSimulator::initialize();
	// Record the initial cell geometry.
	std::string cell_geom_start_filename = data_dir + cell_geom_filename + "-start" + '.' + cell_geom_filename_ext;
	std::cout << "Writing " << cell_geom_start_filename << std::endl;
	OutputFile cell_geom_start_file(cell_geom_start_filename);
	membrane_surface.exportGeometry(cell_geom_start_file.getStream());
	// Initialize the calculator.
	cell_statistics_calculator = new CellStatisticsCalculator(this, ecs_dist, time_moment, data_dir);
	// Determine whether to initialize randome generator randomly.
	ParameterTable::Table& param_table = ParameterTable::instance();
	bool random_seed = strtob(param_table[std::string("random_seed")]);
	if(random_seed) srandom(time(0) * getpid());
	else srandom(1);
	std::cout << "Start simulation..." << std::endl;
}

void MotileCell::finalize()
{
	simulation::DiscreteEventSimulator::finalize();
	// Record the final cell geometry.
	std::string cell_geom_end_filename = data_dir + cell_geom_filename + "-end" + '.' + cell_geom_filename_ext;
	std::cout << "Writing " << cell_geom_end_filename << std::endl;
	OutputFile cell_geom_end_file(cell_geom_end_filename);
	membrane_surface.exportGeometry(cell_geom_end_file.getStream());
	// Print some overall simulation information.
	std::cout << "Simulation is done! " << std::endl;
	std::cout << "Total simulation time is " << time_moment << " seconds and simulation step is " << loop_step << "." << std::endl;
}

UniformMolecularDistribution& MotileCell::getActinDist()
{
	return (*actin_dist);
}

UniformMolecularDistribution& MotileCell::getArp23Dist()
{
	return (*arp23_dist);
}

UniformMolecularDistribution& MotileCell::getCapDist()
{
	return (*cap_dist);
}

UniformMolecularDistribution& MotileCell::getAdfDist()
{
	return (*adf_dist);
}

BranchTrees& MotileCell::getFilamentNetwork()
{
	return filament_network;
}

SurfaceTopology& MotileCell::getMembraneSurface()
{
	return membrane_surface;
}

double MotileCell::computeEnergyChange(FilamentBranchHandle branch_handle, const std::string& type)
{
	double energy_change = 0;
	ParameterTable::Table& param_table = ParameterTable::instance();
	double filament_membrane_resistance_pressure = strtod(param_table[std::string("filament_membrane_resistance_pressure")]);
	double actin_diameter = strtod(param_table[std::string("actin_diameter")]);
	if(type == "GROWING" || type == "CAPPING")
	{
		Vector filament_vector(1, branch_handle->getOrient());
		Vector centered_direct_area = membrane_surface.computeCenteredDirectionalAreaOfLocalSurface(branch_handle->getVertex());
		Vector total_resistance_force = centered_direct_area * filament_membrane_resistance_pressure;
		energy_change = dotProd(total_resistance_force, -filament_vector) * actin_diameter / 2;
	}
	else if(type == "BRANCHING")
	{
		Vector filament_vector(branch_handle->getInitialLength(), branch_handle->getChildBranchOrient());
		Line child_branch_line(branch_handle->getBranchingSiteActinLocation(), filament_vector);
		Vector centered_direct_area = membrane_surface.computeCenteredDirectionalAreaOfLocalSurface(child_branch_line.getEnd(), branch_handle->getChildBranchFacet());
		Vector total_resistance_force = centered_direct_area * filament_membrane_resistance_pressure;
		energy_change = dotProd(total_resistance_force, -filament_vector) * actin_diameter / 2;
	}
	else {}
	return energy_change;
}

double MotileCell::computeFilamentGrowingRate(FilamentBranchHandle branch_handle)
{
	ReactionTypeTable::Table& reac_table = ReactionTypeTable::instance();
	double growing_rate_const = (reac_table[std::string("growing")]).forward_const;
	double growing_rate;
	if(growing_rate_const < DBL_EPSILON) growing_rate = 0;
	else if(growing_rate_const < DBL_INF_POSITIVE)
	{
		// Calculate resistance factor.
		double energy_change = computeEnergyChange(branch_handle, "GROWING");
		ParameterTable::Table& param_table = ParameterTable::instance();
		double kT = strtod(param_table[std::string("kT")]);
		double resistance_factor = 1;
		if(energy_change > DBL_EPSILON) resistance_factor = std::exp(-energy_change / kT);
		double actin_conc = actin_dist->getDensity(branch_handle->getTailEndLocation());
		growing_rate = growing_rate_const * actin_conc * resistance_factor;
	}
	else growing_rate = DBL_INF_POSITIVE;
	return growing_rate;
}

double MotileCell::computeFilamentBranchingRate(FilamentBranchHandle branch_handle)
{
	ReactionTypeTable::Table& reac_table = ReactionTypeTable::instance();
	double branching_rate_const = (reac_table[std::string("branching")]).forward_const;
	double branching_rate;
	if(branching_rate_const < DBL_EPSILON) branching_rate = 0;
	else if(branching_rate_const < DBL_INF_POSITIVE)
	{
		// It is more efficient to determine if a branching reaction is allowed
		// by first checking the length of mother filament and then checking the
		// local surface around the mother filament, than the other sequence of
		// checking.
		if(branch_handle->isBranchingAllowed() && membrane_surface.isBranchingAllowed(*branch_handle))
		{
			// Calculate resistance factor.
			double energy_change = computeEnergyChange(branch_handle, "BRANCHING");
			ParameterTable::Table& param_table = ParameterTable::instance();
			double kT = strtod(param_table[std::string("kT")]);
			double resistance_factor = 1;
			if(energy_change > DBL_EPSILON) resistance_factor = std::exp(-energy_change / kT);
			// Caluclate the rate of filament branching reaction.
			double arp23_conc = arp23_dist->getDensity(branch_handle->getTailEndLocation());
			branching_rate = branching_rate_const * arp23_conc * resistance_factor;
			double actin_conc = actin_dist->getDensity(branch_handle->getTailEndLocation());
			size_t branching_actin_quantity = strtoul(param_table[std::string("branching_actin_quantity")]);
			for(size_t i = 0; i < branching_actin_quantity; ++i) branching_rate *= actin_conc;
		}
		else branching_rate = 0;
	}
	else branching_rate = DBL_INF_POSITIVE;
	return branching_rate;
}

double MotileCell::computeFilamentCappingRate(FilamentBranchHandle branch_handle)
{
	ReactionTypeTable::Table& reac_table = ReactionTypeTable::instance();
	double capping_rate_const = (reac_table[std::string("capping")]).forward_const;
	double capping_rate;
	if(capping_rate_const < DBL_EPSILON) capping_rate = 0;
	else if(capping_rate_const < DBL_INF_POSITIVE)
	{
		// Calculate resistance factor.
		// Add spatial constraint to the calculation of capping rate
		// such that any filaments growing inward will be capped right
		// away. The inclusion of this constaint is to mimic similar
		// constraint imposed on filament growth by membrane surface
		// during membrane protrusion in which any filaments growing
		// inward will cause the clash of cell membrane and therefore
		// raise up membrane energy significantly.
		// Calculate resistance factor.
		double energy_change = computeEnergyChange(branch_handle, "CAPPING");
		ParameterTable::Table& param_table = ParameterTable::instance();
		double kT = strtod(param_table[std::string("kT")]);
		double resistance_factor = 1;
		if(energy_change > DBL_EPSILON) resistance_factor = std::exp(-energy_change / kT);
		// Caluclate the rate of filament capping reaction.
		double cap_conc = cap_dist->getDensity(branch_handle->getTailEndLocation());
		capping_rate = capping_rate_const * cap_conc * resistance_factor;
	}
	else capping_rate = DBL_INF_POSITIVE;
	return capping_rate;
}

void MotileCell::updateCappedFilamentAttachmentToMembrane(VertexHandles& vertices)
{
	/// If the local surface around a capped vertex has a concave shape,
	/// this vertex and its underlying capped filament must be removed
	/// from membrane surface and filament network. All the neighboring
	/// vertices affected by this removal are collected such that the
	/// associated reactions can be updated later.
	/// The filament whose reaction triggers this updating process is
	/// excluded from this list of affected vertices even if it is capped
	/// by its capping reaction because the affected vertices calculated
	/// by capFilament does not include this filament itself. Therefore
	/// this filament remains valid through subsequent update of reacton
	/// interactions.
	VertexHandles affected_vertices;
	VertexHandleHandle vhh = vertices.begin();
	while(vhh != vertices.end())
	{
		bool vertex_removed_flag = false;
		FilamentBranchHandle branch_handle = (*vhh)->getFilament();
		if(!isEqual(ecs_dist->getDensity(branch_handle->getTailEndLocation()), 0) && branch_handle->isCapped())
		{
			double surface_energy_change = computeEnergyChange(branch_handle, "GROWING");
			if(surface_energy_change < -DBL_EPSILON)
			{
				// If no resistance force is imposed on this capped filament,
				// remove corresponding vertex from membrane surface and also
				// remove this filament from branch tree totally. Before this
				// capped filament and its vertex are removed, collect all
				// neighboring filaments of this filament because the removal
				// also affects the geometry of these neighboring filaments.
				// These filaments will be added into the list of affected
				// vertices later.
				FilamentBranchHandle branch_handle_removed = (*vhh)->getFilament();
				VertexHandles neighboring_vertices = membrane_surface.removeVertex(*vhh);
				affected_vertices = merge<VertexHandle>(affected_vertices, neighboring_vertices);
				affected_vertices.remove(*vhh);
				BranchTreeHandle tree_handle_removed = branch_handle_removed->getTreeHandle();
				tree_handle_removed->removeFilamentBranch(branch_handle_removed);
				if(tree_handle_removed->isEmpty()) filament_network.erase(tree_handle_removed);
				VertexHandleHandle vhh_removed = vhh++;
				vertices.erase(vhh_removed);
				vertex_removed_flag = true;
			}
		}
		if(!vertex_removed_flag) ++vhh;
	}
	// Update the list of affected vertices.
	vertices = merge<VertexHandle>(vertices, affected_vertices);
}

VertexHandles MotileCell::growFilament(FilamentBranchHandle branch_handle)
{
	/// This function executes filament growing reaction by adding
	/// one actin monomer to a filament and updating the geometry
	/// of the local surface of the filament.
	FilamentBranch& branch = *branch_handle;
	ParameterTable::Table& param_table = ParameterTable::instance();
	double actin_diameter = strtod(param_table[std::string("actin_diameter")]);
	bool action = branch.addActin(Actin("ATP", actin_diameter));
	assert(action);
	membrane_surface.updateCompositeProperties(branch.getVertex(), true, false);
	VertexHandles affected_vertices = membrane_surface.updateLocalSurface(branch.getVertex());
	updateCappedFilamentAttachmentToMembrane(affected_vertices);
	return affected_vertices;
}

VertexHandles MotileCell::branchFilament(FilamentBranchHandle branch_handle)
{
	/// This function executes filament branching reaction by creating
	/// a new actin filament, attaching it to a filament, and updating
	/// the geometry of the local surface of both filaments.
	// Everytime when a new filament is created, it must be linked to
	// a vertex on cell membrane_surface. This can be done by calling the member
	// function 'setVertex' of FilamentBranch through the member function
	// of Membrane which is responsible for adding new vertex into the
	// geometry of cell membrane_surface.
	//
	// Design principles:
	//
	// To separate filament-related operations and geometry-related
	// operations as much as possible by:
	//
	// 1) Move the filament-related operations to 'FilamentBranch' class.
	// 2) Move the geometry-related operations to 'Membrane' class.
	// 3) This function dispatches various tasks to corresponding classes
	// and only contains the codes interacting with both filament and
	// geometry.
	FilamentBranch& branch = *branch_handle;
	BranchTreeHandle tree_handle = branch.getTreeHandle();
	ParameterTable::Table& param_table = ParameterTable::instance();
	double actin_diameter = strtod(param_table[std::string("actin_diameter")]);
	double arp23_diameter = strtod(param_table[std::string("arp23_diameter")]);
	// 1) Determine the probability distribution of making a child
	// filament branch towards the neighboring membrane facets of
	// its mother filament, by examing the spatial relationship
	// between the mother filament and its neighboring facets.
	// 2) Select the branching direction of the child filament based
	// on the calculated probability distribution.
	assert(branch.getChildBranchingFlag());
	tree_handle->addFilamentBranch(ARP23(0, arp23_diameter), Actin("ATP", actin_diameter), branch.getChildBranchOrient(), branch_handle, tree_handle);
	FilamentBranchHandle child_branch_handle = tree_handle->getLastBranchHandle();
	// Add the child vertex into cell membrane_surface.
	VertexHandle child_vertex_handle = membrane_surface.addVertex(Vertex(child_branch_handle));
	VertexHandles affected_vertices = membrane_surface.updateLocalSurface(child_vertex_handle, branch.getChildBranchFacet());
	updateCappedFilamentAttachmentToMembrane(affected_vertices);
	branch.setChildBranchingFlag(false);
	return affected_vertices;
}

VertexHandles MotileCell::capFilament(FilamentBranchHandle branch_handle)
{
	/// This function executes filament capping reaction by adding
	/// a capping protein to an uncapped filament and updating the
	/// geometry of the local surface of the filament.
	FilamentBranch& branch = *branch_handle;
	VertexHandle vertex_handle = branch.getVertex();
	ParameterTable::Table& param_table = ParameterTable::instance();
	double cap_diameter = strtod(param_table[std::string("cap_diameter")]);
	bool action = branch.addCap(CAP(1, cap_diameter));
	assert(action);
	membrane_surface.updateCompositeProperties(vertex_handle, true, true);
	VertexHandles affected_vertices = membrane_surface.updateLocalSurface(vertex_handle);
	updateCappedFilamentAttachmentToMembrane(affected_vertices);
	return affected_vertices;
}

}
