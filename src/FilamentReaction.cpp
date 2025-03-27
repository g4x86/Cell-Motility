#include <cfloat>
#include <FilamentReaction.hpp>
#include <constants.hpp>
#include <algorithms.hpp>

namespace motility
{

FilamentReaction::FilamentReaction(const std::string& t, FilamentBranchHandle f, MotileCell* c, UniformMolecularDistribution* e) : simulation::DiscreteEvent()
{
	type = t;
	filament = f;
	cell = c;
	ecs = e;
	affected_filaments.clear();
	// Initiate the rate and period of this filament reactions.
	update();
}

FilamentReaction::~FilamentReaction() throw() {}

const std::string& FilamentReaction::getType() const
{
	return type;
}

FilamentBranchHandle FilamentReaction::getFilament()
{
	return filament;
}

double FilamentReaction::compute_rate()
{
	double r;
	bool outside_flag = isEqual(ecs->getDensity(filament->getTailEndLocation()), 0);
	SWITCH(type)
	{
		// If a filament grows out of extracellular signaling region,
		// set its capping rate to +inf and its growing and branching
		// rates to zero, such that it can be capped immediately.
		CASE("GROWING")
		{
			if(!outside_flag) r = cell->computeFilamentGrowingRate(filament);
			else r = 0;
			break;
		}
		CASE("BRANCHING")
		{
			if(!outside_flag) r = cell->computeFilamentBranchingRate(filament);
			else r = 0;
			break;
		}
		CASE("CAPPING")
		{
			if(!outside_flag) r = cell->computeFilamentCappingRate(filament);
			else r = DBL_INF_POSITIVE;
			break;
		}
		DEFAULT()
		{
			r = DBL_INF_POSITIVE;
		}
	}
	SWITCH_END()
	return r;
}

FilamentBranchHandles& FilamentReaction::getAffectedFilaments()
{
	return affected_filaments;
}

void FilamentReaction::action()
{
	VertexHandles affected_vertices;
	SWITCH(type)
	{
		CASE("GROWING")
		{
			affected_vertices = cell->growFilament(filament);
			break;
		}
		CASE("BRANCHING")
		{
			affected_vertices = cell->branchFilament(filament);
			break;
		}
		CASE("CAPPING")
		{
			affected_vertices = cell->capFilament(filament);
			break;
		}
		DEFAULT() {}
	}
	SWITCH_END()
	// Since the action of this reaction may cause the geometry change
	// of the local surface around this filament and therefore affect
	// neighboring filaments attached to local surface, it is important
	// to update the list of affect filaments of this reaction right
	// after its action is executed.
	affected_filaments.clear();
	for(VertexHandleHandle vhh = affected_vertices.begin(); vhh != affected_vertices.end(); ++vhh)
	{
		FilamentBranchHandle branch_handle = (*vhh)->getFilament();
		if(!branch_handle->isCapped()) affected_filaments.push_back(branch_handle);
	}
	// The action of this reaction also affects the filament that
	// this reaction is associated with.
    unique_append<FilamentBranchHandle>(affected_filaments, filament);
}

}
