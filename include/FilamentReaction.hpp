#ifndef REACTION_HPP_
#define REACTION_HPP_

#include <cstdlib>
#include <cmath>
#include <typedefs.hpp>
#include <DiscreteEvent.hpp>
#include <MotileCell.hpp>
#include <UniformMolecularDistribution.hpp>

namespace motility
{

/// This class characterizes how fast a physical event proceeds
/// by two parameters: the rate of this event and the waiting
/// time before this event occurs.

class FilamentReaction : public simulation::DiscreteEvent
{
  private:

	std::string type;

	FilamentBranchHandle filament;

	MotileCell* cell;

	UniformMolecularDistribution* ecs;

	/// The filaments affected by the action of this reaction,
	/// also including itself. MotileCell uses this property to
	/// search for other reactions affected by the occurrence
	/// of this reaction.
	FilamentBranchHandles affected_filaments;

  private:

	double compute_rate();

  public:

	FilamentReaction(const std::string& t, FilamentBranchHandle f, MotileCell* c, UniformMolecularDistribution* e);

	virtual ~FilamentReaction() throw();

	const std::string& getType() const;

	FilamentBranchHandle getFilament();

	FilamentBranchHandles& getAffectedFilaments();

	void action();
};

}

#endif /*REACTION_HPP_*/
