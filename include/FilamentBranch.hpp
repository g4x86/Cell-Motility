#ifndef FILAMENTBRANCH_HPP_
#define FILAMENTBRANCH_HPP_

#include <cassert>
#include <typedefs.hpp>
#include <Actin.hpp>
#include <ARP23.hpp>
#include <CAP.hpp>
#include <VertexEdgeFacet.hpp>
#include <Vector.hpp>
#include <Line.hpp>
#include <DiscreteEvent.hpp>

namespace motility
{

/// This class defines a number of components comprising actin
/// filmanet.
///
/// The data structure of actin filament contains all necessary
/// components to integrate individual filaments into filament
/// tree and to facilitate filament-related biochemical reactions.

class FilamentBranch
{
  private:

	Actins filament;

	/// The pointer to corresponding vertex on cell membrane
	/// surface.
	VertexHandle vertex;

	/// This variable indicates whether current actin filament
	/// is attached to cell membrane.
	bool membrane_attachment;

	/// This variable indicates whether current actin filament
	/// is linked to cell membrane.
	bool membrane_linkage;

	double branching_angle;

	Orientation orient;

	/// This class uses two pointers to manage the presence of
	/// Arp23 and capping protein for current actin filament,
	/// including their memory management.
	ARP23* arp23_handle;

	CAP* cap_handle;

	/// The pointer to the associated tree.
	BranchTreeHandle tree_handle;

	/// The pointer to the parent actin filament.
	FilamentBranchHandle parent_handle;

	/// The position of current actin filament in the list of
	/// the child filaments of the parent filament, starting
	/// from zero.
	size_t nth_child_of_parent;

	FilamentBranchHandles child_branches;

	/// The position to which current actin filament is attached
	/// on the parent filament, starting from zero.
	std::list<size_t> child_locations;

	/// The flag indicating the allowance of branching a child
	/// filament.
	bool child_branching_flag;

	/// The orientation of the child filament branch to be
	/// created.
	Orientation child_branch_orient;

	/// The facet handle which the child filament branch to be
	/// created
	/// will be added into.
	FacetHandle child_branch_facet;

	/// The initial lenght of this filament.
	double initial_length;

	/// These three variables are defined for the member function
	/// perturbTailEndLocation and restoreTailEndLocation.
	bool virtual_tail_end;

	CartesianCoordinate virtual_tail_end_location;

	double virtual_tail_end_diameter;

	simulation::DiscreteEvent_iterators reactions;

  private:

	void initializeMemory(const ARP23& arp23);

	void clearMemory();

	void copyMemory(const FilamentBranch& fb);

	CartesianCoordinate getLocationOfTheOtherEnd(const CartesianCoordinate& loc, double dist, const Orientation& ot) const;

  public:

	FilamentBranch();

	FilamentBranch(const ARP23& arp23, const Actin& actin, const Orientation& ot);

	FilamentBranch(const FilamentBranch& fb);

	~FilamentBranch();

	VertexHandle getVertex();

	void setVertex(VertexHandle vh);

	bool isAttachedToMembrane() const;

	bool isLinkedToMembrane() const;

	bool addActin(const Actin& actin);

	bool removeActin();

	/// This function artificially modifies the location of filament tail end.
	///
	/// \param mol a Molecule instance whose diameter is needed to calculate tail end.
	/// \return The execuation status of this function.
	bool perturbTailEndLocation(const Molecule& mol);

	/// This function retores the location of filament tail end.
	///
	/// \return The execuation status of this function.
	bool restoreTailEndLocation();

	bool addCap(const CAP& cap);

	bool removeCap();

	bool removeArp23();

	bool isArp23ed() const;

	bool isCapped() const;

	ARP23* getArp23Pointer();

	CartesianCoordinate getArp23Location() const;

	CAP* getCapPointer();

	CartesianCoordinate getCapLocation() const;

	bool isEmpty() const;

	bool isFilamentEmpty() const;

	Actins& getFilament();

	size_t length() const;

	BranchTreeHandle getTreeHandle();

	void setTreeHandle(BranchTreeHandle th);

	FilamentBranchHandle getParentHandle();

	void setParentHandle(FilamentBranchHandle ph);

	size_t getNthChildOfParent();

	void setNthChildOfParent(size_t n);

	size_t sizeofChildHandles();

	FilamentBranchHandles& getChildHandles();

	void addChildHandle(FilamentBranchHandle childHandle);

	std::list<size_t>& getChildLocations();

	size_t getLastChildLocation() const;

	void addChildLocation(size_t childLoc);

	bool isMinimalLengthForBranchingReached() const;

	bool isBranchingAllowed() const;

	ActinConstHandle getBranchingSiteActinConstHandle() const;

	CartesianCoordinate getBranchingSiteActinLocation() const;

	CartesianCoordinate getHeadEndLocation() const;

	CartesianCoordinate getTailEndLocation() const;

	double getTailEndDiameter() const;

	const Orientation& getOrient() const;

	bool getChildBranchingFlag() const;

	void setChildBranchingFlag(bool f);

	Orientation getChildBranchOrient() const;

	void setChildBranchOrient(Orientation orient);

	FacetHandle getChildBranchFacet() const;

	void setChildBranchFacet(FacetHandle fh);

	double getInitialLength() const;

	simulation::DiscreteEvent_iterators& getReactions();

	void addReaction(simulation::DiscreteEvent_iterator reac);

	void removeReaction(simulation::DiscreteEvent_iterator reac);

	FilamentBranch& operator=(const FilamentBranch& fb);
};

inline FilamentBranchHandle filament_branch_sentinel {};
inline FilamentBranchConstHandle filament_branch_const_sentinel {};
inline FilamentBranchHandleHandle filament_branch_handle_sentinel {};
inline FilamentBranchHandleConstHandle filament_branch_handle_const_sentinel {};

}

#endif /*FILAMENTBRANCH_HPP_*/
