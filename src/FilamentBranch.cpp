#include <cassert>
#include <FilamentBranch.hpp>
#include <BranchTree.hpp>
#include <VertexEdgeFacet.hpp>
#include <Coordinate.hpp>
#include <algorithms.hpp>
#include <ParameterTable.hpp>

namespace motility
{

FilamentBranch::FilamentBranch()
{
	vertex = static_cast<VertexHandle>(nullptr);
	membrane_attachment = false;
	membrane_linkage = false;
	arp23_handle = nullptr;
	cap_handle = nullptr;
	tree_handle = static_cast<BranchTreeHandle>(nullptr);
	parent_handle = static_cast<FilamentBranchHandle>(nullptr);
	nth_child_of_parent = 0;
	branching_angle = 0;
	initial_length = 0;
	child_branching_flag = false;
	child_branch_facet = static_cast<FacetHandle>(nullptr);
	virtual_tail_end = false;
	virtual_tail_end_diameter = 0;
}

FilamentBranch::FilamentBranch(const ARP23& arp23, const Actin& actin, const Orientation& ot)
{
	initializeMemory(arp23);
	arp23_handle->immobilize();
	orient = ot;
	vertex = static_cast<VertexHandle>(nullptr);
	membrane_attachment = false;
	membrane_linkage = false;
	tree_handle = static_cast<BranchTreeHandle>(nullptr);
	parent_handle = static_cast<FilamentBranchHandle>(nullptr);
	nth_child_of_parent = 0;
	ParameterTable::Table& param_table = ParameterTable::instance();
	double angle = strtod(param_table[std::string("branching_angle")]);
	branching_angle = angle * M_PI / 180;
	virtual_tail_end = false;
	virtual_tail_end_diameter = 0;
	// One Arp23 molecule and 'branching_actin_quantity' number of
	// actin molecules are needed to create a new actin filament.
	size_t branching_actin_quantity = strtoul(param_table[std::string("branching_actin_quantity")]);
	for(size_t i = 0; i < branching_actin_quantity; ++i) assert(addActin(actin));
	initial_length = distance(getHeadEndLocation(), getTailEndLocation());
}

FilamentBranch::FilamentBranch(const FilamentBranch& fb)
{
	copyMemory(fb);
}

FilamentBranch::~FilamentBranch()
{
	clearMemory();
}

void FilamentBranch::initializeMemory(const ARP23& arp23)
{
	arp23_handle = new ARP23;
	*arp23_handle = arp23;
	cap_handle = nullptr;
}

void FilamentBranch::clearMemory()
{
	if(arp23_handle != nullptr)
	{
		delete arp23_handle;
		arp23_handle = nullptr;
	}
	if(cap_handle != nullptr)
	{
		delete cap_handle;
		cap_handle = nullptr;
	}
}

void FilamentBranch::copyMemory(const FilamentBranch& fb)
{
	if(fb.arp23_handle != nullptr)
	{
		arp23_handle = new ARP23;
		*arp23_handle = *(fb.arp23_handle);
	}
	else arp23_handle = nullptr;
	if(fb.cap_handle != nullptr)
	{
		cap_handle = new CAP;
		*cap_handle = *(fb.cap_handle);
	}
	else cap_handle = nullptr;
	orient = fb.orient;
	filament = fb.filament;
	vertex = fb.vertex;
	membrane_attachment = fb.membrane_attachment;
	membrane_linkage = fb.membrane_linkage;
	tree_handle = fb.tree_handle;
	parent_handle = fb.parent_handle;
	nth_child_of_parent = fb.nth_child_of_parent;
	child_branches = fb.child_branches;
	child_locations = fb.child_locations;
	branching_angle = fb.branching_angle;
	child_branching_flag = fb.child_branching_flag;
	child_branch_orient = fb.child_branch_orient;
	child_branch_facet = fb.child_branch_facet;
	initial_length = fb.initial_length;
	virtual_tail_end = fb.virtual_tail_end;
	virtual_tail_end_location = fb.virtual_tail_end_location;
	virtual_tail_end_diameter = fb.virtual_tail_end_diameter;
	reactions = fb.reactions;
}

CartesianCoordinate FilamentBranch::getLocationOfTheOtherEnd(const CartesianCoordinate& loc, double dist, const Orientation& ot) const
{
	Line l(loc, Vector(dist, ot));
	return l.getEnd();
}

VertexHandle FilamentBranch::getVertex()
{
	return vertex;
}

void FilamentBranch::setVertex(VertexHandle vh)
{
	vertex = vh;
	if(vertex != static_cast<VertexHandle>(nullptr)) membrane_attachment = true;
	else membrane_attachment = false;
}

bool FilamentBranch::isAttachedToMembrane() const
{
	return membrane_attachment;
}

bool FilamentBranch::isLinkedToMembrane() const
{
	return membrane_linkage;
}

bool FilamentBranch::addActin(const Actin& actin)
{
	bool action;
	if(cap_handle == nullptr)
	{
		CartesianCoordinate last_tail_end_loc = getTailEndLocation();
		double last_tail_end_diam = getTailEndDiameter();
		filament.push_back(actin);
		Actin& a = filament.back();
		a.setLocation(getLocationOfTheOtherEnd(last_tail_end_loc, last_tail_end_diam, orient));
		a.immobilize();
		action = true;
	}
	else action = false;
	return action;
}

bool FilamentBranch::removeActin()
{
	bool action;
	if(!filament.empty())
	{
		filament.pop_back();
		if(cap_handle != nullptr)
		{
			if(!filament.empty())
			{
				Actin& tail_end_actin = filament.back();
				cap_handle->setLocation(getLocationOfTheOtherEnd(tail_end_actin.getLocation(), tail_end_actin.getDiameter() / 2, orient));
			}
			else
			{
				assert(arp23_handle != nullptr);
				cap_handle->setLocation(getLocationOfTheOtherEnd(arp23_handle->getLocation(), arp23_handle->getDiameter(), orient));
			}
			action = true;
		}
		else action = false;
	}
	else action = false;
	return action;
}

bool FilamentBranch::perturbTailEndLocation(const Molecule& mol)
{
	/// This function is provided for the pertubation method to
	/// calculate the change of membrane surface energy.
	bool action;
	if(!virtual_tail_end)
	{
		virtual_tail_end_diameter = mol.getDiameter();
		virtual_tail_end_location = getLocationOfTheOtherEnd(getTailEndLocation(), virtual_tail_end_diameter, orient);
		virtual_tail_end = true;
		action = true;
	}
	else action = false;
	return action;
}

bool FilamentBranch::restoreTailEndLocation()
{
	/// This function is provided for the pertubation method to
	/// calculate the change of membrane surface energy.
	bool action;
	if(virtual_tail_end)
	{
		virtual_tail_end_diameter = 0;
		virtual_tail_end = false;
		action = true;
	}
	else action = false;
	return action;
}

bool FilamentBranch::addCap(const CAP& cap)
{
	bool action;
	if(cap_handle == nullptr)
	{
		cap_handle = new CAP;
		*cap_handle = cap;
		if(filament.empty())
		{
			assert(arp23_handle != nullptr);
			cap_handle->setLocation(getLocationOfTheOtherEnd(arp23_handle->getLocation(), arp23_handle->getDiameter(), orient));
		}
		else
		{
			Actin& tail_end_actin = filament.back();
			cap_handle->setLocation(getLocationOfTheOtherEnd(tail_end_actin.getLocation(), tail_end_actin.getDiameter() / 2, orient));
		}
		cap_handle->immobilize();
		action = true;
	}
	else action = false;
	return action;
}

bool FilamentBranch::removeCap()
{
	bool action;
	if(cap_handle != nullptr)
	{
		delete cap_handle;
		cap_handle = nullptr;
		action = true;
	}
	else action = false;
	return action;
}

bool FilamentBranch::removeArp23()
{
	bool action;
	if(arp23_handle != nullptr)
	{
		delete arp23_handle;
		arp23_handle = nullptr;
		action = true;
	}
	else action = false;
	return action;
}

bool FilamentBranch::isArp23ed() const
{
	return ((arp23_handle != nullptr) ? true : false);
}

bool FilamentBranch::isCapped() const
{
	return ((cap_handle != nullptr) ? true : false);
}

ARP23* FilamentBranch::getArp23Pointer()
{
	return arp23_handle;
}

CartesianCoordinate FilamentBranch::getArp23Location() const
{
	assert(arp23_handle != nullptr);
	return arp23_handle->getLocation();
}

CAP* FilamentBranch::getCapPointer()
{
	return cap_handle;
}

CartesianCoordinate FilamentBranch::getCapLocation() const
{
	assert(cap_handle != nullptr);
	return cap_handle->getLocation();
}

bool FilamentBranch::isEmpty() const
{
	return (!arp23_handle && !cap_handle && filament.empty());
}

bool FilamentBranch::isFilamentEmpty() const
{
	return filament.empty();
}

Actins& FilamentBranch::getFilament()
{
	return filament;
}

size_t FilamentBranch::length() const
{
	return filament.size();
}

BranchTreeHandle FilamentBranch::getTreeHandle()
{
	return tree_handle;
}

void FilamentBranch::setTreeHandle(BranchTreeHandle th)
{
	tree_handle = th;
}

FilamentBranchHandle FilamentBranch::getParentHandle()
{
	return parent_handle;
}

void FilamentBranch::setParentHandle(FilamentBranchHandle ph)
{
	parent_handle = ph;
}

size_t FilamentBranch::getNthChildOfParent()
{
	return nth_child_of_parent;
}

void FilamentBranch::setNthChildOfParent(size_t n)
{
	nth_child_of_parent = n;
}

size_t FilamentBranch::sizeofChildHandles()
{
	return child_branches.size();
}

FilamentBranchHandles& FilamentBranch::getChildHandles()
{
	return child_branches;
}

void FilamentBranch::addChildHandle(FilamentBranchHandle childHandle)
{
	child_branches.push_back(childHandle);
}

std::list<size_t>& FilamentBranch::getChildLocations()
{
	return child_locations;
}

size_t FilamentBranch::getLastChildLocation() const
{
	size_t cLoc;
	if(!child_locations.empty()) cLoc = child_locations.back();
	else cLoc = 0;
	return cLoc;
}

void FilamentBranch::addChildLocation(size_t childLoc)
{
	child_locations.push_back(childLoc);
}

bool FilamentBranch::isMinimalLengthForBranchingReached() const
{
	bool isReached;
	ParameterTable::Table& param_table = ParameterTable::instance();
	size_t arp23_binding_actins = strtoul(param_table[std::string("arp23_binding_actins")]);
	size_t minimal_filament_length = arp23_binding_actins + 1;
	if(filament.size() < minimal_filament_length) isReached = false;
	else isReached = true;
	return isReached;
}

bool FilamentBranch::isBranchingAllowed() const
{
	ParameterTable::Table& param_table = ParameterTable::instance();
	bool branching = true;
	size_t nChild = child_branches.size();
	if(nChild == 0)
	{
		// Determine if current actin filament is long enough to
		// accommodate child branch.
		if(!isMinimalLengthForBranchingReached()) branching = false;
	}
	else
	{
		// Determine if the child branch to be created is far away
		// enough from previously created neighboring child branch
		// of current mother filament.
		size_t arp23_binding_actins = strtoul(param_table[std::string("arp23_binding_actins")]);
		size_t branching_site_location = filament.size() - arp23_binding_actins / 2 - 1;
		// 'branching_site_location' is the location of ARP23-binding
		// site on mother filament. It is counted from the pointed end
		// of mother filament and it starts from 0.
		long dist = branching_site_location - getLastChildLocation();
		assert(dist >= 0);
		if(static_cast<size_t>(dist) < arp23_binding_actins) branching = false;
	}
	return branching;
}

ActinConstHandle FilamentBranch::getBranchingSiteActinConstHandle() const
{
	ActinConstHandle ah;
	if(isMinimalLengthForBranchingReached())
	{
		ParameterTable::Table& param_table = ParameterTable::instance();
		size_t arp23_binding_actins = strtoul(param_table[std::string("arp23_binding_actins")]);
		size_t branching_site_location = arp23_binding_actins / 2 + 1;
		// 'branching_site_location' is the location of ARP23-binding
		// site on mother filament. It is counted from the barbed end
		// of mother filament and it starts from 1.
		ah = filament.end();
		for(size_t i = 0; i < branching_site_location; i++) --ah;
	}
	else ah = static_cast<ActinConstHandle>(nullptr);
	return ah;
}

CartesianCoordinate FilamentBranch::getBranchingSiteActinLocation() const
{
	ActinConstHandle ah = getBranchingSiteActinConstHandle();
	assert(ah != static_cast<ActinConstHandle>(nullptr));
	return ah->getLocation();
}

CartesianCoordinate FilamentBranch::getHeadEndLocation() const
{
	CartesianCoordinate loc;
	if(arp23_handle != nullptr)
	{
		loc = arp23_handle->getLocation();
	}
	else
	{
		assert(filament.empty() == false);
		loc = filament.front().getLocation();
	}
	return loc;
}

CartesianCoordinate FilamentBranch::getTailEndLocation() const
{
	CartesianCoordinate loc;
	if(!virtual_tail_end)
	{
		if(cap_handle != nullptr)
		{
			loc = getLocationOfTheOtherEnd(cap_handle->getLocation(), cap_handle->getDiameter(), orient);
		}
		else
		{
			if(filament.empty())
			{
				assert(arp23_handle != nullptr);
				loc = getLocationOfTheOtherEnd(arp23_handle->getLocation(), arp23_handle->getDiameter(), orient);
			}
			else
			{
				const Actin& actin = filament.back();
				loc = getLocationOfTheOtherEnd(actin.getLocation(), actin.getDiameter() / 2, orient);
			}
		}
		assert(finite(loc.x) && finite(loc.y) && finite(loc.z));
	}
	else loc = virtual_tail_end_location;
	return loc;
}

double FilamentBranch::getTailEndDiameter() const
{
	double diam;
	if(!virtual_tail_end)
	{
		if(cap_handle != nullptr) diam = cap_handle->getDiameter();
		else
		{
			if(filament.empty())
			{
				assert(arp23_handle != nullptr);
				diam = arp23_handle->getDiameter();
			}
			else
			{
				const Actin& actin = filament.back();
				diam = actin.getDiameter() / 2;
			}
		}
	}
	else diam = virtual_tail_end_diameter;
	return diam;
}

const Orientation& FilamentBranch::getOrient() const
{
	return orient;
}

bool FilamentBranch::getChildBranchingFlag() const
{
	return child_branching_flag;
}

void FilamentBranch::setChildBranchingFlag(bool f)
{
	child_branching_flag = f;
}

Orientation FilamentBranch::getChildBranchOrient() const
{
	return child_branch_orient;
}

void FilamentBranch::setChildBranchOrient(Orientation orient)
{
	child_branch_orient = orient;
}

FacetHandle FilamentBranch::getChildBranchFacet() const
{
	return child_branch_facet;
}

void FilamentBranch::setChildBranchFacet(FacetHandle fh)
{
	child_branch_facet = fh;
}

double FilamentBranch::getInitialLength() const
{
	return initial_length;
}

simulation::DiscreteEvent_iterators& FilamentBranch::getReactions()
{
	return reactions;
}

void FilamentBranch::addReaction(simulation::DiscreteEvent_iterator reac)
{
	unique_append<simulation::DiscreteEvent_iterator>(reactions, reac);
}

void FilamentBranch::removeReaction(simulation::DiscreteEvent_iterator reac)
{
	reactions.remove(reac);
}

FilamentBranch& FilamentBranch::operator=(const FilamentBranch& fb)
{
	if(this != &fb)
	{
		clearMemory();
		copyMemory(fb);
	}
	return *this;
}

}
