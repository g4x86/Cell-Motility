#include <BranchTree.hpp>
#include <VertexEdgeFacet.hpp>
#include <ParameterTable.hpp>
#include <algorithms.hpp>

namespace motility
{

BranchTree::BranchTree() {}

BranchTree::~BranchTree() {}

BranchTree::BranchTree(const ARP23& arp23, const Actin& actin, const Orientation& orient)
{
	branches.push_back(FilamentBranch(arp23, actin, orient));
}

BranchTree::BranchTree(const FilamentBranch& branch)
{
	branches.push_back(branch);
}

void BranchTree::addFilamentBranch(const ARP23& arp23, const Actin& actin, const Orientation& orient, FilamentBranchHandle parent_handle, BranchTreeHandle tree_handle)
{
	ARP23 newArp23(arp23);
	newArp23.setLocation(parent_handle->getBranchingSiteActinLocation());
	addFilamentBranch(FilamentBranch(newArp23, actin, orient), parent_handle, tree_handle);
}

void BranchTree::addFilamentBranch(const FilamentBranch& branch, FilamentBranchHandle parent_handle, BranchTreeHandle tree_handle)
{
	size_t arp23_binding_actins = strtoul((ParameterTable::instance())[std::string("arp23_binding_actins")]);
	branches.push_back(branch);
	FilamentBranchHandle childHandle = --(branches.end());
	childHandle->setParentHandle(parent_handle);
	childHandle->setTreeHandle(tree_handle);
	childHandle->setNthChildOfParent(parent_handle->sizeofChildHandles());
	parent_handle->addChildHandle(childHandle);
	long child_location_on_filament = parent_handle->length() - arp23_binding_actins / 2 - 1;
	// 'child_location_on_filament' is calculated in a way such
	// that an ARP23 molecule binds to 'arp23_binding_actins'
	// number of polymerized actin monomers on mother filament
	// starting from filament barbed end immediately.
	assert(child_location_on_filament >= 0);
	parent_handle->addChildLocation(child_location_on_filament);
}

void BranchTree::removeFilamentBranch(FilamentBranchHandle fbh)
{
	//
	// To be implemented:
	//
	// Update the interlinking relationships between this filament
	// and other filaments before removing it completely.
	//
	//
	// Set the filament pointer of the corresponding vertex, if it
	// exists, to 0 before remove the actin filament.
	if(fbh->getVertex() != static_cast<VertexHandle>(nullptr)) fbh->getVertex()->setFilament(static_cast<FilamentBranchHandle>(nullptr));
	// Now remove the given actin filament from the filament pool.
	branches.erase(fbh);
}

bool BranchTree::isEmpty() const
{
	return branches.empty();
}

FilamentBranches& BranchTree::getBranches()
{
	return branches;
}

FilamentBranchHandle BranchTree::getFirstBranchHandle()
{
	FilamentBranchHandle fbh;
	if(!branches.empty()) fbh = branches.begin();
	else fbh = static_cast<FilamentBranchHandle>(nullptr);
	return fbh;
}

FilamentBranchHandle BranchTree::getLastBranchHandle()
{
	FilamentBranchHandle fbh;
	if(!branches.empty()) fbh = (--branches.end());
	else fbh = static_cast<FilamentBranchHandle>(nullptr);
	return fbh;
}

}
