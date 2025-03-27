#ifndef BRANCHTREE_HPP_
#define BRANCHTREE_HPP_

#include <typedefs.hpp>
#include <FilamentBranch.hpp>

namespace motility
{

class BranchTree
{
  private:

	FilamentBranches branches;

  public:

	BranchTree();

	BranchTree(const ARP23& arp23, const Actin& actin, const Orientation& orient);

	/// This constructor creates a new actin filament tree from
	/// a given filament.
	BranchTree(const FilamentBranch& branch);

	virtual ~BranchTree();

	/// This constructor adds a new actin filament into current
	/// filament tree.
	void addFilamentBranch(const ARP23& arp23, const Actin& actin, const Orientation& orient, FilamentBranchHandle parent_handle, BranchTreeHandle tree_handle);

	/// This constructor adds a new actin filament into current
	/// filament tree.
	void addFilamentBranch(const FilamentBranch& branch, FilamentBranchHandle parent_handle, BranchTreeHandle tree_handle);

	void removeFilamentBranch(FilamentBranchHandle fbh);

	bool isEmpty() const;

	FilamentBranches& getBranches();

	FilamentBranchHandle getFirstBranchHandle();

	FilamentBranchHandle getLastBranchHandle();
};

}

#endif /*BRANCHTREE_HPP_*/
