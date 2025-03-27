#ifndef TYPEDEFS_HPP_
#define TYPEDEFS_HPP_

#include <list>

namespace motility
{

class CartesianCoordinate;
typedef std::list<CartesianCoordinate> CartesianCoordinates;
typedef CartesianCoordinates::iterator CartesianCoordinateHandle;
typedef CartesianCoordinates::const_iterator CartesianCoordinateConstHandle;

class SphericalCoordinate;
typedef std::list<SphericalCoordinate> SphericalCoordinates;
typedef SphericalCoordinates::iterator SphericalCoordinateHandle;
typedef SphericalCoordinates::const_iterator SphericalCoordinateConstHandle;

class CylindricalCoordinate;
typedef std::list<CylindricalCoordinate> CylindricalCoordinates;
typedef CylindricalCoordinates::iterator CylindricalCoordinateHandle;
typedef CylindricalCoordinates::const_iterator CylindricalCoordinateConstHandle;

class Orientation;
typedef std::list<Orientation> Orientations;
typedef Orientations::iterator OrientationHandle;
typedef Orientations::const_iterator OrientationConstHandle;

class Molecule;
typedef std::list<Molecule> Molecules;
typedef Molecules::iterator MoleculeHandle;
typedef Molecules::const_iterator MoleculeConstHandle;

class Actin;
typedef std::list<Actin> Actins;
typedef Actins::iterator ActinHandle;
typedef Actins::const_iterator ActinConstHandle;

class CAP;
typedef std::list<CAP> CAPs;
typedef CAPs::iterator CAPHandle;
typedef CAPs::const_iterator CAPConstHandle;

class ARP23;
typedef std::list<ARP23> ARP23s;
typedef ARP23s::iterator ARP23Handle;
typedef ARP23s::const_iterator ARP23ConstHandle;

class ADF;
typedef std::list<ADF> ADFs;
typedef ADFs::iterator ADFHandle;
typedef ADFs::const_iterator ADFConstHandle;

struct Vertex;
typedef std::list<Vertex> Vertices;
typedef Vertices::iterator VertexHandle;
typedef Vertices::const_iterator VertexConstHandle;
typedef std::list<VertexHandle> VertexHandles;
typedef VertexHandles::iterator VertexHandleHandle;
typedef VertexHandles::const_iterator VertexHandleConstHandle;

struct Edge;
typedef std::list<Edge> Edges;
typedef Edges::iterator EdgeHandle;
typedef Edges::const_iterator EdgeConstHandle;
typedef std::list<EdgeHandle> EdgeHandles;
typedef EdgeHandles::iterator EdgeHandleHandle;
typedef EdgeHandles::const_iterator EdgeHandleConstHandle;

struct Facet;
typedef std::list<Facet> Facets;
typedef Facets::iterator FacetHandle;
typedef Facets::const_iterator FacetConstHandle;
typedef std::list<FacetHandle> FacetHandles;
typedef FacetHandles::iterator FacetHandleHandle;
typedef FacetHandles::const_iterator FacetHandleConstHandle;

class FilamentBranch;
typedef std::list<FilamentBranch> FilamentBranches;
typedef FilamentBranches::iterator FilamentBranchHandle;
typedef FilamentBranches::const_iterator FilamentBranchConstHandle;
typedef std::list<FilamentBranchHandle> FilamentBranchHandles;
typedef FilamentBranchHandles::iterator FilamentBranchHandleHandle;
typedef FilamentBranchHandles::const_iterator FilamentBranchHandleConstHandle;

class BranchTree;
typedef std::list<BranchTree> BranchTrees;
typedef BranchTrees::iterator BranchTreeHandle;
typedef BranchTrees::const_iterator BranchTreeConstHandle;
typedef std::list<BranchTreeHandle> BranchTreeHandles;
typedef BranchTreeHandles::iterator BranchTreeHandleHandle;
typedef BranchTreeHandles::const_iterator BranchTreeHandleConstHandle;

}

#endif /*TYPEDEFS_HPP_*/
