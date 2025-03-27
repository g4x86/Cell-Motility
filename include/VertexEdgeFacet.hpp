#ifndef VERTEXEDGEFACET_HPP_
#define VERTEXEDGEFACET_HPP_

#include <typedefs.hpp>
#include <Coordinate.hpp>
#include <Vector.hpp>
#include <Triangle.hpp>

namespace motility
{

/// The definition of Half-edge data structure
///
/// Half-edge data structure is used to store the topological
/// information of polyhedral surface.
///
/// Ref. Kettner, L. (1998) Designing of a data structure for
///      polyhedral surfaces.

class Vertex
{
	/// The link to actin filament
	FilamentBranchHandle filament;

	/// The link to the incident edges clockwisely oriented
	EdgeHandles edges;

	/// A vertex can have multiple incident edges and these
	/// edges may not be adjacent to each other when holes
	/// are created after the neighboring facet of this vertex
	/// is removed. Therefore this member data is indispensable
	/// for referencing the neighboring facets on cell surface
	/// of a given filament. Edges are inserted and removed by
	/// the member function addFacet() of Membrane which also
	/// sorts these edges clockwised.

  public:

	Vertex();

	Vertex(FilamentBranchHandle branch);

	FilamentBranchHandle getFilament();

	void setFilament(FilamentBranchHandle fbh);

	EdgeHandles& getEdges();

	CartesianCoordinate getLocation() const;

	bool operator==(const Vertex& v) const;

	bool operator!=(const Vertex& v) const;

	friend class Edge;

	friend class Facet;

	friend class SurfaceTopology;
};

class Edge
{
	/// The link to the incident vertex
	VertexHandle vertex;

	/// The link to the previous edge of the facet to which the
	/// current edge belongs
	EdgeHandle prev;

	/// The link to the next edge of the facet to which the
	/// current edge belongs
	EdgeHandle next;

	/// The link to the dual edge
	EdgeHandle dual;

	/// The beginning vertex of the dual edge is the ending one
	/// of the current edge, and rhe ending vertex of the dual
	/// edge is the beginning one of the current edge. Therefore
	/// the dual edge points to the opposite direction of the
	/// current edge. If the current edge is at the boundary of
	/// a polydedral surface, then the dual edge is 0. The link
	/// to the incident facet which is always at the left side
	/// of edge.
	FacetHandle facet;

  public:

	Edge();

	Edge(const VertexHandle& v);

	VertexHandle getVertex() const;

	EdgeHandle getPrev() const;

	EdgeHandle getNext() const;

	EdgeHandle getDual() const;

	FacetHandle getFacet() const;

	bool operator==(const Edge& e) const;

	bool operator!=(const Edge& e) const;

	friend class Vertex;

	friend class Facet;

	friend class SurfaceTopology;
};

class Facet
{
	/// The link to the incident edges counter-clockwisely
	/// oriented
	EdgeHandle edges[3];

	/// The area of this facet
	double area;

	/// The normal vector of this facet
	Vector normal;

  public:

	Facet();

	Facet(const EdgeHandle& e1, const EdgeHandle& e2, const EdgeHandle& e3);

	EdgeHandle* getEdges();

	double getArea() const;

	void update();

	Triangle getTriangle();

	bool operator==(const Facet& f) const;

	bool operator!=(const Facet& f) const;

	friend class Vertex;

	friend class Edge;

	friend class SurfaceTopology;
};

}

#endif /*VERTEXEDGEFACET_HPP_*/
