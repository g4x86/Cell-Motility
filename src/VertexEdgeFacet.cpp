#include <cassert>
#include <VertexEdgeFacet.hpp>
#include <FilamentBranch.hpp>
#include <Triangle.hpp>
#include <algorithms.hpp>

namespace motility
{

// Vertex definitions

Vertex::Vertex() {}

Vertex::Vertex(FilamentBranchHandle branch)
{
	//
	// Everytime when a new filament is added to cell surface as
	// a vertex, this filament and corresponding vertex must be
	// bi-directionally linked to each other.
	//
	// 1) The vertex holds a FilamentBranchHandle-type pointer to the
	// filament and the tree where this filament is. This operation
	// is carried out within Vertex constructor.
	//
	// 2) The filament also holds a VertexHandle-type pointer to
	// the vertex. This operation is carried out by the member
	// function 'addVertex' of Membrane when this vertex is created
	// and added to cell surface.
	//
	filament = branch;
}

FilamentBranchHandle Vertex::getFilament()
{
	return filament;
}

void Vertex::setFilament(FilamentBranchHandle fbh)
{
	filament = fbh;
}

EdgeHandles& Vertex::getEdges()
{
	return edges;
}

CartesianCoordinate Vertex::getLocation() const
{
	return filament->getTailEndLocation();
}

bool Vertex::operator==(const Vertex& v) const
{
	return (filament == v.filament);
}

bool Vertex::operator!=(const Vertex& v) const
{
	return (filament != v.filament);
}

// Edge definitions

Edge::Edge()
{
	vertex = vertex_handle_null;
	prev = edge_handle_null;
	next = edge_handle_null;
	dual = edge_handle_null;
	facet = facet_handle_null;
}

Edge::Edge(const VertexHandle& v)
{
	vertex = v;
	prev = edge_handle_null;
	next = edge_handle_null;
	dual = edge_handle_null;
	facet = facet_handle_null;
}

VertexHandle Edge::getVertex() const
{
	return vertex;
}

EdgeHandle Edge::getPrev() const
{
	return prev;
}

EdgeHandle Edge::getNext() const
{
	return next;
}

EdgeHandle Edge::getDual() const
{
	return dual;
}

FacetHandle Edge::getFacet() const
{
	return facet;
}

bool Edge::operator==(const Edge& e) const
{
	return (vertex == e.vertex);
}

bool Edge::operator!=(const Edge& e) const
{
	return (vertex != e.vertex);
}

// Facet definitions

Facet::Facet() {}

Facet::Facet(const EdgeHandle& e1, const EdgeHandle& e2, const EdgeHandle& e3)
{
	edges[0] = e1;
	edges[1] = e2;
	edges[2] = e3;
	update();
}

EdgeHandle* Facet::getEdges()
{
	return edges;
}

double Facet::getArea() const
{
	return area;
}

void Facet::update()
{
	Triangle t(edges[0]->vertex->getLocation(), edges[1]->vertex->getLocation(), edges[2]->vertex->getLocation());
	area = t.getArea();
	normal = t.getNormal();
}

Triangle Facet::getTriangle()
{
	return Triangle(edges[0]->vertex->getLocation(), edges[1]->vertex->getLocation(), edges[2]->vertex->getLocation());
}

bool Facet::operator==(const Facet& f) const
{
	return (edges[0] == f.edges[0] && edges[1] == f.edges[1] && edges[2] == f.edges[2]);
}

bool Facet::operator!=(const Facet& f) const
{
	return (edges[0] != f.edges[0] || edges[1] != f.edges[1] || edges[2] != f.edges[2]);
}

}
