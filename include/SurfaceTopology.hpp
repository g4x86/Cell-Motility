#ifndef SURFACETOPOLOGY_HPP_
#define SURFACETOPOLOGY_HPP_

#include <string>
#include <typedefs.hpp>
#include <Coordinate.hpp>
#include <Line.hpp>
#include <VertexEdgeFacet.hpp>

namespace motility
{

/// SurfaceTopology class describes the surface geometry of cell membrane.
///
///
/// Purpose
///
/// The data structure SurfaceTopology deals with the geometry of cell surface
/// during the dynamic re-organization of the actin filament networks in
/// cytosol. In order to update surface geometry according to the dynamic
/// changes of filaments, SurfaceTopology needs a Vertex-Edge-Facet data structure
/// which links cell surface geometry to the underlying actin filament 
/// networks.
///
/// And also this class includes the biophysical and biochemical properties
/// of cell membrane to support membrane-based biological processes.
///
///
/// Logical design
///
/// When SurfaceTopology needs to adjust its geometry through its various
/// member functions, it does not deal with underlying filament directly,
/// but only deals with geometric traits, such as vertex, edges and facets.
/// These geometrical traits are linked to actual actin filaments through
/// the FilamentBranchHandle contained in Vertex data structure.
///
///
/// Update local surface 
///
/// During the growing and the branching processes of actin filaments,
/// the relative distances between the growing ends of actin filaments
/// are changed constantly, such that the geometry of the composited
/// triangular facets are also changed constantly. These persistant
/// changes requires the necessary update of appropriate triangulation
/// because a good triangulation of a set of filament ends at current
/// time may become very bad next time if the locations of filament
/// ends change dramatically.
///
/// Therefore a crucial step to maintain the smoothness and the
/// correctness of membrane geometry is to update the triangulation of
/// a local set of filament ends. 3D Delaunay triangulation is one of
/// the possible methods.
///
/// An important assumption of updating the local geometry of cell
/// surface is that only ONE vertex is updated at one time, e.g. a
/// new vertex is inserted, an existing vertex is removed, or the
/// location of an existing vertex is changed.
///
///
/// Thumbrules
///
/// 1) Usually the member functions addVertex(), removeVertex() and
/// updateLocalSurface() are used by the updating function during
/// the reorganization of filament networks. Other member functions,
/// such as addFacet() and removeFacet() are not directly used by the
/// updating function. Instead they will be used internally in this
/// class based on the triangulation rules.
///
/// 2) The capped and detached filament will be removed from cell surface
/// and filament tree.
///
/// 3) Sometimes the member functions, e.g. addFacet() and removeFacet()
/// are used directly by the external function to create cell surface
/// geometry. For example, initialize_membrane_filaments() function uses
/// these member functions directly in order to create an initial sphere
/// cell surface.
///
/// 4) To get the properties of any vertex, edge and facet from this
/// class, the 'Handle'-type arguments are preferred, e.g. VertexHandle,
/// EdgeHandle and FacetHandle.
class SurfaceTopology
{
  private:
  
	/// The pool of the vertices representing the end of actin filaments.
	Vertices vertices;

	/// The pool of the edges connecting all vertices.
	Edges edges;

	/// The pool of the facets comprising cell membrane surface.
	Facets facets;

	/// The number of uncapped actin filaments.
	size_t n_uncapped_vertex;

	/// The total surface area of cell membrane.
	double area;

  private:
  
	/// This function inserts a vertex into a facet.
	///
	/// \param vertex_handle the handle of a vertex.
	/// \param facet_handle the handle of a facet into which the
	/// vertex is added.
	/// \return No returned value.
	void insertVerex(VertexHandle vertex_handle, FacetHandle facet_handle);

	/// This function updates the properties of neighboring facets.
	///
	/// \param vh the handle of a vertex on membrane surface.
	/// \return No returned value.
	void updateNeighboringFacets(VertexHandle vertex_handle);

  public:

	/// SurfaceTopology constructor function.
	SurfaceTopology();

	virtual ~SurfaceTopology();

	/// This function returns the handle of vertex pool on membrane surface.
	Vertices& getVertices();

	/// This function returns the handle of edge pool on membrane surface.
	Edges& getEdges();

	/// This function returns the handle of facet pool on membrane surface.
	Facets& getFacets();

	/// This function returns the number of vertices on membrane surface.
	size_t getVertexSize() const;

	/// This function returns the number of edges on membrane surface.
	size_t getEdgeSize() const;

	/// This function returns the number of facets on membrane surface.
	size_t getFacetSize() const;

	/// This function returns the number of uncapped actin filaments.
	size_t getVolatileVertexSize() const;

	/// This function returns the area of membrane surface.
	double getArea() const;

	/// This function returns the handle of the first vertex in vertex pool.
	VertexHandle getFirstVertexHandle();

	/// This function returns the handle of the last vertex in vertex pool.
	VertexHandle getLastVertexHandle();

	/// This function returns the handle of the first edge in edge pool.
	EdgeHandle getFirstEdgeHandle();

	/// This function returns the handle of the last edge in edge pool.
	EdgeHandle getLastEdgeHandle();

	/// This function returns the handle of the first facet in facet pool.
	FacetHandle getFirstFacetHandle();

	/// This function returns the handle of the last facet in facet pool.
	FacetHandle getLastFacetHandle();

	/// This function returns the position of a vertex in vertex pool.
	size_t getVertexHandlePosition(VertexHandle vh);

	/// This function returns the position of a edge in edge pool.
	size_t getEdgeHandlePosition(EdgeHandle eh);

	/// This function returns the position of a facet in facet pool.
	size_t getFacetHandlePosition(FacetHandle fh);

	/// This function updates the local surface around a dynamic vertex.
	///
	/// \param vertex_handle the handle of a dynamic vertex.
	/// \return The handles of all affected neighboring vertices.
	VertexHandles updateLocalSurface(VertexHandle vertex_handle);

	/// This function updates the local surface of a new vertex.
	///
	/// \param child_vertex_handle the handle of a new vertex.
	/// \param start_facet_handle the handle of the facet from
	/// which this algorithm starts looking for an appropriate
	/// intersecting facet.
	/// \return The handles of all affected neighboring vertices.
	VertexHandles updateLocalSurface(VertexHandle child_vertex_handle, FacetHandle start_facet_handle);

	/// This function perturbs the local surface of a vertex.
	///
	/// \param vertex_handle the handle of a given vertex.
	/// \return the execution status.
	bool perturbLocalSurface(VertexHandle vertex_handle);

	/// This function restores the local surface of a vertex.
	///
	/// \param vertex_handle the handle of a given vertex.
	/// \return the execution status.
	bool restoreLocalSurface(VertexHandle vertex_handle);

	/// This function calculates the centered directional local area around a vertex.
	///
	/// \param vertex_handle the vertex handle of a filament.
	/// \return The centered directional area of local surface around a filament.
	Vector computeCenteredDirectionalAreaOfLocalSurface(VertexHandle vertex_handle);


	/// This function calculates the centered directional local area around a point.
	///
	/// \param v the location of the tail end of a filament.
	/// \param facet_handle the facet handle which a filament intersects with.
	/// \return The centered directional area of local surface around a filament.
	Vector computeCenteredDirectionalAreaOfLocalSurface(CartesianCoordinate v, FacetHandle facet_handle);

	/// This function calculates the integral of surface curvature.
	///
	/// \param vertex_handles the handle of a list of vertices.
	/// \return The curvature integration of the local surface
	/// comprised by vertex_handles.
	double computeIntegralOfSurfaceCurvature(VertexHandles vertex_handles);

	/// This function calculates the area of local surface.
	///
	/// \param vertex_handle the handle of a given vertex.
	/// \return The area of local surface.
	double computeAreaOfLocalSurface(VertexHandle vertex_handle);

	/// This function triangulates the 3D polygon comprised by the given set of vertices.
	///
	/// \param polygon_vertices the handles of a set of counter-
	/// clockwisely oriented vertices that comprise the polygonal
	/// hole.
	/// \return A boolean status indicating whether or not the
	/// triangulation process succeeds.
	bool triangulatePolygonSurface(VertexHandles polygon_vertices);

	/// This function adds a given vertex into the vertex pool of membrane geometry.
	///
	/// \param vertex the vertex needed to be added into membrane geometry.
	/// \return The handle of the added vertex.
	VertexHandle addVertex(const Vertex& vertex);

	/// This function adds a new facet defined by three vertices into the facet pool of membrane geometry.
	///
	/// \param vh1 the handle of the first vertex.
	/// \param vh2 the handle of the second vertex.
	/// \param vh3 the handle of the third vertex.
	/// \return A boolean status indicating whether or not the
	/// facet specified by three given vertices is added into
	/// cell membrane surface.
	bool addFacet(VertexHandle vh1, VertexHandle vh2, VertexHandle vh3);

	/// This function removes a given facet from membrane geometry. 
	///
	/// \param fh the handle of a facet on membrane surface.
	/// \return No returned value.
	void removeFacet(FacetHandle fh);

	/// This function removes a given vertex from membrane geometry.
	///
	/// \param vh the handle of a vertex on membrane surface.
	/// \return The handles of all affected neighboring vertices.
	VertexHandles removeVertex(VertexHandle vh);

	/// This function updates the composite properties of cell membrane.
	///
	/// \param vh the handle of a vertex on membrane surface.
	/// \param pos_change a boolean variable indicating whether or not
	/// the localtion of the given vertex changes.
	/// \param cap_change a boolean variable indicating whether or not
	/// the capping status of corresponding actin filament changes.
	/// \return No returned value.
	void updateCompositeProperties(VertexHandle vh, bool pos_change, bool cap_change);

	/// This function returns the closest neighboring vertices around a vertex.
	///
	/// \param vertex_handle the handle of a vertex on membrane surface.
	/// \return The handles of the closest neighboring vertices.
	VertexHandles getNeighboringVertexHandles(VertexHandle vertex_handle) const;

	/// This function returns the closest local facets of a vertex.
	///
	/// \param vertex_handle the handle of a vertex on membrane surface.
	/// \return The handles of the closest local facets.
	FacetHandles getNeighboringFacetHandles(VertexHandle vertex_handle) const;

	/// This function returns the closest local facets of a facet.
	///
	/// \param facet_handle the handle of a facet on membrane surface.
	/// \return The handles of the closest local facets.
	FacetHandles getNeighboringFacetHandles(FacetHandle facet_handle) const;

	/// This function returns the local facets around a vertex.
	///
	/// \param vertex_handle the handle of a vertex on membrane surface.
	/// \return The handles of the closest neighboring facets.
	FacetHandles getLocalFacetHandles(VertexHandle vertex_handle) const;

	/// This function calculates the orientation of a child filament
	///
	/// \param b the location of the branching site on mother filament.
	/// \param v0 the location of the barbed end of mother filament.
	/// \param t the triangle of the facet that child filament is oriented
	/// towards.
	/// \return The orientation of child filament.
	Orientation computeBranchingOrientation(const CartesianCoordinate& bs, double alpha, const CartesianCoordinate& tip, Triangle& t);

	/// This function determins whether a branching reaction is allowed for an actin filament.
	///
	/// \param branch an actin filament.
	/// \return Whether or not a branching reaction can occur on this
	/// actin filament.
	bool isBranchingAllowed(FilamentBranch& branch);

	/// This function determines whether or not an actin filament intersects with a facet or its local facets.
	///
	/// \param branch an actin filament.
	/// \param start_facet_handle the handle of a facet on membrane
	/// surface.
	/// \return The handle of a facet with which the given actin
	/// filament intersects. If no such facet is found, then 0 is
	/// returned.
	FacetHandle getIntersectingFacetHandle(FilamentBranch& branch, FacetHandle start_facet_handle);

	/// This function calculates the orientation of extracellular signaling
	/// region which will be used as the perfered branching orientation for
	/// child filament.
	///
	/// \param loc the location in space where its extracellular orientation
	/// is calculated.
	/// \return The orientation pointing to extracellular signaling region.
	Orientation computeExtraCellularOrientation(const CartesianCoordinate& loc);

	/// This function calculates the angle between filament growing
	/// direction and the radial direction.
	double computeDeviationAngleOfFilamentGrowth(FilamentBranch& branch);

	double computeDeviationAngleOfFilamentGrowth(FilamentBranch& mother_branch, const Orientation& child_branch_orient);

	/// This function searches the facet which a child branch represented
	/// by a line will be fused with.
	///
	/// \param child_branch_line the line representing a child actin filament.
	/// \param mother_branch_vertex the vertex of the barbed end of a mother
	/// actin filament.
	/// \param child_fusion_facet the found fusion facet.
	/// \return The status of search.
	bool searchFusionFacetForChildBranch(Line& child_branch_line, VertexHandle mother_branch_vertex, FacetHandle& child_fusion_facet);

	/// This function writes the surface geometry of cell membrane into an OFF-format file.
	///
	/// \param output the output stream.
	/// \param delim the delimiter.
	/// \return No value is returned.
	void exportGeometry(std::ofstream& output, char delim = ' ');
};

}

#endif /*SURFACETOPOLOGY_HPP_*/
