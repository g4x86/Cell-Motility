#include <cassert>
#include <cfloat>
#include <fstream>
#include <sstream>
#include <SurfaceTopology.hpp>
#include <constants.hpp>
#include <FilamentBranch.hpp>
#include <Vector.hpp>
#include <Line.hpp>
#include <Triangle.hpp>
#include <ParameterTable.hpp>
#include <algorithms.hpp>

namespace motility
{

SurfaceTopology::SurfaceTopology()
{
	n_uncapped_vertex = 0;
	area = 0;
}

SurfaceTopology::~SurfaceTopology() {}

Vertices& SurfaceTopology::getVertices()
{
	return vertices;
}

Edges& SurfaceTopology::getEdges()
{
	return edges;
}

Facets& SurfaceTopology::getFacets()
{
	return facets;
}

size_t SurfaceTopology::getVertexSize() const
{
	return vertices.size();
}

size_t SurfaceTopology::getEdgeSize() const
{
	return edges.size();
}

size_t SurfaceTopology::getFacetSize() const
{
	return facets.size();
}

size_t SurfaceTopology::getVolatileVertexSize() const
{
	return n_uncapped_vertex;
}

double SurfaceTopology::getArea() const
{
	return area;
}

VertexHandle SurfaceTopology::getFirstVertexHandle()
{
	return vertices.begin();
}

VertexHandle SurfaceTopology::getLastVertexHandle()
{
	return --(vertices.end());
}

EdgeHandle SurfaceTopology::getFirstEdgeHandle()
{
	return edges.begin();
}

EdgeHandle SurfaceTopology::getLastEdgeHandle()
{
	return --(edges.end());
}

FacetHandle SurfaceTopology::getFirstFacetHandle()
{
	return facets.begin();
}

FacetHandle SurfaceTopology::getLastFacetHandle()
{
	return --(facets.end());
}

size_t SurfaceTopology::getVertexHandlePosition(VertexHandle vh)
{
	return distance(vertices.begin(), vh);
}

size_t SurfaceTopology::getEdgeHandlePosition(EdgeHandle eh)
{
	return distance(edges.begin(), eh);
}

size_t SurfaceTopology::getFacetHandlePosition(FacetHandle fh)
{
	return distance(facets.begin(), fh);
}

VertexHandles SurfaceTopology::updateLocalSurface(VertexHandle vertex_handle)
{
	// Step 1. Check if current vertex is located at the center of
	// the polygon comprised by its neighboring vertices.
	//
	// Step 2. If not, find where current vertex is close to.
	//
	// Step 3. Check if the filament intersects with any facets if
	// the vertex is moved to the new location.
	//
	// Step 4. If no intersection occurs, then add current vertex to
	// the neighboring facet which current vertex is close to.
	//
	// Step 5. Disconnect the edge between current vertex and its
	// original neighboring vertices which current vertex is far
	// from.

	VertexHandles affected_vertices = getNeighboringVertexHandles(vertex_handle);
	CartesianCoordinate vertex_location = vertex_handle->getLocation();
	EdgeHandles vertex_edges = vertex_handle->edges;
	EdgeHandles::reverse_iterator ehh = vertex_edges.rbegin();
	while(ehh != vertex_edges.rend())
	{
		bool updated_flag = false;
		EdgeHandle curr_edge = (*ehh);
		EdgeHandle next_edge = curr_edge->next;
		EdgeHandle prev_edge = curr_edge->prev;
		EdgeHandle prev_dual = prev_edge->dual;
		FacetHandle curr_facet = curr_edge->facet;
		// Step #1: Check whether there is a neighboring facet.
		if(prev_dual != static_cast<EdgeHandle>(nullptr))
		{
			// Step #2: Check whether V1 is inside the positive side of
			// the neighboring facet.
			FacetHandle neighboring_facet = prev_dual->facet;
			Triangle t = neighboring_facet->getTriangle();
			if(isInsideTriangularBox(vertex_location, t, DBL_INF_POSITIVE))
			{
				// Step #3: Find four vertices V1, V2, V3 and V4.
				VertexHandle vh1 = vertex_handle;
				VertexHandle vh2 = next_edge->vertex;
				VertexHandle vh3 = prev_edge->vertex;
				VertexHandle vh4 = prev_dual->next->vertex;
				// Step #4: Check whether the vertex V1 and V4 are connected.
				// Step 4.1 Check if there is an edge from v4 to v1.
				bool v1_v4_connection_flag = isContained<VertexHandle>(getNeighboringVertexHandles(vh1), vh4);
				// Step 4.2 Check if there is an edge from v1 to v4.
				if(!v1_v4_connection_flag) v1_v4_connection_flag = isContained<VertexHandle>(getNeighboringVertexHandles(vh4), vh1);
				if(!v1_v4_connection_flag)
				{
					// Step #5: If the vertex V1 and V4 are not connected,
					// then update the local surface comprised by the facet
					// V1-V2-V3 and V4-V3-V2, and add the vertex V2, V3 and
					// v4 into affected_vertices and affected_vertices_2.
					++ehh;
					removeFacet(curr_facet);
					removeFacet(neighboring_facet);
					addFacet(vh2, vh4, vh1);
					addFacet(vh3, vh1, vh4);
					affected_vertices.push_back(vh4);
					updated_flag = true;
				}
			}
		}
		if(!updated_flag) ++ehh;
	}
	// Do NOT include the current filament into the list of affected
	// vertices because it must remain valid through the later update
	// of reaction interactions.
	return affected_vertices;
}

VertexHandles SurfaceTopology::updateLocalSurface(VertexHandle child_vertex_handle, FacetHandle start_facet_handle)
{
	VertexHandles affected_vertices;
	// Step #1: Check whether the child filament intersects with the given
	// starting facet.
	FacetHandle intersecting_facet_handle = getIntersectingFacetHandle(*(child_vertex_handle->filament), start_facet_handle);
	if(intersecting_facet_handle != static_cast<FacetHandle>(nullptr))
	{
		// Step #2: If the child filament intersects with the given starting
		// facet, check whether the child vertex is located in the positive
		// region of the neighboring facets of the intersecting facet:
		CartesianCoordinate child_vertex_location = child_vertex_handle->getLocation();
		Triangle t = intersecting_facet_handle->getTriangle();
		// The distance between the child vertex and the starting facet.
		double intersecting_distance = distance(child_vertex_location, t);
		FacetHandles local_facet_handles = getNeighboringFacetHandles(intersecting_facet_handle);
		// Remember that all candidate facets are the neighbors of the
		// given starting facet with which the child filament intersects.
		FacetHandles candidate_facet_handles;
		for(FacetHandleHandle fhh = local_facet_handles.begin(); fhh != local_facet_handles.end(); ++fhh)
		{
			Triangle t = (*fhh)->getTriangle();
			// 'candidate_facet_handles' contains all the facets within
			// which the child vertex is projected.
			if(isInsideTriangularBox(child_vertex_location, t, DBL_INF_POSITIVE)) candidate_facet_handles.push_back(*fhh);
		}
		// Step #3: If the child vertex is located within the positive region
		// of the neighboring facets, find out the closest neighboring facet.
		FacetHandle closest_facet_handle = static_cast<FacetHandle>(nullptr);
		double closest_distance = DBL_INF_POSITIVE;
		if(!candidate_facet_handles.empty())
		{
			// If there exist neighboring facets that contain the
			// projection of the child vertex, then find the facet
			// closest to the child vertex.
			for(FacetHandleHandle fhh = candidate_facet_handles.begin(); fhh != candidate_facet_handles.end(); ++fhh)
			{
				Triangle t = (*fhh)->getTriangle();
				double dst = distance(child_vertex_location, t);
				if(dst < closest_distance)
				{
					closest_facet_handle = *fhh;
					closest_distance = dst;
				}
			}
		}
		if(closest_facet_handle != static_cast<FacetHandle>(nullptr) && closest_distance <= intersecting_distance)
		{
			// Step #4: If the closest facet is not further from the child vertex
			// than the intersecting facet, the intersecting facet and the closest
			// neighboring facet will be used to determine the methods of actual
			// re-triangulation.
			// Step #5: Check whether the intersecting facet and the closest neighboring
			// facet share the same edge.
			EdgeHandle* intersecting_facet_edges = intersecting_facet_handle->edges;
			EdgeHandle* closest_facet_edges = closest_facet_handle->edges;
			bool shared_edge_flag = false;
			bool shared_vertex_flag = false;
			size_t i, j;
			for(i = 0; i < 3; ++i)
			{
				for(j = 0; j < 3; ++j)
				{
					if(intersecting_facet_edges[i]->dual == closest_facet_edges[j])
					{
						shared_edge_flag = true;
						break;
					}
				}
				if(shared_edge_flag) break;
			}
			if(shared_edge_flag)
			{
				// Step #6: If the intersecting facet and the closest neighboring facet share
				// the same edge, re-triangulate the both facets together by:
				//
				// 1) Removing the shared edge.
				//
				// 2) Connecting the child vertex to all vertices of these two facets.
				VertexHandle shared_vertex_handle_1 = intersecting_facet_edges[i]->vertex;
				VertexHandle shared_vertex_handle_2 = intersecting_facet_edges[i]->dual->vertex;
				VertexHandle intersecting_vertex_handle = intersecting_facet_edges[i]->next->vertex;
				VertexHandle closest_vertex_handle = intersecting_facet_edges[i]->dual->next->vertex;
				removeFacet(intersecting_facet_handle);
				removeFacet(closest_facet_handle);
				assert(addFacet(child_vertex_handle, shared_vertex_handle_1, intersecting_vertex_handle));
				assert(addFacet(child_vertex_handle, intersecting_vertex_handle, shared_vertex_handle_2));
				assert(addFacet(child_vertex_handle, shared_vertex_handle_2, closest_vertex_handle));
				assert(addFacet(child_vertex_handle, closest_vertex_handle, shared_vertex_handle_1));
				affected_vertices.push_back(intersecting_vertex_handle);
				affected_vertices.push_back(closest_vertex_handle);
				affected_vertices.push_back(shared_vertex_handle_1);
				affected_vertices.push_back(shared_vertex_handle_2);
			}
			else
			{
				// Step #7: If the intersecting facet and the closest neighboring facet
				// do not share the same edge, check whether they share the same vertex.
				for(i = 0; i < 3; ++i)
				{
					for(j = 0; j < 3; ++j)
					{
						if(intersecting_facet_edges[i]->vertex == closest_facet_edges[j]->vertex)
						{
							shared_vertex_flag = true;
							break;
						}
					}
					if(shared_vertex_flag) break;
				}
				if(shared_vertex_flag)
				{
					// Step #8: If the intersecting facet and the closest neighboring facet
					// share the same vertex, check whether the neighboring vertices of the
					// shared vertex can be divided into two groups such that:
					//
					// 1) The vertices in one group is closer to the shared vertex.
					//
					// 2) The vertices in the other group is closer to the child vertex.
					//
					// 3) Two groups have only two vertices which can be connected to both
					// the child vertex and the shared vertex.
					VertexHandle shared_vertex_handle = intersecting_facet_edges[i]->vertex;
					EdgeHandles& shared_vertex_edges = shared_vertex_handle->edges;
					VertexHandles polygon_vertices;
					FacetHandles shared_vertex_facets;
					size_t n_shared_vertex_edge = shared_vertex_edges.size();
					// For the child vertex and the shared vertex, compare the distances
					// from them to the two vertices of each edge of the boundary polygon
					// comprised by the direct neighboring vertices of the shared vertex.
					// If the distance for the child vertex case is less than that for the
					// shared case, set corresponding element of 'distance_comparison' to
					// true. Otherwise set it to false.
					bool* distance_comparison = 0;
					distance_comparison = new bool [n_shared_vertex_edge];
					size_t k = 0;
					for(EdgeHandles::reverse_iterator ehh = shared_vertex_edges.rbegin(); ehh != shared_vertex_edges.rend(); ++ehh)
					{
						polygon_vertices.push_back((*ehh)->prev->vertex);
						shared_vertex_facets.push_back((*ehh)->facet);
						double child_vertex_distance = distance(child_vertex_handle->getLocation(), (*ehh)->next->vertex->getLocation()) + distance(child_vertex_handle->getLocation(), (*ehh)->prev->vertex->getLocation());
						double shared_vertex_distance = distance(shared_vertex_handle->getLocation(), (*ehh)->next->vertex->getLocation()) + distance(shared_vertex_handle->getLocation(), (*ehh)->prev->vertex->getLocation());
						if(child_vertex_distance < shared_vertex_distance) distance_comparison[k] = true;
						else distance_comparison[k] = false;
						++k;
					}
					// Determine where the change of distance comparison occurs.
					size_t n_distance_comparison_change = 0;
					bool current_compare = false;
					for(k = 0; k < n_shared_vertex_edge; ++k)
					{
						if(k == 0) current_compare = distance_comparison[k];
						else
						{
							if(distance_comparison[k] != current_compare)
							{
								current_compare = distance_comparison[k];
								++n_distance_comparison_change;
							}
						}
						if(k == n_shared_vertex_edge - 1)
						{
							if(distance_comparison[0] != current_compare)
							{
								current_compare = distance_comparison[0];
								++n_distance_comparison_change;
							}
						}
					}
					if(n_distance_comparison_change == 2)
					{
						// Step #9: If the neighboring vertices of the shared vertex can be divided
						// into two groups, re-triangulate the local surface around the shared vertex
						// by:
						//
						// 1) Connecting two groups of vertices to the shared vertex and the child
						// vertex respectively.
						//
						// 2) Connecting the two common vertices of two groups to both vertices.
						//
						// Now remove the existing facets around the shared vertex
						for(FacetHandleHandle fhh = shared_vertex_facets.begin(); fhh != shared_vertex_facets.end(); ++fhh) removeFacet(*fhh);
						// Then connect the child vertex and the shared vertex to the
						// corresponding vertices of the boundary polygon.
						k = 0;
						for(VertexHandleHandle vhh = polygon_vertices.begin(); vhh != polygon_vertices.end(); ++vhh)
						{
							size_t l;
							if(k < n_shared_vertex_edge - 1) l = k + 1;
							else l = 0;
							VertexHandleHandle vhh3 = vhh;
							VertexHandleHandle vhh2 = vhh;
							if(vhh2 == polygon_vertices.begin()) --vhh2;
							--vhh2;
							if(distance_comparison[k])
							{
								assert(addFacet(child_vertex_handle, *vhh2, *vhh3));
								if(distance_comparison[l] != distance_comparison[k]) assert(addFacet(child_vertex_handle, *vhh3, shared_vertex_handle));
							}
							else
							{
								assert(addFacet(shared_vertex_handle, *vhh2, *vhh3));
								if(distance_comparison[l] != distance_comparison[k]) assert(addFacet(shared_vertex_handle, *vhh3, child_vertex_handle));
							}
							++k;
						}
						if(distance_comparison != 0) delete [] distance_comparison;
						affected_vertices = polygon_vertices;
						affected_vertices.push_front(shared_vertex_handle);
					}
					else
					{
						// Step #10: If the neighboring vertices of the shared vertex can not be
						// divided into two groups, then simply insert the child vertex into the
						// intersecting facet.
						insertVerex(child_vertex_handle, intersecting_facet_handle);
						affected_vertices = getNeighboringVertexHandles(child_vertex_handle);
					}
				}
			}
			// Make sure that either of above two cases occur.
			assert(shared_edge_flag || shared_vertex_flag);
		}
		else
		{
			// Step #11: If the child vertex is not located within the positive region
			// of any neighboring facets of the intersecting facet, or the closest facet
			// is further away from the child vertex than the intersecting facet, then
			// simply insert the child vertex into the intersecting facet.
			insertVerex(child_vertex_handle, intersecting_facet_handle);
			affected_vertices = getNeighboringVertexHandles(child_vertex_handle);
		}
	}
	else
	{
		// Step #12: If the child filament does not intersect with the given
		// starting facet, then simply insert the child vertex into the starting
		// facet.
		insertVerex(child_vertex_handle, start_facet_handle);
		affected_vertices = getNeighboringVertexHandles(child_vertex_handle);
	}
	// This child filament must be added into the list of affected
	// filaments because it is affected by the branching reaction
	// occuriing on the mother filament, and the mother filament
	// may also be included in the list based on the addition of
	// the child vertex to the local surface of the mother filament.
	affected_vertices.push_front(child_vertex_handle);
	return affected_vertices;
}

double SurfaceTopology::computeIntegralOfSurfaceCurvature(VertexHandles vertex_handles)
{
	///
	/// The integral of surface curvature is calculated according to the
	/// method proposed by Itzykson C. 1986 which computes the integral
	/// through the summation of discretized Laplacian bending energy over
	/// a triangulated surface area.
	///
	double curvature_integral = 0;
	for(VertexHandleHandle vhh = vertex_handles.begin(); vhh != vertex_handles.end(); ++vhh)
	{
		Vector Ri(CartesianCoordinate(0, 0, 0), (*vhh)->getLocation());
		double sigma = 0;
		Vector covariant_laplacian;
		EdgeHandles& ehs = (*vhh)->edges;
		for(EdgeHandleHandle ehh = ehs.begin(); ehh != ehs.end(); ++ehh)
		{
			//
			//                    i
			//                    O
			//                   *|*
			//                 *  | *
			//               *    |   *
			//             *      |    *
			//           *        |      *
			//         *          |       *
			//       *            |         *
			//   k1 O) theta1     |  theta2 (O k2
			//       *            |         *
			//         *          |       *
			//           *        |      *
			//             *      |    *
			//               *    |   *
			//                 *  | *
			//                   *|*
			//                    O
			//                    j
			//
			VertexHandle vhj = (*ehh)->prev->vertex;
			Vector Rj(CartesianCoordinate(0, 0, 0), vhj->getLocation());
			VertexHandle vhk1 = (*ehh)->next->vertex;
			VertexHandle vhk2 = (*ehh)->dual->next->vertex;
			Vector Vk1i(vhk1->getLocation(), (*vhh)->getLocation());
			Vector Vk1j(vhk1->getLocation(), vhj->getLocation());
			double theta1 = angle(Vk1i, Vk1j);
			Vector Vk2i(vhk2->getLocation(), (*vhh)->getLocation());
			Vector Vk2j(vhk2->getLocation(), vhj->getLocation());
			double theta2 = angle(Vk2i, Vk2j);
			Vector Rij = Ri - Rj;
			double lij = abs(Rij);
			double theta1_cot = cot(theta1); 
			double theta2_cot = cot(theta2);
			assert(::finite(theta1_cot) && ::finite(theta2_cot));
			covariant_laplacian += (Rij * (0.5 * (theta1_cot + theta2_cot)));
			sigma += (0.5 * (theta1_cot + theta2_cot) * lij * lij);
		}
		double covariant_laplacian_abs = abs(covariant_laplacian);
		sigma /= 4;
		curvature_integral += (covariant_laplacian_abs * covariant_laplacian_abs / sigma);
	}
	return curvature_integral;
}

bool SurfaceTopology::perturbLocalSurface(VertexHandle vertex_handle)
{
	///
	/// The local surface around the given vertex is perturbed by
	/// adding one actin monomer to corresponding filament.
	///
	ParameterTable::Table& param_table = ParameterTable::instance();
	double actin_diameter = strtod(param_table[std::string("actin_diameter")]);
	FilamentBranchHandle branch_handle = vertex_handle->getFilament();
	bool result = branch_handle->perturbTailEndLocation(Actin("ATP", actin_diameter / 2));
	return result;
}

bool SurfaceTopology::restoreLocalSurface(VertexHandle vertex_handle)
{
	///
	/// The local surface around the given vertex is restored by
	/// removing one actin monomer from the filament.
	///
	FilamentBranchHandle branch_handle = vertex_handle->getFilament();
	bool result = branch_handle->restoreTailEndLocation();
	return result;
}

double SurfaceTopology::computeAreaOfLocalSurface(VertexHandle vertex_handle)
{
	double surface_area = 0;
	FacetHandles fhs = 	getNeighboringFacetHandles(vertex_handle);
	for(FacetHandleHandle fhh = fhs.begin(); fhh != fhs.end(); ++fhh) surface_area += (*fhh)->area;
	return surface_area;
}

Vector SurfaceTopology::computeCenteredDirectionalAreaOfLocalSurface(VertexHandle vertex_handle)
{
	///
	/// The centered directional area of a neighboring triangular
	/// facet of a given vertex is defined as a vector whose
	/// direction is the normal vector from the given vertex to
	/// the center of this facet and whose magnitude is the area
	/// of this triangular facet. For example, R(V->C1) is the
	/// vector pointing from the vertex V and the center C1.
	///
	///
	///                       V
	///                       *
	///                     * /* \*
	///                   *  / *  \  *
	///                 *   /   *  \    * .
	///               *    /    *   *C2 *
	///             *     *      *    *
	///           *       C1     *  *
	///         ******************* 
	///
	/// Therefore the total centered directional area of the local
	/// facets is the summation of each centered diretional area
	/// of all these facets. So the centered directional area of
	/// the vertex V is calculated by summing over all neighboring
	/// facets.
	///
	/// This algorithm needs to be double checked!
	///
	Vector total_direct_area;
	EdgeHandles& ehs = vertex_handle->edges;
	for(EdgeHandleHandle ehh = ehs.begin(); ehh != ehs.end(); ++ehh)
	{
		FacetHandle& fh = (*ehh)->facet;
		Triangle t = fh->getTriangle();
		Vector direct_area(vertex_handle->getLocation(), t.getCenter());
		direct_area.setMag(fh->area);
		total_direct_area += direct_area;
	}
	return total_direct_area;
}

Vector SurfaceTopology::computeCenteredDirectionalAreaOfLocalSurface(CartesianCoordinate v, FacetHandle facet_handle)
{
	Vector total_direct_area;
	Triangle facet_triangle = facet_handle->getTriangle();
	const CartesianCoordinate* vs = facet_triangle.getVertices();
	for(size_t i = 0; i < 3; ++i)
	{
		size_t j;
		if(i < 2) j = i + 1;
		else j = 0;
		Triangle t(v, vs[i], vs[j]);
		Vector direct_area(v, t.getCenter());
		direct_area.setMag(t.getArea());
		total_direct_area += direct_area;
	}
	return total_direct_area;
}

bool SurfaceTopology::triangulatePolygonSurface(VertexHandles polygon_vertices)
{
	///
	/// After a vertex is removed from cell surface, a 3D polygonal
	/// hole is created. This hole must be triangulated into smaller
	/// triangular facets that make up a new local surface.
	///
	///
	/// Rationale
	///
	/// When a vertex is removed, a 3D polygon hole is left on cell surface.
	/// This function re-triangulates this 3D polygon into a set of triangular
	/// facets which represent the smooth surface enclosed by the 3D polygon.
	///
	///
	/// Prerequisite
	///
	/// The requirements for the eligibility of the given polygon to be
	/// triangulated:
	///
	/// 1) The given polygon must have three vertices at least.
	///
	/// 2) All polygon edges must be boundary edges.
	///
	///
	/// Algorithm
	///
	/// The re-triangulation of a 3D polygon includes two parts: the traverse
	/// process and the decision-making process.
	///
	/// 1) The traverse process is an edge-based process that goes through all
	/// edges linking the neighboring vertices of the polygon and form smaller
	/// polygon. This process is outlined in the following steps:
	///
	///     a) Start from an edge of the given 3D polygon and traverse in a
	///     counter-clockwise direction.
	///
	///     b) Connect all qualified two neighboring edges of the polygon and
	///     create non-overlapping new facets.
	///
	///     c) Traverse back to the starting edge and the remaining boundary
	///     edges form a new smaller 3D polygon.
	///
	///     d) Repeat the re-triangulation process until the new 3D polygon
	///     becomes a triangle and then create a new facet from this triangle.
	///
	/// 2) The decision-make process is to select a qulified facet triangle that
	/// satisfies the following empirical requirements:
	///
	///     a) Triangle tends to be as regular as possible:
	///
	///     Choice #1 - compare side length (current implementation).
	///
	///     If a(i) (i=1,2,3) denotes the length of the three side of the
	///     triangle, then always select a triangle with the minimum value of:
	///
	///                       max{a(i)} - min{a(i)}
	///                      -----------------------
	///                             min{a(i)}
	///
	///     Choice #2 - use 2D Delaunay triangulation.
	///
	///     b) If n is the normal vector of the new triangular facet, and n1
	///     and n2 are the normal vectors of two side-neighboring triangular
	///     facets of the new facet, then the following requirements must be
	///     satisfied:
	///
	///                 n * n1 >= 0   and   n * n2 >= 0
	///
	///
	/// Attention
	///
	/// Since this function keeps changing the list of boundary vertices
	/// constantly, it is important to use the value-passing scheme to
	/// accept the input argument.
	///
	bool eligibility_flag = true;
	// Prerequisite #1
	if(polygon_vertices.size() < 3) eligibility_flag = false;
	else
	{
		for(VertexHandleHandle vhh = polygon_vertices.begin(); vhh != polygon_vertices.end(); ++vhh)
		{
			VertexHandleHandle vhh1 = vhh;
			if(++vhh1 == polygon_vertices.end()) ++vhh1;
			bool neighboring_edge_flag = false;
			for(EdgeHandleHandle ehh = (*vhh)->edges.begin(); ehh != (*vhh)->edges.end(); ++ehh)
			{
				VertexHandle begin_vertex_handle = (*ehh)->prev->vertex;
				// Prerequisite #2
				if(begin_vertex_handle == *vhh1 && (*ehh)->dual == static_cast<EdgeHandle>(nullptr)) neighboring_edge_flag = true;
			}
			if(!neighboring_edge_flag)
			{
				eligibility_flag = false;
				break;
			}
		}
	}
	if(eligibility_flag)
	{
		VertexHandleHandle vhh, vhh0, vhh1, vhh2, vhh3;
		VertexHandles new_polygon_vertices;
		bool re_triangulation_flag = false;
		while(polygon_vertices.size() >= 3)
		{
			bool traverse_stop_flag = false;
			vhh = polygon_vertices.begin();
			while(vhh != polygon_vertices.end())
			{
				// This variable indicates the type of candidate triangulation
				// facet to be selected:
				//
				// 0: No triangulation facet is added.
				// 1: Type-1 triangulation facet is added.
				// 2: Type-2 triangulation facet is added.
				size_t facet_addition_type = 0;
				vhh1 = vhh;
				vhh0 = vhh1;
				if(--vhh0 == polygon_vertices.end()) --vhh0;
				vhh2 = vhh1;
				if(++vhh2 == polygon_vertices.end()) ++vhh2;
				vhh3 = vhh2;
				if(++vhh3 == polygon_vertices.end()) ++vhh3;
				if(vhh3 != vhh0)
				{
					// The current polygon has at least four vertices.
					//
					// The following steps are taken to select an eligible triangle:
					//
					// 1) Construct two candidate facets from four neighboring vertices
					// on the boundary of the given 3D polygon:
					//
					//      (vhh0, vhh1, vhh2) and (vhh1, vhh2, vhh3).
					//
					// 2) Occupation eligibility
					//
					// a) If the facets contructed from these sets of vertices already exist,
					// say, the existing neighboring edges (vhh1 --> vhh0) and (vhh2 --> vhh1)
					// belong to the same existing facet, then the candidate facet (vhh0,
					// vhh1, vhh2) is invalid.
					//
					// b) If the edge connecting the non-neighboring vertices of the candidate
					// facet already exists, say the edge (vhh2 --> vhh0), then candidate facet
					// (vhh0, vhh1, vhh2) is invalid.
					//
					// This eligibility is required by the topology of 2-manifold surface,
					// which must be satisfied.
					//
					// 3) Orientation eligibility
					//
					// a) If there is only one candidate facet, say, (vhh0, vhh1, vhh2),
					// check its orientation relationships with its two side-neighboring
					// facets which have the following boundary edges:
					//
					//              edge (v1 --> v0) and edge (v2 --> v1)
					//
					// b) If there are two candidate facets, check their orientation
					// relationships with their three side-neighboring facets which have
					// the following boundary edges:
					//
					//     edge (v1 --> v0), edge (v2 --> v1) and edge (v3 --> v2)
					//
					// Each candidate facet is eligible only if its orientation is not
					// againt to both of its two side-neighboring facets.
					// edges:
					//
					// This eligibility is required by surface smoothness, which can be
					// compromised under certain circumstances.
					//
					// 4) If both candiate facets are orientationally eligible, compare
					// the regularity of corresponding triangles.
					//
					// 5) If neither of candiate facets is orientationally eligible, then
					// move forward conditionally along the polygon boundary.
					//
					//
					// Step #1
					// Corresponding candidate triangles are constructed after their
					// occupation eligibilities are examined.
					VertexHandle candidate_vertices[4];
					candidate_vertices[0] = *vhh0;
					candidate_vertices[1] = *vhh1;
					candidate_vertices[2] = *vhh2;
					candidate_vertices[3] = *vhh3;
					// Step #2
					EdgeHandle candidate_edges[3];
					for(size_t i = 0; i < 3; ++i)
					{
						for(EdgeHandleHandle ehh = candidate_vertices[i]->edges.begin(); ehh != candidate_vertices[i]->edges.end(); ++ehh)
						{
							if((*ehh)->prev->vertex == candidate_vertices[i + 1]) candidate_edges[i] = *ehh;
						}
					}
					bool non_occupation_flags[2] = {true, true};
					// Step #2 (a)
					for(size_t i = 0; i < 2; ++i)
					{
						if(candidate_edges[i]->facet != candidate_edges[i + 1]->facet) non_occupation_flags[i] = true;
						else non_occupation_flags[i] = false;
					}
					// Step #2 (b)
					for(size_t i = 0; i < 2; ++i)
					{
						for(EdgeHandleHandle ehh = candidate_vertices[i]->edges.begin(); ehh != candidate_vertices[i]->edges.end(); ++ehh)
						{
							if((*ehh)->prev->vertex == candidate_vertices[i + 2]) non_occupation_flags[i] = false;
						}
					}
					// Construct candidate triangles.
					Triangle candidate_triangles[2];
					for(size_t i = 0; i < 2; ++i)
					{
						if(non_occupation_flags[i]) candidate_triangles[i] = Triangle(candidate_vertices[i]->getLocation(), candidate_vertices[i + 1]->getLocation(), candidate_vertices[i + 2]->getLocation());
					}
					// Step #3
					FacetHandle neighboring_facets[3];
					for(size_t i = 0; i < 3; ++i)
					{
						for(EdgeHandleHandle ehh = candidate_vertices[i]->edges.begin(); ehh != candidate_vertices[i]->edges.end(); ++ehh)
						{
							if((*ehh)->prev->vertex == candidate_vertices[i + 1]) neighboring_facets[i] = (*ehh)->facet;
						}
					}
					double orientation_consistencies[2][2];
					for(size_t i = 0; i < 2; ++i)
					{
						for(size_t j = 0; j < 2; ++j)
						{
							if(non_occupation_flags[i]) orientation_consistencies[i][j] = dotProd(candidate_triangles[i].getNormal(), neighboring_facets[i + j]->getTriangle().getNormal());
							else orientation_consistencies[i][j] = -1;
						}
					}
					bool consistent_orientation_flags[2] = {true, true};
					for(size_t i = 0; i < 2; ++i)
					{
						for(size_t j = 0; j < 2; ++j)
						{
							if(orientation_consistencies[i][j] < 0)
							{
								consistent_orientation_flags[i] = false;
								break;
							}
						}
					}
					// Step #4
					//
					// Only the satisfying candidate facets are subject to the
					// calculation of triangular regularity. Two criteria are
					// used to make decision:
					//
					// 1) Occupation eligibility has the highest priority and
					// it must be satisfied.
					//
					// 2) If occupation eligibility is satisfied, then check
					// the orientation eligibility.
					double regularities[2];
					for(size_t i = 0; i < 2; ++i)
					{
						if(non_occupation_flags[i])
						{
							if(!re_triangulation_flag)
							{
								// In the normal triangulation of current polygon, the check for
								// surface smoothness is required.
								if(consistent_orientation_flags[i]) regularities[i] = candidate_triangles[i].getRegularity();
								else regularities[i] = 0;
							}
							else
							{
								// In the re-triangulation of current polygon, the check for
								// surface smoothness is not required.
								regularities[i] = candidate_triangles[i].getRegularity();
							}
						}
						else regularities[i] = 0;
					}
					double max_regularity = 0;
					for(size_t i = 0; i < 2; ++i)
					{
						if(regularities[i] > max_regularity)
						{
							max_regularity = regularities[i];
							facet_addition_type = i + 1;
						}
					}
					if(facet_addition_type > 0)
					{
						// Add the facet of corresponding type into cell surface.
						assert(addFacet(candidate_vertices[facet_addition_type - 1], candidate_vertices[facet_addition_type], candidate_vertices[facet_addition_type + 1]));
						if(candidate_vertices[facet_addition_type - 1] != new_polygon_vertices.front() && candidate_vertices[facet_addition_type - 1] != new_polygon_vertices.back()) new_polygon_vertices.push_back(candidate_vertices[facet_addition_type - 1]);
						if(candidate_vertices[facet_addition_type + 1] != new_polygon_vertices.front() && candidate_vertices[facet_addition_type + 1] != new_polygon_vertices.back()) new_polygon_vertices.push_back(candidate_vertices[facet_addition_type + 1]);
					}
					else
					{
						// Do not add any candidate facets.
						if(candidate_vertices[1] != new_polygon_vertices.front() && candidate_vertices[1] != new_polygon_vertices.back()) new_polygon_vertices.push_back(candidate_vertices[1]);
					}
					// Step #5
					for(size_t i = 0; i < facet_addition_type + 1; ++i)
					{
						if(++vhh == polygon_vertices.end()) ++vhh;
						if(new_polygon_vertices.front() == polygon_vertices.back())
						{
							// Finish traversing at the last second vertex of current polygon.
							VertexHandleHandle pv_last_second_vhh = polygon_vertices.end();
							--pv_last_second_vhh; --pv_last_second_vhh;
							if(vhh == pv_last_second_vhh)
							{
								if(*vhh != candidate_vertices[facet_addition_type] && *vhh != new_polygon_vertices.front() && *vhh != new_polygon_vertices.back()) new_polygon_vertices.push_back(*vhh);
								traverse_stop_flag = true;
							}
						}
						else if(new_polygon_vertices.front() == polygon_vertices.front())
						{
							VertexHandleHandle npv_second_vhh = ++(new_polygon_vertices.begin());
							VertexHandleHandle pv_second_vhh = ++(polygon_vertices.begin());
							VertexHandleHandle pv_third_vhh = pv_second_vhh;
							++pv_third_vhh;
							if(*npv_second_vhh == *pv_third_vhh)
							{
								// Finish traversing at the last vertex of current polygon.
								VertexHandleHandle pv_last_vhh = --(polygon_vertices.end());
								if(vhh == pv_last_vhh)
								{
									if(*vhh != candidate_vertices[facet_addition_type] && *vhh != new_polygon_vertices.front() && *vhh != new_polygon_vertices.back()) new_polygon_vertices.push_back(*vhh);
									traverse_stop_flag = true;
								}
							}
							else if(*npv_second_vhh == *pv_second_vhh)
							{
								// Finish traversing at the first vertex of current polygon.
								if(vhh == polygon_vertices.begin()) traverse_stop_flag = true;
							}
						}
					}
				}
				else
				{
					// The current polygon is a triangle.
					assert(addFacet(*vhh0, *vhh1, *vhh2));
					new_polygon_vertices.clear();
					traverse_stop_flag = true;
				}
				if(traverse_stop_flag) break;
			}
			if(new_polygon_vertices.size() < polygon_vertices.size())
			{
				// The triangulation of a given polygon into a smaller polygon
				// is successful.
				polygon_vertices.clear();
				polygon_vertices = new_polygon_vertices;
				re_triangulation_flag = false;
			}
			else
			{
				// The triangulation of a given polygon into a smaller polygon
				// is not successful. Therefore loosen the orientation eligibility
				// and triangulate current polygon again. 
				//
				// Attention
				//
				// Re-triangulation can only occurs once because if the polygon
				// after the first re-triangulation is the same as the given one,
				// it means that loosening the orientation eligibility does not
				// help triangulating this polygon. In other words, the given
				// polygon is topologically faulty.
				assert(re_triangulation_flag == false);
				re_triangulation_flag = true;
			}
			new_polygon_vertices.clear();
		}
	}
	return eligibility_flag;
}

void SurfaceTopology::insertVerex(VertexHandle vertex_handle, FacetHandle facet_handle)
{
	/// A vertex is inserted into a given facet by partitioning
	/// the facet into three smaller triangular facets and adding
	/// these facets into cell membrane geometry.
	VertexHandle selected_facet_vertices[3];
	for(size_t i = 0; i < 3; ++i) selected_facet_vertices[i] = facet_handle->edges[i]->vertex;
	removeFacet(facet_handle);
	for(size_t i = 0; i < 3; ++i)
	{
		size_t j = i + 1;
		if(j > 2) j = 0;
		assert(addFacet(vertex_handle, selected_facet_vertices[i], selected_facet_vertices[j]));
	}
}

void SurfaceTopology::updateNeighboringFacets(VertexHandle vertex_handle)
{
	FacetHandles fhs = 	getNeighboringFacetHandles(vertex_handle);
	for(FacetHandleHandle fhh = fhs.begin(); fhh != fhs.end(); ++fhh)
	{
		area -= (*fhh)->area;
		(*fhh)->update();
		area += (*fhh)->area;
	}
}

VertexHandle SurfaceTopology::addVertex(const Vertex& vertex)
{
	///
	/// When a new actin filament is created, the next step is
	/// to add a vertex representing this filament into the
	/// vertex pool of cell membrane. During the addition of
	/// a vertex, the vertex and the actin filament it represents
	/// are directionally linked to each other.
	///
	/// 1) The vertex holds a FilamentBranchHandle-type pointer to
	/// the filament and the tree to which this filament belongs.
	/// This operation is performed within Vertex().
	///
	/// 2) The filament also holds a VertexHandle-type pointer
	/// to the vertex. This operation is performed in this
	/// function.
	///
	/// However this function does not connect this vertex with
	/// other existing vertices in the vertex pool. Therefore it
	/// is important to make sure that the end of the new filament
	/// does form a vertex on membrane surface, which is determined
	/// by the mechanical property of cell membrane, before using
	/// this function to add a filament into the vertex pool of
	/// membrane surface.
	///
	vertices.push_back(vertex);
	VertexHandle vertex_handle = --(vertices.end());
	// Add the vertex handle to corresponding filament.
	vertex_handle->filament->setVertex(vertex_handle);
	++n_uncapped_vertex;
	return vertex_handle;
}

bool SurfaceTopology::addFacet(VertexHandle vh1, VertexHandle vh2, VertexHandle vh3)
{
	///
	/// Three given vertices are used to create a new facet on membrane
	/// surface by connecting the new facet to all existing neighboring
	/// facets in the local surface. This function also checks the
	/// eligibility of the new facet to maintain the two-manifold
	/// topological requirement of 3D surface.
	///
	///
	/// Rationale
	///
	/// Incremental facet addition
	///
	/// A new facet is added into existing cell surface, which means that
	/// this new facet must share at least one edge with the existing cell
	/// surface. The reason why the incremental facet addition is required
	/// is that the edges of the new facet must be sorted clockwisely around
	/// the edge list of its incident vertex in order to ensure the correct
	/// construction of two-manifold surface.
	///
	///
	/// Prerequisite
	///
	/// The following requirements must be satisfied in order to add a facet
	///
	///                        vh1
	///                        * ^
	///                       *   *
	///                      *     *
	///                     *_      *
	///                    v_        *
	///                   vh2------->vh3
	///
	/// 1) The edges, vh1-->vh2, vh2-->vh3 and vh3--vh1, must not exist.
	///
	/// 2) If the edges, vh1-->vh3, vh3-->vh2 and vh2--vh1, exist, their
	/// dual edges must 0, which means that the existing edge must be
	/// boundary edge. This is the requirement of two-manifold surface
	/// topology which says that an edge can only be shared by two facets.
	///
	///
	/// Algorithm
	///
	/// The facet is created from vh1, vh2 and vh3:
	///
	/// For each edge of this new facet, there are three possible cases:
	///
	/// 1) The edge is linked to the existing surface through its dual
	/// edge which points to an edge of the existing surface.
	///
	/// 2) The edge is linked to the existing surface through its incident
	/// vertex whose edge list contains the incident edges of the existing
	/// surface.
	///
	/// 3) The edge is linked to the existing surface through its beginning
	/// vertex and its incident vertex is a newly created vertex which does
	/// not have any other incident edges.
	///
	bool eligibility = true;
	VertexHandle vhs[3];
	vhs[0] = vh1;
	vhs[1] = vh2;
	vhs[2] = vh3;
	for(size_t i = 0; i < 3; ++i)
	{
		size_t j = i;
		if(++j > 2) j = 0;
		// Requirement #1
		for(EdgeHandleHandle ehh = ((vhs[j])->edges).begin(); ehh != ((vhs[j])->edges).end(); ++ehh)
		{
			if((*ehh)->prev->vertex == vhs[i])
			{
				eligibility = false;
				break;
			}
		}
		// Requirement #2
		if(eligibility)
		{
			for(EdgeHandleHandle ehh = ((vhs[i])->edges).begin(); ehh != ((vhs[i])->edges).end(); ++ehh)
			{
				if((*ehh)->prev->vertex == vhs[j] && (*ehh)->dual != static_cast<EdgeHandle>(nullptr))
				{
					eligibility = false;
					break;
				}
			}
		}
		else break;
		if(!eligibility) break;
	}
	// Construct a new facet from vh1, vh2 and vh2.
	if(eligibility)
	{
		edges.push_back(Edge(vh1));
		edges.push_back(Edge(vh2));
		edges.push_back(Edge(vh3));
		EdgeHandle eh = edges.end();
		EdgeHandle eh3 = (--eh);
		EdgeHandle eh2 = (--eh);
		EdgeHandle eh1 = (--eh);
		eh1->prev = eh3; eh1->next = eh2;
		eh2->prev = eh1; eh2->next = eh3;
		eh3->prev = eh2; eh3->next = eh1;
		Facet f(eh1, eh2, eh3);
		facets.push_back(f);
		FacetHandle fh = (--facets.end());
		eh1->facet = fh;
		eh2->facet = fh;
		eh3->facet = fh;
		// End of constructing the new facet.
		//
		// Integrate the new facet into the existing cell surface.
		EdgeHandle ehs[3];
		ehs[0] = eh1;
		ehs[1] = eh2;
		ehs[2] = eh3;
		// Go through all the vertices of the newly created facet.
		for(size_t i = 0; i < 3; ++i)
		{
			// For each edge, there are three tasks to accomplish:
			//
			// 1) If the current edge shares an edge with the local surface
			// of the incident vertex, then:
			//
			//     a) link the current edge to the local surface by dual edge
			//     information.
			//
			//     b) insert the current edge to the incident edge list of
			//     the incident vertex clockwisely.
			//
			// 2) If the current edge shares a vertex with the local surface
			// of the incident vertex, then:
			//
			//     Insert the current edge to the incident edge list of the
			//     incident vertex clockwisely.
			//
			// 3) If the incident vertex of the current edge does not have
			// local surface except for the newly created facet, then:
			//
			//     Insert the current edge to the incident edge list of the
			//     incident vertex.
			//
			// Task 1 has higher priority than task 2 and then task 3.
			EdgeHandles& vertex_edges = ehs[i]->vertex->edges;
			// The local surface around the vertex to which the current
			// edge points exists.
			if(!vertex_edges.empty())
			{
				bool share_edge_flag = false;
				bool share_vertex_flag = false;
				for(EdgeHandleHandle ehh = vertex_edges.begin(); ehh != vertex_edges.end(); ++ehh)
				{
					// Case 1
					//
					// The current edge, ehs[i], shares an edge with the local
					// surface around the vertex to which the current edge points.
					if((*ehh)->next->vertex == ehs[i]->prev->vertex)
					{
						(*ehh)->next->dual = ehs[i];
						(ehs[i])->dual = (*ehh)->next;
						vertex_edges.insert(++ehh, ehs[i]);
						// Insert the current edge to the incident edge list of
						// the incident vertex clockwisely.
						share_edge_flag = true;
						break;
					}
				}
				if(!share_edge_flag)
				{
					for(EdgeHandleHandle ehh = vertex_edges.begin(); ehh != vertex_edges.end(); ++ehh)
					{
						// Case 2
						//
						// The current edge, ehs[i], shares a vertex with the local
						// surface around the vertex to which the current edge points.
						if((*ehh)->prev->vertex == ehs[i]->next->vertex)
						{
							vertex_edges.insert(ehh, ehs[i]);
							// Insert the current edge to the incident edge list of
							// the incident vertex clockwisely.
							share_vertex_flag = true;
							break;
						}
					}
				}
				assert(share_edge_flag || share_vertex_flag);
				// Ensure the incremental facet addition from existing surface.
			}
			// Case 3
			//
			// The local surface around the incident vertex does not exist.
			else vertex_edges.push_back(ehs[i]);
		}
		// End of integration.
		area += fh->area;
	}
	return eligibility;
}

void SurfaceTopology::removeFacet(FacetHandle fh)
{
	///
	/// Both the given facet and its connections with the rest of
	/// membrane surface are removed:
	///
	/// 1) Remove the dual connection between the edges comprising
	/// the given facet and their neigboring edges.
	///
	/// 2) Remove the edges comprising the given facet from the edge
	/// list of the vertices that these edges point to.
	///
	/// Then remove the edges and the given facet from the edge and
	/// the facet list of SurfaceTopology.
	///
	EdgeHandle* facet_edges = fh->edges;
	for(size_t i = 0; i < 3; ++i)
	{
		// Remove the pointer to current edge from the 'dual'
		// field of its dual edge.
		if(facet_edges[i]->dual != static_cast<EdgeHandle>(nullptr))
		{
			facet_edges[i]->dual->dual = static_cast<EdgeHandle>(nullptr);
		}
		// Remove the pointer to current edge from the 'edges'
		// field of the pointed vertex.
		VertexHandle& edge_vertex = facet_edges[i]->vertex;
		EdgeHandles& vertex_edges = edge_vertex->edges;
		for(EdgeHandleHandle ehh = vertex_edges.begin(); ehh != vertex_edges.end(); ++ehh)
		{
			if(*ehh == facet_edges[i])
			{
				vertex_edges.erase(ehh);
				break;
			}
		}
	}
	// Remove the edge from the edge list of 'SurfaceTopology'
	for(size_t i = 0; i < 3; ++i) edges.erase(facet_edges[i]);
	// Update membrane area
	area -= fh->area;
	// Finally remove this facet.
	facets.erase(fh);
}

VertexHandles SurfaceTopology::removeVertex(VertexHandle vh)
{
	///
	/// Both the given vertex and the facets containing this vertex
	/// are removed from membrane surface. This function also calls
	/// triangulatePolygonSurface() to make up the polygonal hole
	/// left by vertex removal:
	///
	/// 1) Store its neighboring vertices.
	///
	/// 2) Remove all incident facets.
	///
	/// 3) Disconnect the vertex and corresponding actin filament.
	///
	/// 4) Remove this vertex
	///
	/// 5) Triangulate the polygonal hole left after the vertex is
	/// removed.
	///
	/// Attention:
	///
	/// The corresponding actin filament of the vertex to be removed
	/// must exist when this function is called.
	///
	// Make sure that the actin filament that the given vertex to be
	// removed is based on does exist.
	assert(vh->filament != static_cast<FilamentBranchHandle>(nullptr));
	// Store all neighboring vertices of the given vertex.
	VertexHandles polygon_vertices = getNeighboringVertexHandles(vh);
	// Remove all incident facets of the given vertex.
	EdgeHandles& vertex_edges = vh->edges;
	EdgeHandleHandle ehh = vertex_edges.begin();
	while(ehh != vertex_edges.end())
	{
		// Since all the edges comprising the given facet are also
		// removed from the edge list of the vertex to which these
		// edges are incident, when the facet is removed, it is
		// important to move the handle of current edge forward
		// before removing the incident facet of current edge.
		EdgeHandleHandle ehh_removed = (ehh++);
		removeFacet((*ehh_removed )->facet);
	}
	// Remove the vertex link from the corresponding actin filament.
	vh->filament->setVertex(static_cast<VertexHandle>(nullptr));
	// If the vertex to be removed is not capped, subtract the
	// number of uncapped actin filaments by 1.
	if(!vh->filament->isCapped()) --n_uncapped_vertex;
	// Remove the vertex from the the vertex pool of cell surface.
	vertices.erase(vh);
	// The polygonal hole left by the removal of current vertex
	// must be triangulated.
	assert(triangulatePolygonSurface(polygon_vertices));
	return polygon_vertices;
}

void SurfaceTopology::updateCompositeProperties(VertexHandle vh, bool pos_change, bool cap_change)
{
	///
	/// When filament reactions occur, such as filament growing
	/// and filament capping, some composite properties of cell
	/// membrane, e.g. the area of membrane surface and the
	/// number of dynamic vertex, need to be updated.
	///
	if(pos_change) updateNeighboringFacets(vh);
	if(cap_change)
	{
		if(vh->filament->isCapped()) --n_uncapped_vertex;
		else ++n_uncapped_vertex;
	}
}

VertexHandles SurfaceTopology::getNeighboringVertexHandles(VertexHandle vertex_handle) const
{
	///
	/// The closest neighboring vertices of the given vertex are
	/// the direct neighboring vertices of the given vertex. The
	/// result vertices are oriented counter-clockwisely.
	///
	/// This algorithm also works if holes exist around the local
	/// surface of the given vertex.
	///
	VertexHandles neighboring_vertices;
	EdgeHandles vertex_edges = vertex_handle->edges;
	for(EdgeHandles::reverse_iterator ehh = vertex_edges.rbegin(); ehh != vertex_edges.rend(); ++ehh)
	{
		EdgeHandle curr_edge = (*ehh);
		EdgeHandle next_edge = curr_edge->next;
		EdgeHandle next_dual = next_edge->dual;
		if(next_dual == static_cast<EdgeHandle>(nullptr)) neighboring_vertices.push_back(next_edge->vertex);
		neighboring_vertices.push_back(curr_edge->prev->vertex);
	}
	return neighboring_vertices;
}

FacetHandles SurfaceTopology::getNeighboringFacetHandles(VertexHandle vertex_handle) const
{
	///
	/// The local facets of the given vertex is the facets
	/// surrounding the vertex in its local surface. The result
	/// facets are oriented clockwisely.
	///
	/// This algorithm also works if holes exist around the local
	/// surface of the given vertex.
	///
	FacetHandles vertex_facets;
	EdgeHandles vertex_edges = vertex_handle->edges;
	for(EdgeHandleHandle ehh = vertex_edges.begin(); ehh != vertex_edges.end(); ++ehh)
	{
		vertex_facets.push_back((*ehh)->facet);
	}
	return vertex_facets;
}

FacetHandles SurfaceTopology::getNeighboringFacetHandles(FacetHandle facet_handle) const
{
	///
	/// The local facets of the given facet is the facets surrounding
	/// the facet in its local surface, including both the facets
	/// sharing common edges and the facets sharing common vertices
	/// with the given facet. The result facets are oriented clockwisely.
	///
	/// This algorithm also works if holes exist around the local surface
	/// of the given facet.
	///
	FacetHandles neighboring_facets;
	EdgeHandle* facet_edges = facet_handle->edges;
	for(size_t vertex_index = 0; vertex_index < 3; ++vertex_index)
	{
		EdgeHandles vertex_edges = facet_edges[2 - vertex_index]->vertex->edges;
		EdgeHandleHandle ehh = vertex_edges.begin();
		/// Find the edge to start with.
		while((*ehh)->facet != facet_handle) if((++ehh) == vertex_edges.end()) ++ehh;
		if((++ehh) == vertex_edges.end()) ++ehh;
		/// Start adding the local facets of the given facet.
		FacetHandle current_facet_handle = (*ehh)->facet;
		while(current_facet_handle != facet_handle)
		{
			if(current_facet_handle != neighboring_facets.front() && current_facet_handle != neighboring_facets.back()) neighboring_facets.push_back(current_facet_handle);
			if((++ehh) == vertex_edges.end()) ++ehh;
			current_facet_handle = (*ehh)->facet;
		}
	}
	return neighboring_facets;
}

FacetHandles SurfaceTopology::getLocalFacetHandles(VertexHandle vertex_handle) const
{
	FacetHandles local_facets;
	EdgeHandles& vertex_edges = vertex_handle->edges;
	// Add the neighboring facets of the given vertex.
	for(EdgeHandleHandle ehh = vertex_edges.begin(); ehh != vertex_edges.end(); ++ehh)
	{
		local_facets.push_back((*ehh)->facet);
	}
	// Add the local facets that share edges with the neighboring
	// facets of the given vertex.
	for(EdgeHandleHandle ehh = vertex_edges.begin(); ehh != vertex_edges.end(); ++ehh)
	{
		EdgeHandle prev_dual_edge = (*ehh)->prev->dual;
		FacetHandle local_facet = prev_dual_edge->facet;
		if(local_facet != local_facets.back()) local_facets.push_back(local_facet);
		else local_facets.push_back(prev_dual_edge->prev->dual->facet);
	}
	// Add the neighboring facets of the neighboring vertices of
	// the given vertex and these facets don't overlap with the
	// facets that have already been added above.
	for(EdgeHandleHandle ehh = vertex_edges.begin(); ehh != vertex_edges.end(); ++ehh)
	{
		VertexHandle curr_vertex = (*ehh)->prev->vertex;
		VertexHandle prev_vertex = (*ehh)->dual->next->vertex;
		VertexHandle next_vertex = (*ehh)->next->vertex;
		EdgeHandles& curr_vertex_edges = curr_vertex->edges;
		for(EdgeHandleHandle curr_vertex_ehh = curr_vertex_edges.begin(); curr_vertex_ehh != curr_vertex_edges.end(); ++curr_vertex_ehh)
		{
			VertexHandle curr_vertex_prev = (*curr_vertex_ehh)->prev->vertex;
			VertexHandle curr_vertex_next = (*curr_vertex_ehh)->next->vertex;
			if(curr_vertex_prev != prev_vertex && curr_vertex_prev != next_vertex && curr_vertex_prev != vertex_handle && \
			   curr_vertex_next != prev_vertex && curr_vertex_next != next_vertex && curr_vertex_next != vertex_handle) \
			{
				local_facets.push_back((*curr_vertex_ehh)->facet);
			}
		}
	}
	return local_facets;
}

Orientation SurfaceTopology::computeBranchingOrientation(const CartesianCoordinate& bs, double alpha, const CartesianCoordinate& tip, Triangle& t)
{
	Orientation orient;
	Line line_b_tip(bs, tip);
	Vector vector_b_tip = line_b_tip.getVector();
	double dist_b_tip = vector_b_tip.getMag();
	Vector n = t.getNormal();
	if(!isEqual(dotProd(vector_b_tip, n), 0))
	{
		//
		// If the moether filament is not parallel to the given facet:
		//
		// The child filament is located in the plane determined by
		// the mother filament and the perpendicular line from the
		// branching site on the mother filament to the plane of
		// selected neighboring facet.
		//
		// inside cell
		//                        B
		//                        *
		//                      * |*
		//                    *   | *
		//     mother       *     |  *   child
		//     filament   *       |   *  filament
		//              *TIP      |    *
		//            *         \ | /   *
		//          *            \|/     *
		//        *************************   --> plane of the
		//       V0               P        D      selected facet
		//
		// outside cell
		//
		// v0 is the intersection point between the orientation defined
		// by the Line (B, TIP) and the plane defined by the triangle.
		//
		CartesianCoordinate v0 = intersection(line_b_tip, t);
		double dist_b_v0 = distance(bs, v0);
		CartesianCoordinate p = projection(bs, t);
		double dist_b_p = distance(p, bs);
		double angle_v0_b_p = 0.5 * M_PI - std::asin(dist_b_p / dist_b_v0);
		// angle_p_b_d may be less than zero.
		double angle_p_b_d = alpha - angle_v0_b_p;
		// dist_p_d may be less than zero.
		double dist_p_d = dist_b_p * std::tan(angle_p_b_d);
		double dist_v0_p = distance(p, tip);
		double dist_v0_d = dist_v0_p + dist_p_d;
		Line line_v0_p(tip, p);
		CartesianCoordinate d = line_v0_p.getLocation(dist_v0_d, false);
		orient = Vector(bs, d).getOrient();
	}
	else
	{
		// If the moether filament is parallel to the given facet:
		if(!isEqual(alpha, M_PI / 2))
		{
			// If the branching angle of the child filament is not equal
			// to 90 degree:
			//
			// the child filament is located in the plane determined by
			// the mother filament and the perpendicular line from the
			// branching site on the mother filament to the plane of
			// selected neighboring facet.
			//
			// inside cell
			//
			//      mother      TIP            B
			//      filamen      **************   --> plane of the
			//                   |           *        selected facet
			//                   |          *
			//                   |         *
			//                   |        *
			//                   |       *
			// normal direction  |      *
			// of the selected   |     *  child filament
			// facet             |    *
			//                 \ | / *
			//                  \|/ *
			//                   | *
			//                   *
			//                   D
			// outside cell
			//
			double dist_tip_d = dist_b_tip * std::tan(alpha);
			Line line_tip_d(tip, n);
			CartesianCoordinate d = line_tip_d.getLocation(dist_tip_d, false);
			orient = Vector(bs, d).getOrient();
		}
		else
		{
			// If the branching angle of the child filament is equal to
			// 90 degree
			orient = n.getOrient();
		}
	}
	assert(finite(orient.theta) && finite(orient.phi));
	return orient;
}

bool SurfaceTopology::isBranchingAllowed(FilamentBranch& branch)
{
	// Two rules need to be satisfied to create a new filament.
	// Rule 1. The branching site on mother filament must be located
	// within the cortical region of neighboring facet.
	// Rule 2. The child filament must intersect with neighboring
	// facet.
	bool branching_flag = false;
	FacetHandle selected_branching_facet = static_cast<FacetHandle>(nullptr);
	Orientation selected_branching_orient;
	if(branch.getBranchingSiteActinConstHandle() != static_cast<ActinConstHandle>(nullptr))
	{
		ParameterTable::Table& param_table = ParameterTable::instance();
		double cortical_region_thickness = strtod(param_table[std::string("cortical_region_thickness")]);
		double branching_angle = strtod(param_table[std::string("branching_angle")]) * M_PI / 180;
		CartesianCoordinate branching_site_location = branch.getBranchingSiteActinLocation();
		CartesianCoordinate filament_tip_location = branch.getTailEndLocation();
		FacetHandles cortical_facets;
		FacetHandles searched_facets = getNeighboringFacetHandles(branch.getVertex());
		// If the neighboring facets are not enough to find an acceptable orientation
		// for child filament branch, then try the local facets which has larger pool
		// of facets.
		FacetHandleHandle fhh;
		for(fhh = searched_facets.begin(); fhh != searched_facets.end(); ++fhh)
		{
			Triangle t = (*fhh)->getTriangle();
			// Since the normal vector of a surface facet points towards extra-
			// cellular region, the thickness of the cortical region of the
			// facet should have negative sign.
			if(isInsideTriangularBox(branching_site_location, t, -cortical_region_thickness)) cortical_facets.push_back(*fhh);
		}
		// For each of these candidate facet, calculate the orientation of
		// child filament to see child filament intersects with these facets.
		FacetHandles intersect_facets;
		Orientations intersect_orients;
		for(fhh = cortical_facets.begin(); fhh != cortical_facets.end(); ++fhh)
		{
			Triangle facet_triangle = (*fhh)->getTriangle();
			Orientation child_branch_orient = computeBranchingOrientation(branching_site_location, branching_angle, filament_tip_location, facet_triangle);
			Line child_filament_line(branching_site_location, Vector(branch.getInitialLength(), child_branch_orient));
			if(isIntersecting(child_filament_line, facet_triangle))
			{
				intersect_facets.push_back(*fhh);
				intersect_orients.push_back(child_branch_orient);
			}
		}
		double deviation_angle_min = DBL_INF_POSITIVE;
		OrientationHandle oh = intersect_orients.begin();
		for(fhh = intersect_facets.begin(); fhh != intersect_facets.end(); ++fhh)
		{
			double deviation_angle = computeDeviationAngleOfFilamentGrowth(branch, *oh);
			if(deviation_angle < deviation_angle_min)
			{
				selected_branching_facet = *fhh;
				selected_branching_orient = *oh;
				deviation_angle_min = deviation_angle;
			}
			++oh;
		}
		double max_deviation_angle = strtod(param_table[std::string("max_deviation_angle")]) * M_PI / 180;
		if(deviation_angle_min < max_deviation_angle)
		{
			branch.setChildBranchOrient(selected_branching_orient);
			branch.setChildBranchFacet(selected_branching_facet);
			branching_flag = true;
		}
		else branching_flag = false;
	}
	else branching_flag = false;
	branch.setChildBranchingFlag(branching_flag);
	return branching_flag;
}

FacetHandle SurfaceTopology::getIntersectingFacetHandle(FilamentBranch& branch, FacetHandle start_facet_handle)
{
	FacetHandle intersecting_facet_handle = static_cast<FacetHandle>(nullptr);
	Line l(branch.getHeadEndLocation(), branch.getTailEndLocation());
	Triangle t;
	t = start_facet_handle->getTriangle();
	if(isIntersecting(l, t)) intersecting_facet_handle = start_facet_handle;
	else
	{
		FacetHandles fhs = getNeighboringFacetHandles(start_facet_handle);
		for(FacetHandleHandle fhh = fhs.begin(); fhh != fhs.end(); ++fhh)
		{
			Triangle t = (*fhh)->getTriangle();
			if(isIntersecting(l, t))
			{
				intersecting_facet_handle = *fhh;
				break;
			}
		}
	}
	return intersecting_facet_handle;
}

Orientation SurfaceTopology::computeExtraCellularOrientation(const CartesianCoordinate& loc)
{
	double theta = CylindricalCoordinate(loc).theta;
	Orientation eco(theta, M_PI/2);
	return eco;
}


double SurfaceTopology::computeDeviationAngleOfFilamentGrowth(FilamentBranch& branch)
{
	Vector prefered_growth_vector(1, computeExtraCellularOrientation(branch.getHeadEndLocation()));
	Vector filament_vector(1, branch.getOrient());
	double deviation_angle = angle(filament_vector, prefered_growth_vector);
	return deviation_angle;
}


double SurfaceTopology::computeDeviationAngleOfFilamentGrowth(FilamentBranch& mother_branch, const Orientation& child_branch_orient)
{
	Vector child_branch_orient_vector(1, child_branch_orient);
	CartesianCoordinate branching_site_location = mother_branch.getBranchingSiteActinLocation();
	Vector prefered_growth_vector(1, computeExtraCellularOrientation(branching_site_location));
	double deviation_angle = angle(child_branch_orient_vector, prefered_growth_vector);
	return deviation_angle;
}

bool SurfaceTopology::searchFusionFacetForChildBranch(Line& child_branch_line, VertexHandle mother_branch_vertex, FacetHandle& child_fusion_facet)
{
	bool found_flag = false;
	CartesianCoordinate child_branch_end = child_branch_line.getEnd();
	Vector child_branch_vector = child_branch_line.getVector();
	// Step 1. Find the facets with the same direction as the given line
	// and such facet must exist.
	FacetHandles same_dir_facets;
	FacetHandles local_facets = getLocalFacetHandles(mother_branch_vertex);
	for(FacetHandleHandle fhh = local_facets.begin(); fhh != local_facets.end(); ++fhh)
	{
		Triangle triangle = (*fhh)->getTriangle();
		if(dotProd(child_branch_vector, triangle.getNormal()) > 0) same_dir_facets.push_back(*fhh);
	}
	if(same_dir_facets.size() > 0)
	{
		// Step 2. Find the facets which the end point the given line
		// is located inside.
		FacetHandles inside_box_facets, outside_box_facets;
		for(FacetHandleHandle fhh = same_dir_facets.begin(); fhh != same_dir_facets.end(); ++fhh)
		{
			Triangle triangle = (*fhh)->getTriangle();
			if(isProjectionInsideTriangle(child_branch_end, triangle)) inside_box_facets.push_back(*fhh);
			else outside_box_facets.push_back(*fhh);
		}
		// Step 3. Select the inside-box facet to which the end point is
		// closest.
		if(inside_box_facets.size() > 0)
		{
			double dist_min = DBL_INF_POSITIVE;
			for(FacetHandleHandle fhh = inside_box_facets.begin(); fhh != inside_box_facets.end(); ++fhh)
			{
				Triangle triangle = (*fhh)->getTriangle();
				double dist = distance(child_branch_end, triangle);
				if(dist < dist_min)
				{
					child_fusion_facet = *fhh;
					dist_min = dist;
					found_flag = true;
				}
			}
		}
	}
	return found_flag;
}

void SurfaceTopology::exportGeometry(std::ofstream& output, char delim)
{
	output << "OFF" << std::endl;
	output << vertices.size() << delim << facets.size() << " 0" << std::endl;
	VertexHandle vh;
	for(vh = vertices.begin(); vh != vertices.end(); ++vh)
	{
		CartesianCoordinate p = vh->filament->getTailEndLocation();
		output << p.x << delim << p.y << delim << p.z << std::endl;
	}
	FacetHandle fh;
	for(fh = facets.begin(); fh != facets.end(); ++fh)
	{
		output << '3';
		for(size_t i = 0; i < 3; ++i) output << delim << getVertexHandlePosition(fh->edges[i]->vertex);
		// Set facet color to purple.
		output << delim << "0.666" << delim << "0.666" << delim << "0.888";
		output << std::endl;
	}
}

}
