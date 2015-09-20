/*
 * QuickHull.hpp
 *
 *  Created on: May 4, 2014
 *      Author: anttiku
 */

#ifndef QUICKHULL_HPP_
#define QUICKHULL_HPP_

#include <vector>
#include <array>
#include <limits>
#include "Structs/Vector3.hpp"
#include "Structs/Plane.hpp"
#include "Structs/HalfEdgeMesh.hpp"
#include "ConvexHull.hpp"

/*
 * Implementation of the 3D QuickHull algorithm by Antti Kuukka
 *
 * No copyrights. What follows is 100% Public Domain.
 *
 *
 *
 * INPUT:  a list of points in 3D space (for example, vertices of a 3D mesh)
 *
 * OUTPUT: a ConvexHull object which provides vertex and index buffers of the generated convex hull as a triangle mesh.
 *
 *
 *
 * The implementation is thread-safe if each thread is using its own QuickHull object.
 *
 *
 * SUMMARY OF THE ALGORITHM:
 *         - Create initial simplex (tetrahedron) using extreme points. We have four faces now and they form a convex mesh M.
 *         - For each point, assign them to the first face for which they are on the positive side of (so each point is assigned to at most
 *           one face). Points inside the initial tetrahedron are left behind now and no longer affect the calculations.
 *         - Add all faces that have points assigned to them to Face Stack.
 *         - Iterate until Face Stack is empty:
 *              - Pop topmost face F from the stack
 *              - From the points assigned to F, pick the point P that is farthest away from the plane defined by F.
 *              - Find all faces of M that have P on their positive side. Let us call these the "visible faces".
 *              - Because of the way M is constructed, these faces are connected. Solve their horizon edge loop.
 *				- "Extrude to P": Create new faces by connecting P with the points belonging to the horizon edge. Add the new faces to M and remove the visible
 *                faces from M.
 *              - Each point that was assigned to visible faces is now assigned to at most one of the newly created faces.
 *              - Those new faces that have points assigned to them are added to the top of Face Stack.
 *          - M is now the convex hull.
 *
 * TO DO:
 *  - Implement a proper 2D QuickHull and use that to solve the degenerate 2D case (when all the points lie on the same plane in 3D space).
 *  - Make the public interface more flexible (accept vertex data as const float*, const double*, etc)
 *  - Investigate possibility to use unnormalized triangle normals to gain some speed
 * */

namespace quickhull {

	template<typename T>
	class QuickHull {
		static const T Epsilon;

		T m_epsilon;
		const std::vector<Vector3<T>>* m_vertexData;
		Mesh<T> m_mesh;
		std::array<IndexType,6> m_extremeValues;

		// Temporary variables used during iteration process
		std::vector<IndexType> m_newFaceIndices;
		std::vector<IndexType> m_newHalfEdgeIndices;

		// Detect degenerate cases
		ConvexHull<T> checkDegenerateCase0D();
		ConvexHull<T> checkDegenerateCase1D();
		ConvexHull<T> checkDegenerateCase2D();
		
		// Create a half edge mesh representing the base tetrahedron from which the QuickHull iteration proceeds. m_extremeValues must be properly set up when this is called.
		Mesh<T> getInitialTetrahedron();

		// Given a list of half edges, try to rearrange them so that they form a loop. Return true on success.
		bool reorderHorizonEdges(std::vector<IndexType>& horizonEdges);
		
		// Find indices of extreme values (max x, min x, max y, min y, max z, min z) for the given point cloud
		std::array<IndexType,6> findExtremeValues(const std::vector<Vector3<T>>& vPositions);

		// This will update m_mesh from which we create the ConvexHull object that getConvexHull function returns
		void createConvexHalfEdgeMesh();
	public:
		// Computes convex hull for a given point cloud.
		// Params:
		//   pointCloud: a list of 3D points
		//   CCW: whether the output mesh triangles should have CCW orientation
		ConvexHull<T> getConvexHull(const std::vector<Vector3<T>>& pointCloud, bool CCW);
	};

}


#endif /* QUICKHULL_HPP_ */
