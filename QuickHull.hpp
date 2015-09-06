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
#include "Structs/Vector3.hpp"
#include "Structs/Plane.hpp"
#include "Structs/Face.hpp"

/*
 * Implementation of the 3D QuickHull algorithm by Antti Kuukka
 *
 * No copyrights. What follows is 100% Public Domain.
 *
 *
 *
 * INPUT:  a list of points in 3D space (for example, vertices of a 3D mesh)
 *
 * OUTPUT: a 3D mesh which is the convex hull for the given list of points. The output mesh is a list of vertex positions where
 *         three consecutive positions represent a triangle.
 *
 *
 * The implementation is thread-safe.
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
 *              - Because of the way M is constructed, these faces are connected. Solve their horizon edge.
 *				- Create new faces by connecting P with the points belonging to the horizon edge. Add the new faces to M and remove the visible
 *                faces from M.
 *              - Each point that was assigned to visible faces is now assigned to at most one of the newly created faces.
 *              - Those new faces that have points assigned to them are added to the top of Face Stack.
 *          - M is now the convex hull.
 *
 *          Notes:
 *          - When faces of M are removed during the iteration, they are not actually removed from the std::vector holding face data. Instead,
 *            we "disable" the face and when new faces are added to the mesh, disabled faces get replaced. This reduces the number of
 *            allocations and makes the algorithm run faster.
 *
 *
 *
 * TO DO:
 *  - Implement a proper 2D QuickHull and use that to solve the degenerate 2D case (when all the points lie on the same plane in 3D space).
 *  - Make the public interface more flexible (accept vertex data as const float*, const double*, etc)
 *  - Use a half-edge data structure. Should improve performance!
 * */

namespace quickhull {
	
	template<typename T>
	class QuickHull {
		static const T Epsilon;
		
		// Detect degenerate cases
		static std::vector<Vector3<T>> checkDegenerateCase0D(const std::vector<Vector3<T>>& pointCloud, T epsilon);
		static std::vector<Vector3<T>> checkDegenerateCase1D(const std::vector<Vector3<T>>& pointCloud, T epsilon);
		static std::vector<Vector3<T>> checkDegenerateCase2D(const std::vector<Vector3<T>>& pointCloud, T epsilon);
		
		static std::vector<Face<T>> getInitialTetrahedron(const std::vector<Vector3<T>>& vPositions, const std::array<IndexType,6>& extremePoints, T epsilon);
		
		// Find indices of extreme values (max x, min x, max y, min y, max z, min z) for the given point cloud
		static std::array<IndexType,6> findExtremeValues(const std::vector<Vector3<T>>& vPositions);
		
		// Creates a 3D mesh from a list of Face structs
		static std::vector<Vector3<T>> createMesh(const std::vector<Face<T>>& faces, const std::vector<Vector3<T>>& vPositions, bool CCW);
	public:
		// Computes convex hull for a given point cloud.
		// Params:
		//   pointCloud: a list of 3D points
		//   CCW: whether the output mesh triangles should have CCW orientation
		static std::vector<Vector3<T>> getConvexHull(const std::vector<Vector3<T>>& pointCloud, bool CCW);
	};

}


#endif /* QUICKHULL_HPP_ */
