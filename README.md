This implementation is 100% Public Domain.

Feel free to use.
 
C++11 is needed to compile it.

Usage:

	#include "quickhull/quickhull.hpp"

	using namespace quickhull;
	QuickHull<float> qh; // Could be double as well
	std::vector<Vector3<float>> pointCloud;
	// Add points to point cloud
	...
	auto hull = qh.GetConvexHull(pointCloud, true); // Change true to false to get non-CCW triangles
	auto indexBuffer = hull.getIndexBuffer();
	auto vertexBuffer = hull.getVertexBuffer();
	// Do what you want with the convex triangle mesh
