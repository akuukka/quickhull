This implementation is 100% Public Domain.

Feel free to use.
 
C++11 is needed to compile it.

Basic usage:

	#include "quickhull/quickhull.hpp"

	using namespace quickhull;
	QuickHull<float> qh; // Could be double as well
	std::vector<Vector3<float>> pointCloud;
	// Add points to point cloud
	...
	auto hull = qh.getConvexHull(pointCloud, true); // Change true to false to get non-CCW triangles
	auto indexBuffer = hull.getIndexBuffer();
	auto vertexBuffer = hull.getVertexBuffer();
	// Do what you want with the convex triangle mesh

Vertex data can be passed as a pointer to float/double as long as the data is in X_0,Y_0,Z_0,X_1,Y_1,Z_1,...,X_N,Y_N_Z_N format:

	auto hull = qh.getConvexHull(&pointCloud[0].x, pointCloud.size(), true);
