#include "QuickHull.hpp"
#include "Structs/Plane.hpp"
#include "Structs/EdgeLoop.hpp"
#include "Structs/Face.hpp"
#include "MathUtils.hpp"
#include <cmath>
#include <cassert>
#include <functional>
#include <iostream>
#include <algorithm>

namespace quickhull {
	
	template<>
	const float QuickHull<float>::Epsilon = 0.0001f;
	
	template<>
	const double QuickHull<double>::Epsilon = 0.00001;
	
	/*
	 * Implementation of the algorithm
	 */
	
	template<typename T>
	std::vector<Vector3<T>> QuickHull<T>::getConvexHull(const std::vector<Vector3<T>>& pointCloud, bool CCW) {
		if (pointCloud.size()==0) {
			return std::vector<Vector3<T>>();
		}
		
		// Very first: find extreme values and use them to compute the scale of the point cloud.
		auto extremes = findExtremeValues(pointCloud);
		const Vector3<T> maxs(pointCloud[extremes[0]].x,pointCloud[extremes[2]].y,pointCloud[extremes[4]].z);
		const Vector3<T> mins(pointCloud[extremes[1]].x,pointCloud[extremes[3]].y,pointCloud[extremes[5]].z);
		const T scale = std::max(mins.getLength(),maxs.getLength());
		
		// Epsilon we use depends on the scale
		const T eps = Epsilon*scale;
		
		// Check for degenerate cases before proceeding
		std::function<std::vector<Vector3<T>>(const std::vector<Vector3<T>>&,T)> degenerateCaseCheckers[] = {checkDegenerateCase0D,checkDegenerateCase1D,checkDegenerateCase2D};
		for (auto& f : degenerateCaseCheckers) {
			auto degenerateMesh = f(pointCloud,eps);
			if (degenerateMesh.size()) {
				return degenerateMesh;
			}
		}
		
		// Temporary variables used during iteration
		std::vector<IndexType> visibleFaces;
		std::vector<IndexType> disabledFaces;
		std::vector<IndexType> newFaceIndices;
		EdgeLoop<T> el;
		
		// Compute base tetrahedron
		std::vector<Face<T>> mesh = getInitialTetrahedron(pointCloud, extremes, eps);
		assert(mesh.size()==4);
		
		// Init face stack with those faces that have points assigned to them
		std::vector<IndexType> faceStack;
		for (size_t i=0;i < mesh.size();i++) {
			if (mesh[i].pointsOnPositiveSide.size()>0) {
				faceStack.push_back(i);
			}
		}
		
		// Process faces until face stack is empty.
		while (faceStack.size() > 0) {
			// Pop topmost face from the stack
			const IndexType topFaceIndex = *(faceStack.end()-1);
			faceStack.erase(faceStack.end()-1);
			Face<T>& tf = mesh[topFaceIndex];
			if (tf.pointsOnPositiveSide.size()==0 || tf.isDisabled()) {
				continue;
			}
			
			// Find point farthest away from the triangle plane
			T maxD=0;
			size_t maxI=0;
			for (auto pointIndex : tf.pointsOnPositiveSide) {
				const T d = mathutils::getSignedDistanceToPlane(pointCloud[pointIndex],tf.P);
				if (d > maxD) {
					maxI=pointIndex;
					maxD=d;
				}
			}
			
			const Vector3<T>& activePoint = pointCloud[maxI];
			const size_t activePointIndex = maxI;
			
			// Get all faces that are visible from this point
			visibleFaces.clear();
			const size_t s = mesh.size();
			for (size_t i=0;i<s;i++) {
				if (!mesh[i].isDisabled()) {
					Plane<T>& P = mesh[i].P;
					T d = P.m_N.dotProduct(activePoint)+P.m_D;
					if (d>=0) {
						visibleFaces.push_back((IndexType)i);
					}
				}
			}
			
			// Create edge loop from the visible faces
			el.createFrom(visibleFaces,mesh);
			if (el.getEdgeCount()) {
				newFaceIndices.clear();
				
				// Create new triangles from the edges
				const size_t newFaceCount = el.getEdgeCount();
				for (IndexType i=0;i<newFaceCount;i++) {
					const auto& e = el.getEdge(i);
					const Vector3<T> triangleNormal = (pointCloud[e.first]-activePoint).crossProduct(pointCloud[e.second]-activePoint).getNormalized();
					Plane<T> trianglePlane(triangleNormal,activePoint);
					if (disabledFaces.size()) {
						auto it = (disabledFaces.end()-1);
						IndexType newIndex = *it;
						disabledFaces.erase(it);
						auto& f = mesh[newIndex];
						f.pointsOnPositiveSide.clear();
						f.a = activePointIndex;
						f.b = e.first;
						f.c = e.second;
						f.P = trianglePlane;
						newFaceIndices.push_back(newIndex);
					}
					else {
						newFaceIndices.push_back(mesh.size());
						mesh.emplace_back(activePointIndex,e.first,e.second,trianglePlane);
					}
				}
				
				// Disable the faces and assign points that were visible from them to the new faces
				for (auto faceIndex : visibleFaces) {
					auto& disabledFace = mesh[faceIndex];
					disabledFace.disable();
					disabledFaces.push_back(faceIndex);
					for (const auto& point : disabledFace.pointsOnPositiveSide) {
						if (point == activePointIndex) {
							continue;
						}
						for (size_t j=0;j<newFaceIndices.size();j++) {
							if (mathutils::getSignedDistanceToPlane(pointCloud[ point ],mesh[newFaceIndices[j]].P)>eps) {
								mesh[newFaceIndices[j]].pointsOnPositiveSide.push_back( point );
								break;
							}
						}
					}
				}
				
				// Increase face stack size if needed
				for (auto newFaceIndex : newFaceIndices) {
					if (mesh[newFaceIndex].pointsOnPositiveSide.size()) {
						faceStack.push_back(newFaceIndex);
					}
				}
			}
			else {
				// Due to limitations of floating point math, we failed to create edge loop from the visible faces.
#ifdef DEBUG
				std::cout << "ERROR: Failed to solve horizon edge." << std::endl;
#endif
				auto it = std::find(tf.pointsOnPositiveSide.begin(),tf.pointsOnPositiveSide.end(),activePointIndex);
				assert(it != tf.pointsOnPositiveSide.end());
				tf.pointsOnPositiveSide.erase(it);
				// Put active face back to the stack. Perhaps we can solve the horizon edge for some other points.
				faceStack.push_back(topFaceIndex);
			}
			
		}
		return createMesh(mesh, pointCloud, CCW);
	}
	
	/*
	 * Implementations of degenerate case handler functions
	 */
	 
	template<typename T>
	std::vector<Vector3<T>> QuickHull<T>::checkDegenerateCase0D(const std::vector<Vector3<T>>& pointCloud, T epsilon) {
		// 0D degenerate case: all points are at the same location
		const Vector3<T>& v0 = *pointCloud.begin();
		for (const auto& v : pointCloud) {
			T d = (v-v0).getLengthSquared();
			if (d>epsilon*epsilon) {
				return std::vector<Vector3<T>>();
			}
		}
#ifdef DEBUG
		std::cout << "Detected 0D degenerate case: all points are at " << pointCloud[0] << "\n";
#endif
		return std::vector<Vector3<T>>(3,pointCloud[0]);
	}

	template <typename T>
	std::vector<Vector3<T>> QuickHull<T>::checkDegenerateCase1D(const std::vector<Vector3<T>>& pointCloud, T epsilon) {
		// 1D degenerate case: the points form a line in 3D space. To keep things simple, we translate the points so that the first point resides at the origin. This way, if the points do actually form a line, we have a line passing through the origin.
		const Vector3<T>* firstPoint = nullptr;
		T firstPointLengthSquared=0;
		const Vector3<T>* maxPoint = nullptr;
		const Vector3<T>* minPoint = nullptr;
		const T epsilonSquared = epsilon*epsilon;
		T maxDot = -1.0f;
		T minDot = 1.0f;

		// First find a point which does not reside at the origin. Such a point exists, for otherwise we would have the 0D case.
		for (const auto& v : pointCloud) {
			const auto v2 = v-pointCloud[0];
			if (firstPoint == nullptr) {
				if (v2.getLengthSquared() > epsilonSquared) {
					firstPoint = &v;
					firstPointLengthSquared = v2.getLengthSquared();
					break;
				}
			}
		}
		assert(firstPoint != nullptr);

		// Then check if all other translated points point to the same direction as (firstPoint-pointCloud[0]) - or lie at the origin. If not, proceed to checking the 2D degenerate case. While looping, keep track of min and max dot product because if the point cloud turns out to be a line, the min and max points are its end points.
		for (const auto& v : pointCloud) {
			const auto v2 = v-pointCloud[0];
			const T dot = v2.dotProduct(*firstPoint-pointCloud[0]);
			if (v2.getLengthSquared()>epsilon*epsilon) {
				const T V = dot*dot/(v2.getLengthSquared()*firstPointLengthSquared);
				const T d = std::abs(V-1);
				if (d > Epsilon) {
					// The points do not form a line!
					return std::vector<Vector3<T>>();
				}
			}
			if (maxPoint == nullptr || dot > maxDot) {
				maxDot = dot;
				maxPoint = &v;
			}
			if (minPoint == nullptr || dot < minDot) {
				minDot = dot;
				minPoint = &v;
			}
		}

		// We have a degenerate 1D case. Find a third point so that we can construct an infinitely thin triangle.
#ifdef DEBUG
		std::cout << "Detected 1D degenerate case: the point cloud forms a line between " << *minPoint << " and " << *maxPoint << std::endl;
#endif
		const Vector3<T>* thirdPoint = nullptr;
		for (const auto& v : pointCloud) {
			if (&v != minPoint && &v != maxPoint) {
				thirdPoint = &v;
				break;
			}
		}
		if (thirdPoint==nullptr) {
			thirdPoint = minPoint;
		}
		std::vector<Vector3<T>> triangle{*minPoint,*thirdPoint,*maxPoint};
		return triangle;
	}

	template<typename T>
	std::vector<Vector3<T>> QuickHull<T>::checkDegenerateCase2D(const std::vector<Vector3<T>>& pointCloud, T epsilon) {
		// 2D degenerate case: all points lie on the same plane. Just like in the 1D case, we translate the points so that the first point is located at the origin.
		const T epsilonSquared = epsilon*epsilon;

		// Find two points not lying at the origin and not pointing to the same direction (there must be at least two, for otherwise the 0D or 1D cases would have generated the convex hull)
		const Vector3<T>* firstPoint = nullptr;
		const Vector3<T>* secondPoint = nullptr;
		T firstPointLengthSquared = 0.0f;
		for (const auto& v : pointCloud) {
			const auto& vt = v-pointCloud[0];
			if (vt.getLengthSquared()>epsilonSquared) {
				if (firstPoint == nullptr) {
					firstPoint = &v;
					firstPointLengthSquared = vt.getLengthSquared();
					continue;
				}

				const T dot = vt.dotProduct(*firstPoint - pointCloud[0]);
				const T V = dot*dot/(vt.getLengthSquared()*firstPointLengthSquared);
				const T d = std::abs(V-1);
				if (d > Epsilon) {
					secondPoint = &v;
					break;
				}

			}
		}
		assert(firstPoint != nullptr && secondPoint != nullptr);

		// Now firstPoint and secondPoint define a plane. Its normal is their cross product.
		const Vector3<T> cross = (*firstPoint-pointCloud[0]).crossProduct(*secondPoint-pointCloud[0]).getNormalized();
		for (const auto& v : pointCloud) {
			const auto& vt = v-pointCloud[0];
			const T d = std::abs(vt.dotProduct(cross));
			if (d > epsilon) {
				// We have a proper 3D point cloud and the QuickHull algorithm can be applied.
				return std::vector<Vector3<T>>();
			}
		}

		// We have encountered the degenerate 2D case. We solve the problem by adding one extra point to the cloud, applying QuickHull, and then removing extra faces. TODO: implement proper 2D QuickHull, because this solution is not good performance-wise.
#ifdef DEBUG
		std::cout << "Degenerate 2D case detected." << std::endl;
#endif
		auto newPoints = pointCloud;
		const T M = epsilon/Epsilon;
		auto extraPoint = M*cross + pointCloud[0];
		newPoints.push_back(extraPoint);
		auto newHull = QuickHull<T>::getConvexHull(newPoints, true);
		assert(newHull.size()>0);
		size_t faceCount = newHull.size()/3;
		const auto* vData = &newHull[0];
		std::vector<Vector3<T>> finalHull;
		for (size_t i=0;i<faceCount;i++) {
			const auto& v1 = vData[i*3+0];
			const auto& v2 = vData[i*3+1];
			const auto& v3 = vData[i*3+2];
			if (!((v1-extraPoint).getLengthSquared()<M/4 || (v2-extraPoint).getLengthSquared()<M/4 || (v3-extraPoint).getLengthSquared()<M/4)) {
				finalHull.push_back(v1);
				finalHull.push_back(v2);
				finalHull.push_back(v3);
			}
		}

		return finalHull;
	}
	
	/*
	 * Private helper functions
	 */

	template <typename T>
	std::array<IndexType,6> QuickHull<T>::findExtremeValues(const std::vector<Vector3<T>>& vPositions) {
		std::array<IndexType,6> outIndices{0,0,0,0,0,0};
		T extremeVals[6] = {vPositions[0].x,vPositions[0].x,vPositions[0].y,vPositions[0].y,vPositions[0].z,vPositions[0].z};
		const size_t vCount = vPositions.size();
		for (size_t i=1;i<vCount;i++) {
			const Vector3<T>& pos = vPositions[i];
			if (pos.x>extremeVals[0]) {
				extremeVals[0]=pos.x;
				outIndices[0]=(IndexType)i;
			}
			if (pos.x<extremeVals[1]) {
				extremeVals[1]=pos.x;
				outIndices[1]=(IndexType)i;
			}
			if (pos.y>extremeVals[2]) {
				extremeVals[2]=pos.y;
				outIndices[2]=(IndexType)i;
			}
			if (pos.y<extremeVals[3]) {
				extremeVals[3]=pos.y;
				outIndices[3]=(IndexType)i;
			}
			if (pos.z>extremeVals[4]) {
				extremeVals[4]=pos.z;
				outIndices[4]=(IndexType)i;
			}
			if (pos.z<extremeVals[5]) {
				extremeVals[5]=pos.z;
				outIndices[5]=(IndexType)i;
			}
		}
		return outIndices;
	}

	template <typename T>
	std::vector<Face<T>> QuickHull<T>::getInitialTetrahedron(const std::vector<Vector3<T>>& vPositions, const std::array<IndexType,6>& extremePoints, T epsilon) {
		// Find two most distant extreme points.
		T maxD = 0.0f;
		std::pair<IndexType,IndexType> selectedPoints;
		for (IndexType i=0;i<6;i++) {
			for (IndexType j=i+1;j<6;j++) {
				const T d = vPositions[ extremePoints[i] ].getSquaredDistanceTo( vPositions[ extremePoints[j] ] );
				if (d > maxD) {
					maxD=d;
					selectedPoints=std::pair<IndexType,IndexType>(extremePoints[i],extremePoints[j]);
				}
			}
		}
		assert(maxD > 0.0);

		// Find the most distant point to the line between the two chosen extreme points.
		Ray<T> r(vPositions[selectedPoints.first], (vPositions[selectedPoints.second] - vPositions[selectedPoints.first]).getNormalized());
		maxD=0.0f;
		IndexType maxI=-1;
		const size_t vCount = vPositions.size();
		for (size_t i=0;i<vCount;i++) {
			T distToRay = mathutils::getSquaredDistanceBetweenPointAndRay(vPositions[i],r);
			if (distToRay > maxD) {
				maxD=distToRay;
				maxI=(IndexType)i;
			}
		}
		assert(maxI!=-1);

		// These three points form the base triangle for our tetrahedron.
		std::vector<Face<T>> faces;
		const Face<T> baseTriangle( selectedPoints.first, selectedPoints.second, maxI);
		const Vector3<T> baseTriangleVertices[]={ vPositions[baseTriangle.a], vPositions[baseTriangle.b],  vPositions[baseTriangle.c] };
		faces.push_back(baseTriangle);
		
		// Next step is to find the 4th vertex of the tetrahedron. We naturally choose the point farthest away from the triangle plane.
		maxD=0.0f;
		maxI = 0;
		const Vector3<T> N = (baseTriangleVertices[1]-baseTriangleVertices[0]).crossProduct(baseTriangleVertices[2]-baseTriangleVertices[0]).getNormalized();
		Plane<T> trianglePlane(N,vPositions[faces[0].a]);
		for (IndexType i=0;i<vCount;i++) {
			const T d = std::abs(mathutils::getSignedDistanceToPlane(vPositions[i],trianglePlane));
			if (d>maxD) {
				maxD=d;
				maxI=i;
			}
		}

		// Now that we have the 4th point, we can create the tetrahedron
		Plane<T> triPlane(N,baseTriangleVertices[0]);
		if (triPlane.isPointOnPositiveSide(vPositions[maxI])) {
			// Enforce CCW orientation
			IndexType t = faces[0].a;
			faces[0].a = faces[0].b;
			faces[0].b = t;
		}
		Face<T> f = faces[0];
		faces.emplace_back(f.a,f.c,maxI);
		faces.emplace_back(f.c,f.b,maxI);
		faces.emplace_back(f.b,f.a,maxI);

		// Compute the plane defined by each face
		for (auto& face : faces) {
			Vector3<T> triangleNormal = (vPositions[face.b]-vPositions[face.a]).crossProduct(vPositions[face.c]-vPositions[face.a]).getNormalized();
			Plane<T> trianglePlane(triangleNormal,vPositions[face.a]);
			face.P = trianglePlane;
		}

		// Finally we assign a face for each vertex outside the tetrahedron
		for (IndexType i=0;i<vCount;i++) {
			for (auto& face : faces) {
				if (mathutils::getSignedDistanceToPlane(vPositions[i],face.P) > epsilon) {
					face.pointsOnPositiveSide.push_back( i );
					break;
				}
			}
		}
		return faces;
	}
	
	template <typename T>
	std::vector<Vector3<T>> QuickHull<T>::createMesh(const std::vector<Face<T>>& faces, const std::vector<Vector3<T>>& vPositions, bool CCW) {
		std::vector<Vector3<T>> mesh;
		for (const auto& face : faces) {
			if (!face.isDisabled()) {
				mesh.push_back(vPositions[face.a]);
				if (CCW) {
					mesh.push_back(vPositions[face.b]);
					mesh.push_back(vPositions[face.c]);
				}
				else {
					mesh.push_back(vPositions[face.c]);
					mesh.push_back(vPositions[face.b]);
				}
			}
		}
		return mesh;
	}
	
	/*
	 * Explicit template specifications for float and double
	 */

	template class QuickHull<float>;
	template class QuickHull<double>;

}

