#include "QuickHull.hpp"
#include "Structs/Plane.hpp"
#include "MathUtils.hpp"
#include <cmath>
#include <cassert>
#include <functional>
#include <iostream>
#include <algorithm>
#include <limits>
#include "Structs/HalfEdgeMesh.hpp"

namespace quickhull {
	
	template<>
	const float QuickHull<float>::Epsilon = 0.00001f;
	
	template<>
	const double QuickHull<double>::Epsilon = 0.000001;

	/*
	 * Implementation of the algorithm
	 */

	template<typename T>
	ConvexHull<T> QuickHull<T>::getConvexHull(const std::vector<Vector3<T>>& pointCloud, bool CCW) {
		if (pointCloud.size()==0) {
			return ConvexHull<T>();
		}
		m_vertexData = &pointCloud;

		// Very first: find extreme values and use them to compute the scale of the point cloud.
		m_extremeValues = findExtremeValues(pointCloud);
		const Vector3<T> maxs(pointCloud[m_extremeValues[0]].x,pointCloud[m_extremeValues[2]].y,pointCloud[m_extremeValues[4]].z);
		const Vector3<T> mins(pointCloud[m_extremeValues[1]].x,pointCloud[m_extremeValues[3]].y,pointCloud[m_extremeValues[5]].z);
		const T scale = std::max(mins.getLength(),maxs.getLength());
		
		// Epsilon we use depends on the scale
		m_epsilon = Epsilon*scale;

		// Check for degenerate cases before proceeding to the 3D quickhull iteration phase
		std::function<ConvexHull<T>(void)> degenerateCaseCheckers[] = {
			std::bind(&QuickHull<T>::checkDegenerateCase0D,this),
			std::bind(&QuickHull<T>::checkDegenerateCase1D,this),
			std::bind(&QuickHull<T>::checkDegenerateCase2D,this)
		};
		for (auto& f : degenerateCaseCheckers) {
			auto degenerateHull = f();
			if (degenerateHull.getVertexBuffer().size()) {
				return degenerateHull;
			}
		}
		
		createConvexHalfEdgeMesh();
		return ConvexHull<T>(m_mesh,*m_vertexData, CCW);
	}

	template<typename T>
	void QuickHull<T>::createConvexHalfEdgeMesh() {
		const auto& pointCloud = *m_vertexData;
		const T epsilonSquared = m_epsilon*m_epsilon;
		
		// Temporary variables used during iteration
		std::vector<IndexType> visibleFaces;
		std::vector<IndexType> horizonEdges;
		struct FaceData {
			IndexType m_faceIndex;
			IndexType m_enteredFromHalfEdge; // If the face turns out not to be visible, this half edge will be marked as horizon edge
			FaceData(IndexType fi, IndexType he) : m_faceIndex(fi),m_enteredFromHalfEdge(he) {}
		};
		std::vector<FaceData> possiblyVisibleFaces;

		// Compute base tetrahedron
		m_mesh = getInitialTetrahedron();
		assert(m_mesh.m_faces.size()==4);

		// Init face stack with those faces that have points assigned to them
		std::vector<IndexType> faceStack;
		for (size_t i=0;i < 4;i++) {
			auto& f = m_mesh.m_faces[i];
			if (f.m_pointsOnPositiveSide && f.m_pointsOnPositiveSide->size()>0) {
				faceStack.push_back(i);
				f.m_inFaceStack = 1;
			}
		}

		// Process faces until face stack is empty.
		size_t iter = 0;
		while (faceStack.size() > 0) {
			iter++;
			if (iter == std::numeric_limits<size_t>::max()) {
				// Visible face traversal marks visited faces with iteration counter (to mark that the face has been visited on this iteration) and the max value represents unvisited faces. At this point we have to reset iteration counter. This shouldn't be an
				// issue on 64 bit machines.
				iter = 0;
			}

			// Pop topmost face from the stack
			const IndexType topFaceIndex = *(faceStack.end()-1);
			faceStack.erase(faceStack.end()-1);
			auto& tf = m_mesh.m_faces[topFaceIndex];
			tf.m_inFaceStack = 0;

			assert(!tf.m_pointsOnPositiveSide || tf.m_pointsOnPositiveSide->size()>0);
			if (!tf.m_pointsOnPositiveSide || tf.isDisabled()) {
				continue;
			}
			
			// Pick the most distant point to this triangle plane as the point to which we extrude
			const Vector3<T>& activePoint = pointCloud[tf.m_mostDistantPoint];
			const size_t activePointIndex = tf.m_mostDistantPoint;

			// Find out the faces that have our active point on their positive side (these are the "visible faces"). The face on top of the stack of course is one of them. At the same time, we create a list of horizon edges.
			horizonEdges.clear();
			possiblyVisibleFaces.clear();
			visibleFaces.clear();
			tf.m_visibilityCheckedOnIteration = std::numeric_limits<size_t>::max();
			possiblyVisibleFaces.emplace_back(topFaceIndex,std::numeric_limits<size_t>::max());
			while (possiblyVisibleFaces.size()) {
				auto it = (possiblyVisibleFaces.end()-1);
				const auto& faceData = *it;
				possiblyVisibleFaces.erase(it);
				auto& pvf = m_mesh.m_faces[faceData.m_faceIndex];
				assert(!pvf.isDisabled());

				if (pvf.m_visibilityCheckedOnIteration == iter && pvf.m_isVisibleFaceOnCurrentIteration) {
					continue;
				}
				else {
					Plane<T>& P = pvf.m_P;
					pvf.m_visibilityCheckedOnIteration = iter;
					T d = P.m_N.dotProduct(activePoint)+P.m_D;
					if (d>=0) {
						pvf.m_isVisibleFaceOnCurrentIteration = 1;
						pvf.m_horizonEdgesOnCurrentIteration = 0;
						visibleFaces.push_back(faceData.m_faceIndex);
						for (auto heIndex : m_mesh.getHalfEdgeIndicesOfFace(pvf)) {
							if (m_mesh.m_halfEdges[heIndex].m_opp != faceData.m_enteredFromHalfEdge) {
								possiblyVisibleFaces.emplace_back( m_mesh.m_halfEdges[m_mesh.m_halfEdges[heIndex].m_opp].m_face,heIndex );
							}
						}
						continue;
					}
					assert(faceData.m_faceIndex != topFaceIndex);
				}

				// The face is not visible. Therefore, the halfedge we came from is part of the horizon edge.
				pvf.m_isVisibleFaceOnCurrentIteration = 0;
				horizonEdges.push_back(faceData.m_enteredFromHalfEdge);
				// Store which half edge is the horizon edge. The other half edges of the face will not be part of the final mesh so their data slots can by recycled.
				std::int8_t ind = -1;
				const auto halfEdges = m_mesh.getHalfEdgeIndicesOfFace(m_mesh.m_faces[m_mesh.m_halfEdges[faceData.m_enteredFromHalfEdge].m_face]);
				for (std::int8_t k=0;k<3;k++) {
					if (halfEdges[k] == faceData.m_enteredFromHalfEdge) {
						ind = k;
						break;
					}
				}
				assert(ind>=0);
				m_mesh.m_faces[m_mesh.m_halfEdges[faceData.m_enteredFromHalfEdge].m_face].m_horizonEdgesOnCurrentIteration |= (1<<ind);
			}
			const size_t horizonEdgeCount = horizonEdges.size();

			// Order horizon edges so that they form a loop. This may fail due to numerical instability in which case we give up trying to solve horizon edge for this point and accept a minor degeneration in the convex hull.
			if (!reorderHorizonEdges(horizonEdges)) {
				std::cerr << "Failed to solve horizon edge." << std::endl;
				auto it = std::find(tf.m_pointsOnPositiveSide->begin(),tf.m_pointsOnPositiveSide->end(),activePointIndex);
				tf.m_pointsOnPositiveSide->erase(it);
				continue;
			}

			// Except for the horizon edges, all half edges of the visible faces can be marked as disabled. Their data slots will be reused.
			// The faces will be disabled as well, but we need to remember the points that were on the positive side of them - therefore
			// we save pointers to them.
			m_newFaceIndices.clear();
			m_newHalfEdgeIndices.clear();
			m_disabledFacePointVectors.clear();
			size_t disableCounter = 0;
			for (auto faceIndex : visibleFaces) {
				auto& disabledFace = m_mesh.m_faces[faceIndex];
				auto halfEdges = m_mesh.getHalfEdgeIndicesOfFace(disabledFace);
				for (size_t j=0;j<3;j++) {
					if ((disabledFace.m_horizonEdgesOnCurrentIteration & (1<<j)) == 0) {
						if (disableCounter < horizonEdgeCount*2) {
							// Use on this iteration
							m_newHalfEdgeIndices.push_back(halfEdges[j]);
							disableCounter++;
						}
						else {
							// Mark for reusal on later iteration step
							m_mesh.disableHalfEdge(halfEdges[j]);
						}
					}
				}
				// Disabled the face, but retain pointer to the points that were on the positive side of it. We need to assign those points
				// to the new faces we create shortly.
				m_disabledFacePointVectors.push_back(std::move(m_mesh.disableFace(faceIndex)));
			}
			if (disableCounter < horizonEdgeCount*2) {
				const size_t newHalfEdgesNeeded = horizonEdgeCount*2-disableCounter;
				for (size_t i=0;i<newHalfEdgesNeeded;i++) {
					m_newHalfEdgeIndices.push_back(m_mesh.addHalfEdge());
				}
			}
			
			// Create new faces using the edgeloop
			for (size_t i = 0; i < horizonEdgeCount; i++) {
				const IndexType AB = horizonEdges[i];

				auto horizonEdgeVertexIndices = m_mesh.getVertexIndicesOfHalfEdge(m_mesh.m_halfEdges[AB]);
				IndexType A,B,C;
				A = horizonEdgeVertexIndices[0];
				B = horizonEdgeVertexIndices[1];
				C = activePointIndex;

				const IndexType newFaceIndex = m_mesh.addFace();
				m_newFaceIndices.push_back(newFaceIndex);

				const IndexType CA = m_newHalfEdgeIndices[2*i+0];
				const IndexType BC = m_newHalfEdgeIndices[2*i+1];

				m_mesh.m_halfEdges[AB].m_next = BC;
				m_mesh.m_halfEdges[BC].m_next = CA;
				m_mesh.m_halfEdges[CA].m_next = AB;

				m_mesh.m_halfEdges[BC].m_face = newFaceIndex;
				m_mesh.m_halfEdges[CA].m_face = newFaceIndex;
				m_mesh.m_halfEdges[AB].m_face = newFaceIndex;

				m_mesh.m_halfEdges[CA].m_endVertex = A;
				m_mesh.m_halfEdges[BC].m_endVertex = C;

				auto& newFace = m_mesh.m_faces[newFaceIndex];

				const Vector3<T> planeNormal = mathutils::getTriangleNormal(pointCloud[A],pointCloud[B],activePoint);
				newFace.m_P = Plane<T>(planeNormal,activePoint);
				newFace.m_he = AB;

				m_mesh.m_halfEdges[CA].m_opp = m_newHalfEdgeIndices[i>0 ? i*2-1 : 2*horizonEdgeCount-1];
				m_mesh.m_halfEdges[BC].m_opp = m_newHalfEdgeIndices[((i+1)*2) % (horizonEdgeCount*2)];
			}

			// Assign points that were on the positive side of the disabled faces to the new faces.
			for (auto& disabledPoints : m_disabledFacePointVectors) {
				if (!disabledPoints)
					continue;
				for (const auto& point : *(disabledPoints)) {
					if (point == activePointIndex) {
						continue;
					}
					for (size_t j=0;j<horizonEdgeCount;j++) {
						auto& newFace = m_mesh.m_faces[m_newFaceIndices[j]];
						const T D = mathutils::getSignedDistanceToPlane(pointCloud[ point ],newFace.m_P);
						if (D>0 && D*D > epsilonSquared*newFace.m_P.m_sqrNLength) {
							if (!newFace.m_pointsOnPositiveSide) {
								newFace.m_pointsOnPositiveSide = std::move(m_mesh.getIndexVectorFromPool());
							}
							newFace.m_pointsOnPositiveSide->push_back( point );
							if (D > newFace.m_mostDistantPointDist) {
								newFace.m_mostDistantPointDist = D;
								newFace.m_mostDistantPoint = point;
							}
							break;
						}
					}
				}
				// The points are no longer needed: we can move them to the vector pool for reuse.
				m_mesh.reclaimToIndexVectorPool(std::move(disabledPoints));
			}

			// Increase face stack size if needed
			for (auto newFaceIndex : m_newFaceIndices) {
				auto& newFace = m_mesh.m_faces[newFaceIndex];
				if (newFace.m_pointsOnPositiveSide && newFace.m_pointsOnPositiveSide->size()) {
					if (!newFace.m_inFaceStack) {
						faceStack.push_back(newFaceIndex);
						newFace.m_inFaceStack = 1;
					}
				}
			}
		}
		
	}
	
	/*
	 * Implementations of degenerate case handler functions
	 */
	 
	template<typename T>
	ConvexHull<T> QuickHull<T>::checkDegenerateCase0D() {
		const std::vector<Vector3<T>>& pointCloud = *m_vertexData;

		// 0D degenerate case: all points are at the same location
		const Vector3<T>& v0 = *pointCloud.begin();
		for (const auto& v : pointCloud) {
			T d = (v-v0).getLengthSquared();
			if (d>m_epsilon*m_epsilon) {
				return ConvexHull<T>();
			}
		}
#ifdef DEBUG
		std::cout << "Detected 0D degenerate case: all points are at " << pointCloud[0] << "\n";
#endif

		ConvexHull<T> simpleHull;
		auto& indices = simpleHull.getIndexBuffer();
		auto& vertices = simpleHull.getVertexBuffer();
		indices.push_back(0);
		indices.push_back(std::min((size_t)1,pointCloud.size()-1));
		indices.push_back(std::min((size_t)2,pointCloud.size()-1));
		vertices.push_back(pointCloud[indices[0]]);
		vertices.push_back(pointCloud[indices[1]]);
		vertices.push_back(pointCloud[indices[2]]);
		return simpleHull;
	}

	template <typename T>
	ConvexHull<T> QuickHull<T>::checkDegenerateCase1D() {
		const std::vector<Vector3<T>>& vertices = *m_vertexData;

		// 1D degenerate case: the points form a line in 3D space. To keep things simple, we translate the points so that the first point resides at the origin. This way, if the points do actually form a line, we have a line passing through the origin.
		const Vector3<T>* firstPoint = nullptr;
		T firstPointLengthSquared=0;
		const Vector3<T>* maxPoint = nullptr;
		const Vector3<T>* minPoint = nullptr;
		const T epsilonSquared = m_epsilon*m_epsilon;
		T maxDot = -1.0f;
		T minDot = 1.0f;

		// First find a point which does not reside at the origin. Such a point exists, for otherwise we would have the 0D case.
		for (const auto& v : vertices) {
			const auto v2 = v-vertices[0];
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
		for (const auto& v : vertices) {
			const auto v2 = v-vertices[0];
			const T dot = v2.dotProduct(*firstPoint-vertices[0]);
			if (v2.getLengthSquared()>epsilonSquared) {
				const T V = dot*dot/(v2.getLengthSquared()*firstPointLengthSquared);
				const T d = std::abs(V-1);
				if (d > Epsilon) {
					// The points do not form a line!
					return ConvexHull<T>();
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
		for (const auto& v : vertices) {
			if (&v != minPoint && &v != maxPoint) {
				thirdPoint = &v;
				break;
			}
		}
		if (thirdPoint==nullptr) {
			thirdPoint = minPoint;
		}
		ConvexHull<T> hull;
		auto& indexBuffer = hull.getIndexBuffer();
		auto& vertexBuffer = hull.getVertexBuffer();
		vertexBuffer = {*minPoint,*thirdPoint,*maxPoint};
		indexBuffer = {0,1,2};
		return hull;
	}

	template<typename T>
	ConvexHull<T> QuickHull<T>::checkDegenerateCase2D() {
		const std::vector<Vector3<T>>& pointCloud = *m_vertexData;

		// 2D degenerate case: all points lie on the same plane. Just like in the 1D case, we translate the points so that the first point is located at the origin.
		const T epsilonSquared = m_epsilon*m_epsilon;

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
			if (d > m_epsilon) {
				// We have a proper 3D point cloud and the QuickHull algorithm can be applied.
				return ConvexHull<T>();
			}
		}

		// We have encountered the degenerate 2D case. We solve the problem by adding one extra point to the cloud so that we have a shape with volume, then applying the 3D QuickHull, and finally removing faces connected to the extra point.
		// TODO: implement proper 2D QuickHull, because this solution is not good performance-wise.
#ifdef DEBUG
		std::cout << "Degenerate 2D case detected." << std::endl;
#endif
		auto newPoints = pointCloud;
		const T M = m_epsilon/Epsilon;
		auto extraPoint = M*cross + pointCloud[0];
		newPoints.push_back(extraPoint);

		m_vertexData = &newPoints;
		createConvexHalfEdgeMesh();

		std::vector<size_t> disableList;
		for (size_t i = 0; i< m_mesh.m_faces.size();i++) {
			auto& face = m_mesh.m_faces[i];
			if (!face.isDisabled()) {
				auto v = m_mesh.getVertexIndicesOfFace(face);
				if ((newPoints[v[0]]-extraPoint).getLengthSquared()<M/4 || (newPoints[v[1]]-extraPoint).getLengthSquared()<M/4 || (newPoints[v[2]]-extraPoint).getLengthSquared()<M/4) {
					disableList.push_back(i);
				}
			}
		}
		for (auto i: disableList) {
			m_mesh.m_faces[i].disable();
		}
		return ConvexHull<T>(m_mesh,newPoints, true);
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

	template<typename T>
	bool QuickHull<T>::reorderHorizonEdges(std::vector<IndexType>& horizonEdges) {
		const size_t horizonEdgeCount = horizonEdges.size();
		for (size_t i=0;i<horizonEdgeCount-1;i++) {
			IndexType endVertex = m_mesh.m_halfEdges[ horizonEdges[i] ].m_endVertex;
			bool foundNext = false;
			for (size_t j=i+1;j<horizonEdgeCount;j++) {
				const IndexType beginVertex = m_mesh.m_halfEdges[ m_mesh.m_halfEdges[horizonEdges[j]].m_opp ].m_endVertex;
				if (beginVertex == endVertex) {
					std::swap(horizonEdges[i+1],horizonEdges[j]);
					foundNext = true;
					break;
				}
			}
			if (!foundNext) {
				return false;
			}
		}
		assert(m_mesh.m_halfEdges[ horizonEdges[horizonEdges.size()-1] ].m_endVertex == m_mesh.m_halfEdges[ m_mesh.m_halfEdges[horizonEdges[0]].m_opp ].m_endVertex);
		return true;
	}

	template <typename T>
	Mesh<T> QuickHull<T>::getInitialTetrahedron() {
		const auto& vertices = *m_vertexData;
		const T epsilonSquared = m_epsilon*m_epsilon;

		// Find two most distant extreme points.
		T maxD = 0.0f;
		std::pair<IndexType,IndexType> selectedPoints;
		for (size_t i=0;i<6;i++) {
			for (size_t j=i+1;j<6;j++) {
				const T d = vertices[ m_extremeValues[i] ].getSquaredDistanceTo( vertices[ m_extremeValues[j] ] );
				if (d > maxD) {
					maxD=d;
					selectedPoints=std::pair<IndexType,IndexType>(m_extremeValues[i],m_extremeValues[j]);
				}
			}
		}
		assert(maxD > 0.0);
		
		// Find the most distant point to the line between the two chosen extreme points.
		const Ray<T> r(vertices[selectedPoints.first], (vertices[selectedPoints.second] - vertices[selectedPoints.first]));
		maxD=0.0f;
		size_t maxI=std::numeric_limits<size_t>::max();
		const size_t vCount = vertices.size();
		for (size_t i=0;i<vCount;i++) {
			const T distToRay = mathutils::getSquaredDistanceBetweenPointAndRay(vertices[i],r);
			if (distToRay > maxD) {
				maxD=distToRay;
				maxI=i;
			}
		}
		assert(maxI!=std::numeric_limits<IndexType>::max());

		// These three points form the base triangle for our tetrahedron.
		std::array<IndexType,3> baseTriangle{selectedPoints.first, selectedPoints.second, maxI};
		const Vector3<T> baseTriangleVertices[]={ vertices[baseTriangle[0]], vertices[baseTriangle[1]],  vertices[baseTriangle[2]] };
		
		// Next step is to find the 4th vertex of the tetrahedron. We naturally choose the point farthest away from the triangle plane.
		maxD=0.0f;
		maxI = 0;
		const Vector3<T> N = mathutils::getTriangleNormal(baseTriangleVertices[0],baseTriangleVertices[1],baseTriangleVertices[2]);
		Plane<T> trianglePlane(N,baseTriangleVertices[0]);
		for (size_t i=0;i<vCount;i++) {
			const T d = std::abs(mathutils::getSignedDistanceToPlane(vertices[i],trianglePlane));
			if (d>maxD) {
				maxD=d;
				maxI=i;
			}
		}

		// Now that we have the 4th point, we can create the tetrahedron
		const Plane<T> triPlane(N,baseTriangleVertices[0]);
		if (triPlane.isPointOnPositiveSide(vertices[maxI])) {
			// Enforce CCW orientation (if user prefers clockwise orientation, swap two vertices in each triangle when final mesh created)
			std::swap(baseTriangle[0],baseTriangle[1]);
		}

		// Create a tetrahedron half edge mesh and compute planes defined by each triangle
		Mesh<T> mesh(baseTriangle[0],baseTriangle[1],baseTriangle[2],maxI);
		for (auto& f : mesh.m_faces) {
			auto v = mesh.getVertexIndicesOfFace(f);
			const Vector3<T>& va = vertices[v[0]];
			const Vector3<T>& vb = vertices[v[1]];
			const Vector3<T>& vc = vertices[v[2]];
			const Vector3<T> N = mathutils::getTriangleNormal(va, vb, vc);
			const Plane<T> trianglePlane(N,va);
			f.m_P = trianglePlane;
		}

		// Finally we assign a face for each vertex outside the tetrahedron
		for (size_t i=0;i<vCount;i++) {
			for (auto& face : mesh.m_faces) {
				const T D = mathutils::getSignedDistanceToPlane(vertices[i],face.m_P);
				if (D>0 && D*D > epsilonSquared*face.m_P.m_sqrNLength) {
					if (!face.m_pointsOnPositiveSide) {
						face.m_pointsOnPositiveSide.reset(new std::vector<IndexType>());
					}
					face.m_pointsOnPositiveSide->push_back(i);
					if (D > face.m_mostDistantPointDist) {
						face.m_mostDistantPointDist = D;
						face.m_mostDistantPoint = i;
					}
					break;
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

