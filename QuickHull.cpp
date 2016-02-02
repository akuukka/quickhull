#include "QuickHull.hpp"
#include "Structs/Plane.hpp"
#include "MathUtils.hpp"
#include <cmath>
#include <cassert>
#include <functional>
#include <iostream>
#include <algorithm>
#include <deque>
#include <limits>
#include "Structs/Mesh.hpp"

namespace quickhull {
	
	template<>
	const float QuickHull<float>::Epsilon = 0.00001f;
	
	template<>
	const double QuickHull<double>::Epsilon = 0.0000001;
	
	/*
	 * Implementation of the algorithm
	 */
	
	template<typename T>
	ConvexHull<T> QuickHull<T>::getConvexHull(const std::vector<Vector3<T>>& pointCloud, bool CCW, bool useOriginalIndices, T epsilon) {
		VertexDataSource<T> vertexDataSource(pointCloud);
		return getConvexHull(vertexDataSource,CCW,useOriginalIndices,epsilon);
	}
	
	template<typename T>
	ConvexHull<T> QuickHull<T>::getConvexHull(const Vector3<T>* vertexData, size_t vertexCount, bool CCW, bool useOriginalIndices, T epsilon) {
		VertexDataSource<T> vertexDataSource(vertexData,vertexCount);
		return getConvexHull(vertexDataSource,CCW,useOriginalIndices,epsilon);
	}
	
	template<typename T>
	ConvexHull<T> QuickHull<T>::getConvexHull(const T* vertexData, size_t vertexCount, bool CCW, bool useOriginalIndices, T epsilon) {
		VertexDataSource<T> vertexDataSource((const Vector3<T>*)vertexData,vertexCount);
		return getConvexHull(vertexDataSource,CCW,useOriginalIndices,epsilon);
	}

	template<typename T>
	ConvexHull<T> QuickHull<T>::getConvexHull(const VertexDataSource<T>& pointCloud, bool CCW, bool useOriginalIndices, T epsilon) {
		if (pointCloud.size()==0) {
			return ConvexHull<T>();
		}
		m_vertexData = pointCloud;

		// Very first: find extreme values and use them to compute the scale of the point cloud.
		m_extremeValues = getExtremeValues();
		m_scale = getScale(m_extremeValues);
		
		// Epsilon we use depends on the scale
		m_epsilon = epsilon*m_scale;
		
		// Reset diagnostics
		m_diagnostics = DiagnosticsData();

		// Check for degenerate cases before proceeding to the 3D quickhull iteration phase
		std::function<ConvexHull<T>(bool)> degenerateCaseCheckers[] = {
			std::bind(&QuickHull<T>::checkDegenerateCase0D,this,std::placeholders::_1),
			std::bind(&QuickHull<T>::checkDegenerateCase1D,this,std::placeholders::_1),
			std::bind(&QuickHull<T>::checkDegenerateCase2D,this,std::placeholders::_1)
		};
		for (auto& f : degenerateCaseCheckers) {
			auto degenerateHull = f(useOriginalIndices);
			if (degenerateHull.getVertexBuffer().size()) {
				return degenerateHull;
			}
		}
		
		createConvexHalfEdgeMesh();
		return ConvexHull<T>(m_mesh,m_vertexData, CCW, useOriginalIndices);
	}

	template<typename T>
	void QuickHull<T>::createConvexHalfEdgeMesh() {
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
		std::deque<IndexType> faceList;
		for (size_t i=0;i < 4;i++) {
			auto& f = m_mesh.m_faces[i];
			if (f.m_pointsOnPositiveSide && f.m_pointsOnPositiveSide->size()>0) {
				faceList.push_back(i);
				f.m_inFaceStack = 1;
			}
		}
		if (faceList.empty()) {
			return;
		}

		// Process faces until the face list is empty.
		size_t iter = 0;
		while (!faceList.empty()) {
			iter++;
			if (iter == std::numeric_limits<size_t>::max()) {
				// Visible face traversal marks visited faces with iteration counter (to mark that the face has been visited on this iteration) and the max value represents unvisited faces. At this point we have to reset iteration counter. This shouldn't be an
				// issue on 64 bit machines.
				iter = 0;
			}
			
			auto faceIter = faceList.begin();
			const IndexType topFaceIndex = *faceIter;
			faceList.erase(faceIter);
			
			auto& tf = m_mesh.m_faces[topFaceIndex];
			tf.m_inFaceStack = 0;

			assert(!tf.m_pointsOnPositiveSide || tf.m_pointsOnPositiveSide->size()>0);
			if (!tf.m_pointsOnPositiveSide || tf.isDisabled()) {
				continue;
			}
			
			// Pick the most distant point to this triangle plane as the point to which we extrude
			const Vector3<T>& activePoint = m_vertexData[tf.m_mostDistantPoint];
			const size_t activePointIndex = tf.m_mostDistantPoint;

			// Find out the faces that have our active point on their positive side (these are the "visible faces"). The face on top of the stack of course is one of them. At the same time, we create a list of horizon edges.
			horizonEdges.clear();
			possiblyVisibleFaces.clear();
			visibleFaces.clear();
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
				const auto halfEdges = m_mesh.getHalfEdgeIndicesOfFace(m_mesh.m_faces[m_mesh.m_halfEdges[faceData.m_enteredFromHalfEdge].m_face]);
				const std::int8_t ind = (halfEdges[0]==faceData.m_enteredFromHalfEdge) ? 0 : (halfEdges[1]==faceData.m_enteredFromHalfEdge ? 1 : 2);
				m_mesh.m_faces[m_mesh.m_halfEdges[faceData.m_enteredFromHalfEdge].m_face].m_horizonEdgesOnCurrentIteration |= (1<<ind);
			}
			const size_t horizonEdgeCount = horizonEdges.size();

			// Order horizon edges so that they form a loop. This may fail due to numerical instability in which case we give up trying to solve horizon edge for this point and accept a minor degeneration in the convex hull.
			if (!reorderHorizonEdges(horizonEdges)) {
				m_diagnostics.m_failedHorizonEdges++;
				std::cerr << "Failed to solve horizon edge." << std::endl;
				auto it = std::find(tf.m_pointsOnPositiveSide->begin(),tf.m_pointsOnPositiveSide->end(),activePointIndex);
				tf.m_pointsOnPositiveSide->erase(it);
				if (tf.m_pointsOnPositiveSide->size()==0) {
					reclaimToIndexVectorPool(tf.m_pointsOnPositiveSide);
				}
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
				// Disable the face, but retain pointer to the points that were on the positive side of it. We need to assign those points
				// to the new faces we create shortly.
				auto t = std::move(m_mesh.disableFace(faceIndex));
				if (t) {
					assert(t->size()); // Because we should not assign point vectors to faces unless needed...
					m_disabledFacePointVectors.push_back(std::move(t));
				}
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

				const Vector3<T> planeNormal = mathutils::getTriangleNormal(m_vertexData[A],m_vertexData[B],activePoint);
				newFace.m_P = Plane<T>(planeNormal,activePoint);
				newFace.m_he = AB;

				m_mesh.m_halfEdges[CA].m_opp = m_newHalfEdgeIndices[i>0 ? i*2-1 : 2*horizonEdgeCount-1];
				m_mesh.m_halfEdges[BC].m_opp = m_newHalfEdgeIndices[((i+1)*2) % (horizonEdgeCount*2)];
			}

			// Assign points that were on the positive side of the disabled faces to the new faces.
			for (auto& disabledPoints : m_disabledFacePointVectors) {
				assert(disabledPoints);
				for (const auto& point : *(disabledPoints)) {
					if (point == activePointIndex) {
						continue;
					}
					for (size_t j=0;j<horizonEdgeCount;j++) {
						auto& newFace = m_mesh.m_faces[m_newFaceIndices[j]];
						const T D = mathutils::getSignedDistanceToPlane(m_vertexData[ point ],newFace.m_P);
						if (D>0 && D*D > epsilonSquared*newFace.m_P.m_sqrNLength) {
							if (!newFace.m_pointsOnPositiveSide) {
								newFace.m_pointsOnPositiveSide = std::move(getIndexVectorFromPool());
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
				reclaimToIndexVectorPool(disabledPoints);
			}

			// Increase face stack size if needed
			for (auto newFaceIndex : m_newFaceIndices) {
				auto& newFace = m_mesh.m_faces[newFaceIndex];
				if (newFace.m_pointsOnPositiveSide) {
					assert(newFace.m_pointsOnPositiveSide->size()>0);
					if (!newFace.m_inFaceStack) {
						faceList.push_back(newFaceIndex);
						newFace.m_inFaceStack = 1;
					}
				}
			}
		}
		
		// Cleanup
		m_indexVectorPool.clear();
	}
	
	/*
	 * Implementations of degenerate case handler functions
	 */
	 
	template<typename T>
	ConvexHull<T> QuickHull<T>::checkDegenerateCase0D(bool useOriginalIndices) {
		// 0D degenerate case: all points are at the same location
		const T epsilonSquared = m_epsilon*m_epsilon;
		const Vector3<T>& v0 = *m_vertexData.begin();
		const auto it = std::find_if(m_vertexData.begin(),m_vertexData.end(),[&](const Vector3<T>& v) { return ((v-v0).getLengthSquared())>epsilonSquared; });
		if (it!=m_vertexData.end()) {
			return ConvexHull<T>();
		}
#ifdef DEBUG
		std::cout << "Detected 0D degenerate case: all points are at " << m_vertexData[0] << "\n";
#endif
		return ConvexHull<T>({0,std::min((size_t)1,m_vertexData.size()),std::min((size_t)2,m_vertexData.size())}, m_vertexData, true, useOriginalIndices);
	}

	template <typename T>
	ConvexHull<T> QuickHull<T>::checkDegenerateCase1D(bool useOriginalIndices) {
		// 1D degenerate case: the points form a line in 3D space. To keep things simple, we translate the points so that the first point resides at the origin. This way, if the points do actually form a line, we have a line passing through the origin.
		size_t firstPoint = -1;
		T firstPointLengthSquared=0;
		size_t maxPoint = -1;
		size_t minPoint = -1;
		const T epsilonSquared = m_epsilon*m_epsilon;
		T maxDot = -1.0f;
		T minDot = 1.0f;

		// First find a point which does not reside at the origin. Such a point exists, for otherwise we would have the 0D case.
		for (size_t i=0;i<m_vertexData.size();i++) {
			const auto v = m_vertexData[i]-m_vertexData[0];
			if (firstPoint == -1) {
				if (v.getLengthSquared() > epsilonSquared) {
					firstPoint = i;
					firstPointLengthSquared = v.getLengthSquared();
					break;
				}
			}
		}
		assert(firstPoint != -1);

		// Then check if all other translated points point to the same direction as V=(m_vertexData[firstPoint]-m_vertexData[0]) - or lie at the origin. If not, proceed to checking the 2D degenerate case. While looping, keep track of min and max dot product because if the point cloud turns out to be a line, the min and max points are its end points.
		const auto V = m_vertexData[firstPoint]-m_vertexData[0];
		for (size_t i=0;i<m_vertexData.size();i++) {
			const auto v = m_vertexData[i]-m_vertexData[0];
			const auto proj = v.projection(V);
			const auto ortho = v - proj;
			const auto D = ortho.getLengthSquared();
			if (D > epsilonSquared) {
				// The points do not form a line!
				return ConvexHull<T>();
			}
			const T dot = v.dotProduct(V);
			if (maxPoint == -1 || dot > maxDot) {
				maxDot = dot;
				maxPoint = i;
			}
			if (minPoint == -1 || dot < minDot) {
				minDot = dot;
				minPoint = i;
			}
		}

		// We have a degenerate 1D case. Find a third point so that we can construct an infinitely thin triangle.
#ifdef DEBUG
		std::cout << "Detected 1D degenerate case: the point cloud forms a line between " << m_vertexData[minPoint] << " and " << m_vertexData[maxPoint] << std::endl;
#endif
		size_t thirdPoint = -1;
		for (size_t i=0;i<m_vertexData.size();i++) {
			if (i != minPoint && i != maxPoint) {
				thirdPoint = i;
				break;
			}
		}
		if (thirdPoint==-1) {
			thirdPoint = minPoint;
		}
		return ConvexHull<T>({minPoint,maxPoint,thirdPoint},m_vertexData,true,useOriginalIndices);
	}

	template<typename T>
	ConvexHull<T> QuickHull<T>::checkDegenerateCase2D(bool useOriginalIndices) {
		// 2D degenerate case: all points lie on the same plane. Just like in the 1D case, we translate the points so that the first point is located at the origin.
		const T epsilonSquared = m_epsilon*m_epsilon;
		
		// Find first point not lying at the origo (such must exist, because we have already checked for the 0D case)
		const Vector3<T>* firstPoint = nullptr;
		for (size_t i=1;i<m_vertexData.size();i++) {
			const auto v = m_vertexData[i]-m_vertexData[0];
			if (v.getLengthSquared() > epsilonSquared) {
				firstPoint = &m_vertexData[i];
				break;
			}
		}
		assert(firstPoint != nullptr);

		// Find two points not lying at the origin and not pointing to the same direction (there must be at least two, for otherwise the 0D or 1D cases would have generated the convex hull)
		const auto V = *firstPoint-m_vertexData[0];
		const Vector3<T>* secondPoint = nullptr;
		for (size_t i=1;i<m_vertexData.size();i++) {
			const auto v = m_vertexData[i]-m_vertexData[0];
			const auto proj = v.projection(V);
			const auto ortho = v - proj;
			const auto D = ortho.getLengthSquared();
			if (D > epsilonSquared) {
				// The points do not form a line!
				secondPoint = &m_vertexData[i];
				break;
			}
		}
		assert(firstPoint != nullptr && secondPoint != nullptr);
		// Now firstPoint and secondPoint define a plane. Its normal is their cross product.
		const auto N = mathutils::getTriangleNormal(*firstPoint,*secondPoint,m_vertexData[0]);
		const auto P = Plane<T>(N,m_vertexData[0]);
		const auto limit = epsilonSquared*P.m_sqrNLength;
		if (std::isinf(limit)) {
			throw std::runtime_error("Reached infinity. Your point cloud is too large.");
		}
		for (const auto& v : m_vertexData) {
			const T D = mathutils::getSignedDistanceToPlane(v,P);
			const auto A = D*D;
			if (std::isinf(A)) {
				throw std::runtime_error("Reached infinity. Your point cloud is too large.");
			}
			if (A > limit) {
				// We have a proper 3D point cloud and the QuickHull algorithm can be applied.
				return ConvexHull<T>();
			}
		}

		// We have encountered the degenerate 2D case. We solve the problem by adding one extra point to the cloud so that we have a shape with volume, then applying the 3D QuickHull, and finally removing faces connected to the extra point.
		// TODO: implement proper 2D QuickHull, because this solution is not good performance-wise.
#ifdef DEBUG
		std::cout << "Degenerate 2D case detected." << std::endl;
#endif
		std::vector<Vector3<T>> newPoints;
		newPoints.insert(newPoints.begin(),m_vertexData.begin(),m_vertexData.end());
		auto extraPoint = N + m_vertexData[0];
		newPoints.push_back(extraPoint);
		const auto extraPointIndex = newPoints.size()-1;

		const auto origVertexSource = m_vertexData;
		m_vertexData = VertexDataSource<T>(newPoints);
		createConvexHalfEdgeMesh();

		std::vector<size_t> disableList;
		for (size_t i = 0; i< m_mesh.m_faces.size();i++) {
			auto& face = m_mesh.m_faces[i];
			if (!face.isDisabled()) {
				auto v = m_mesh.getVertexIndicesOfFace(face);
				if (v[0]==extraPointIndex || v[1]==extraPointIndex || v[2]==extraPointIndex) {
					disableList.push_back(i);
				}
			}
		}
		for (auto i: disableList) {
			m_mesh.m_faces[i].disable();
		}
		return ConvexHull<T>(m_mesh,origVertexSource, true, useOriginalIndices);
	}
	
	/*
	 * Private helper functions
	 */

	template <typename T>
	std::array<IndexType,6> QuickHull<T>::getExtremeValues() {
		std::array<IndexType,6> outIndices{0,0,0,0,0,0};
		T extremeVals[6] = {m_vertexData[0].x,m_vertexData[0].x,m_vertexData[0].y,m_vertexData[0].y,m_vertexData[0].z,m_vertexData[0].z};
		const size_t vCount = m_vertexData.size();
		for (size_t i=1;i<vCount;i++) {
			const Vector3<T>& pos = m_vertexData[i];
			if (pos.x>extremeVals[0]) {
				extremeVals[0]=pos.x;
				outIndices[0]=(IndexType)i;
			}
			else if (pos.x<extremeVals[1]) {
				extremeVals[1]=pos.x;
				outIndices[1]=(IndexType)i;
			}
			if (pos.y>extremeVals[2]) {
				extremeVals[2]=pos.y;
				outIndices[2]=(IndexType)i;
			}
			else if (pos.y<extremeVals[3]) {
				extremeVals[3]=pos.y;
				outIndices[3]=(IndexType)i;
			}
			if (pos.z>extremeVals[4]) {
				extremeVals[4]=pos.z;
				outIndices[4]=(IndexType)i;
			}
			else if (pos.z<extremeVals[5]) {
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
			const IndexType endVertex = m_mesh.m_halfEdges[ horizonEdges[i] ].m_endVertex;
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
		const T epsilonSquared = m_epsilon*m_epsilon;

		// Find two most distant extreme points.
		T maxD = 0;
		std::pair<IndexType,IndexType> selectedPoints;
		for (size_t i=0;i<6;i++) {
			for (size_t j=i+1;j<6;j++) {
				const T d = m_vertexData[ m_extremeValues[i] ].getSquaredDistanceTo( m_vertexData[ m_extremeValues[j] ] );
				if (d > maxD) {
					maxD=d;
					selectedPoints={m_extremeValues[i],m_extremeValues[j]};
				}
			}
		}
		assert(maxD > 0);
		
		// Find the most distant point to the line between the two chosen extreme points.
		const Ray<T> r(m_vertexData[selectedPoints.first], (m_vertexData[selectedPoints.second] - m_vertexData[selectedPoints.first]));
		maxD=0;
		size_t maxI=std::numeric_limits<size_t>::max();
		const size_t vCount = m_vertexData.size();
		for (size_t i=0;i<vCount;i++) {
			const T distToRay = mathutils::getSquaredDistanceBetweenPointAndRay(m_vertexData[i],r);
			if (distToRay > maxD) {
				maxD=distToRay;
				maxI=i;
			}
		}
		assert(maxI!=std::numeric_limits<IndexType>::max());

		// These three points form the base triangle for our tetrahedron.
		std::array<IndexType,3> baseTriangle{selectedPoints.first, selectedPoints.second, maxI};
		const Vector3<T> baseTriangleVertices[]={ m_vertexData[baseTriangle[0]], m_vertexData[baseTriangle[1]],  m_vertexData[baseTriangle[2]] };
		
		// Next step is to find the 4th vertex of the tetrahedron. We naturally choose the point farthest away from the triangle plane.
		maxD=0;
		maxI=0;
		const Vector3<T> N = mathutils::getTriangleNormal(baseTriangleVertices[0],baseTriangleVertices[1],baseTriangleVertices[2]);
		Plane<T> trianglePlane(N,baseTriangleVertices[0]);
		for (size_t i=0;i<vCount;i++) {
			const T d = std::abs(mathutils::getSignedDistanceToPlane(m_vertexData[i],trianglePlane));
			if (d>maxD) {
				maxD=d;
				maxI=i;
			}
		}

		// Enforce CCW orientation (if user prefers clockwise orientation, swap two vertices in each triangle when final mesh is created)
		const Plane<T> triPlane(N,baseTriangleVertices[0]);
		if (triPlane.isPointOnPositiveSide(m_vertexData[maxI])) {
			std::swap(baseTriangle[0],baseTriangle[1]);
		}

		// Create a tetrahedron half edge mesh and compute planes defined by each triangle
		Mesh<T> mesh(baseTriangle[0],baseTriangle[1],baseTriangle[2],maxI);
		for (auto& f : mesh.m_faces) {
			auto v = mesh.getVertexIndicesOfFace(f);
			const Vector3<T>& va = m_vertexData[v[0]];
			const Vector3<T>& vb = m_vertexData[v[1]];
			const Vector3<T>& vc = m_vertexData[v[2]];
			const Vector3<T> N = mathutils::getTriangleNormal(va, vb, vc);
			const Plane<T> trianglePlane(N,va);
			f.m_P = trianglePlane;
		}

		// Finally we assign a face for each vertex outside the tetrahedron (vertices inside the tetrahedron have no role anymore)
		for (size_t i=0;i<vCount;i++) {
			for (auto& face : mesh.m_faces) {
				const T D = mathutils::getSignedDistanceToPlane(m_vertexData[i],face.m_P);
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

