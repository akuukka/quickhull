/*
 * HalfEdgeMesh.hpp
 *
 *  Created on: Sep 7, 2015
 *      Author: anttiku
 */

#ifndef HALFEDGEMESH_HPP_
#define HALFEDGEMESH_HPP_

#include <vector>
#include "Vector3.hpp"
#include "Plane.hpp"
#include "../Types.hpp"
#include <string>
#include <array>
#include <cassert>
#include <limits>

namespace quickhull {

	template <typename T>
	class Mesh {
	public:
		struct HalfEdge {
			IndexType m_endVertex;
			IndexType m_opp;
			IndexType m_face;
			IndexType m_next;
		};

		struct Face {
			IndexType m_he;
			Plane<T> m_P;
			std::vector<IndexType> m_pointsOnPositiveSide;
			T m_mostDistantPointDist;
			IndexType m_mostDistantPoint;
			size_t m_visibilityCheckedOnIteration;
			std::uint8_t m_isVisibleFaceOnCurrentIteration : 1;
			std::uint8_t m_inFaceStack : 1;
			std::uint8_t m_horizonEdgesOnCurrentIteration : 3; // Bit for each half edge assigned to this face, each being 0 or 1 depending on whether the edge belong to horizon edge

			Face() : m_he(std::numeric_limits<IndexType>::max()),
					 m_mostDistantPointDist(0),
					 m_mostDistantPoint(0),
					 m_visibilityCheckedOnIteration(0),
					 m_isVisibleFaceOnCurrentIteration(0),
					 m_inFaceStack(0),
					 m_horizonEdgesOnCurrentIteration(0)
			{

			}

			void disable() {
				m_he = std::numeric_limits<IndexType>::max();
			}

			bool isDisabled() const {
				return m_he == std::numeric_limits<IndexType>::max();
			}
		};

		std::vector<Face> m_faces;
		std::vector<HalfEdge> m_halfEdges;
		std::vector<IndexType> m_disabledFaces,m_disabledHalfEdges;

		IndexType addFace() {
			if (m_disabledFaces.size()) {
				auto it = m_disabledFaces.end()-1;
				IndexType index = *it;
				auto& f = m_faces[index];
				assert(f.isDisabled());
				if ((f.m_pointsOnPositiveSide.size()+1)*256 < f.m_pointsOnPositiveSide.capacity()) {
					// Reduce memory usage! Huge vectors are needed at the beginning of iteration when faces have many points on their positive side. Since we recycle the face objects, vector buffers should be reset when they are clearly too large.
					// There is performance penalty, naturally, so we don't want to do this too often.
					f.m_pointsOnPositiveSide = std::vector<IndexType>();
				}
				else {
					f.m_pointsOnPositiveSide.clear();
				}
				f.m_mostDistantPointDist = 0;
				m_disabledFaces.erase(it);
				return index;
			}
			m_faces.push_back(Face());
			return m_faces.size()-1;
		}

		IndexType addHalfEdge()	{
			if (m_disabledHalfEdges.size()) {
				auto it = m_disabledHalfEdges.end()-1;
				IndexType index = *it;
				m_disabledHalfEdges.erase(it);
				return index;
			}
			m_halfEdges.push_back(HalfEdge());
			return m_halfEdges.size()-1;
		}

		void disableFace(IndexType faceIndex) {
			m_faces[faceIndex].disable();
			m_disabledFaces.push_back(faceIndex);
		}

		void disableHalfEdge(IndexType heIndex) {
			m_disabledHalfEdges.push_back(heIndex);
		}

		Mesh() {}

		// Create a mesh with initial tetrahedron ABCD. Dot product of AB with the normal of triangle ABC should be negative.
		Mesh(IndexType a, IndexType b, IndexType c, IndexType d) {
			// Create halfedges
			HalfEdge AB;
			AB.m_endVertex = b;
			AB.m_opp = 6;
			AB.m_face = 0;
			AB.m_next = 1;
			m_halfEdges.push_back(AB);

			HalfEdge BC;
			BC.m_endVertex = c;
			BC.m_opp = 9;
			BC.m_face = 0;
			BC.m_next = 2;
			m_halfEdges.push_back(BC);

			HalfEdge CA;
			CA.m_endVertex = a;
			CA.m_opp = 3;
			CA.m_face = 0;
			CA.m_next = 0;
			m_halfEdges.push_back(CA);

			HalfEdge AC;
			AC.m_endVertex = c;
			AC.m_opp = 2;
			AC.m_face = 1;
			AC.m_next = 4;
			m_halfEdges.push_back(AC);

			HalfEdge CD;
			CD.m_endVertex = d;
			CD.m_opp = 11;
			CD.m_face = 1;
			CD.m_next = 5;
			m_halfEdges.push_back(CD);

			HalfEdge DA;
			DA.m_endVertex = a;
			DA.m_opp = 7;
			DA.m_face = 1;
			DA.m_next = 3;
			m_halfEdges.push_back(DA);

			HalfEdge BA;
			BA.m_endVertex = a;
			BA.m_opp = 0;
			BA.m_face = 2;
			BA.m_next = 7;
			m_halfEdges.push_back(BA);

			HalfEdge AD;
			AD.m_endVertex = d;
			AD.m_opp = 5;
			AD.m_face = 2;
			AD.m_next = 8;
			m_halfEdges.push_back(AD);

			HalfEdge DB;
			DB.m_endVertex = b;
			DB.m_opp = 10;
			DB.m_face = 2;
			DB.m_next = 6;
			m_halfEdges.push_back(DB);

			HalfEdge CB;
			CB.m_endVertex = b;
			CB.m_opp = 1;
			CB.m_face = 3;
			CB.m_next = 10;
			m_halfEdges.push_back(CB);

			HalfEdge BD;
			BD.m_endVertex = d;
			BD.m_opp = 8;
			BD.m_face = 3;
			BD.m_next = 11;
			m_halfEdges.push_back(BD);

			HalfEdge DC;
			DC.m_endVertex = c;
			DC.m_opp = 4;
			DC.m_face = 3;
			DC.m_next = 9;
			m_halfEdges.push_back(DC);

			// Create faces
			Face ABC;
			ABC.m_he = 0;
			m_faces.push_back(ABC);

			Face ACD;
			ACD.m_he = 3;
			m_faces.push_back(ACD);

			Face BAD;
			BAD.m_he = 6;
			m_faces.push_back(BAD);

			Face CBD;
			CBD.m_he = 9;
			m_faces.push_back(CBD);
		}

		std::array<IndexType,3> getVertexIndicesOfFace(const Face& f) const {
			std::array<IndexType,3> v;
			const HalfEdge* he = &m_halfEdges[f.m_he];
			v[0] = he->m_endVertex;
			he = &m_halfEdges[he->m_next];
			v[1] = he->m_endVertex;
			he = &m_halfEdges[he->m_next];
			v[2] = he->m_endVertex;
			return v;
		}

		std::array<IndexType,2> getVertexIndicesOfHalfEdge(const HalfEdge& he) const {
			return {m_halfEdges[he.m_opp].m_endVertex,he.m_endVertex};
		}

		std::array<IndexType,3> getHalfEdgeIndicesOfFace(const Face& f) const {
			return {f.m_he,m_halfEdges[f.m_he].m_next,m_halfEdges[m_halfEdges[f.m_he].m_next].m_next};
		}
	};

}



#endif /* HALFEDGEMESH_HPP_ */
