#ifndef QuickHull_EdgeLoop_hpp
#define QuickHull_EdgeLoop_hpp

#include <vector>
#include "Vector3.hpp"
#include "Face.hpp"
#include <algorithm>

namespace quickhull {
	
	/*
	 
	  This helper struct is used to find edge loop for a list of connected triangles. The following picture should make this clear:
	 

	       A-----H
	      / \   / \
	     /   \ /   \
	    B-----C-----G    =======>   A,B,D,F,E,G,H
	     \   / \   /
	      \ /   \ /
	       D-----E
	        \   /
	         \ /
	          F
	 */
	
	template<typename T>
	struct EdgeLoop {
		std::vector<IndexType> m_vertices;
		
		EdgeLoop() {
			
		}
		
		// Construct an edgeloop from a vector of faces. If this fails (happens when the faces are not connected), vertex vector will be empty.
		void createFrom(std::vector<IndexType> facesToJoin, const std::vector<Face<T>>& faces) {
			m_vertices.clear();
			IndexType i=0;
			IndexType failCounter=0;
			
			while (facesToJoin.size()>0) {
				auto it = facesToJoin.end()-1-i;
				if (joinWith(faces[*it])) {
					facesToJoin.erase(it);
					failCounter=0;
					if (i>=facesToJoin.size()) {
						i=0;
					}
				}
				else {
					failCounter++;
					i++;
					if (i==facesToJoin.size()) {
						i=0;
					}
					if (failCounter==facesToJoin.size()) {
						// None of the remaining faces can be added to the edge loop. We must give up.
						m_vertices.clear();
						return;
					}
				}
			}
		}
		
		size_t getEdgeCount() const {
			return m_vertices.size();
		}
		
		std::pair<IndexType,IndexType> getEdge(IndexType index) const {
			return std::pair<IndexType,IndexType>(m_vertices[index],m_vertices[(index+1)%m_vertices.size()]);
		}
		
		// Try to add a face to the edge loop
		bool joinWith(const Face<T>& f) {
			if (m_vertices.size()==0) {
				// First face
				m_vertices.push_back(f.a);
				m_vertices.push_back(f.b);
				m_vertices.push_back(f.c);
				return true;
			}
			
			// Check which vertices of the face are already included in the edge loop
			bool notHaveA = true;
			bool notHaveB = true;
			bool notHaveC = true;
			for (auto v : m_vertices) {
				if (v == f.a) notHaveA = false;
				if (v == f.b) notHaveB = false;
				if (v == f.c) notHaveC = false;
			}
			
			// If the edge loop shares no vertices with the face, we of course can not add the vertices to the edge loop.
			if (notHaveA && notHaveB && notHaveC) {
				return false;
			}

			const size_t s = m_vertices.size();
			for (IndexType i=0;i<s;i++) {
				const IndexType v1=i;
				const IndexType v2=(i+1)%s;
				const IndexType v3=(i+2)%s;
				if ( (m_vertices[v1] == f.b && m_vertices[v2] == f.a) && notHaveC ) {
					m_vertices.insert(m_vertices.begin()+i+1,f.c);
					return true;
				}
				if ( (m_vertices[v1] == f.a && m_vertices[v2] == f.c) && notHaveB  ) {
					m_vertices.insert(m_vertices.begin()+i+1,f.b);
					return true;
				}
				if ( (m_vertices[v1] == f.c && m_vertices[v2] == f.b) && notHaveA ) {
					m_vertices.insert(m_vertices.begin()+i+1,f.a);
					return true;
				}
				if (f.a==m_vertices[v1] && f.c==m_vertices[v2] && f.b==m_vertices[v3]) {
					m_vertices.erase(m_vertices.begin()+v2);
					return true;
				}
				else if (f.b==m_vertices[v1] && f.a==m_vertices[v2] && f.c==m_vertices[v3]) {
					m_vertices.erase(m_vertices.begin()+v2);
					return true;
				}
				else if (f.c==m_vertices[v1] && f.b==m_vertices[v2] && f.a==m_vertices[v3]) {
					m_vertices.erase(m_vertices.begin()+v2);
					return true;
				}
			}
			return false;
		}
	};
	
}

#endif
