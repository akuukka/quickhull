/*
 * ConvexHull.hpp
 *
 *  Created on: Sep 19, 2015
 *      Author: anttiku
 */

#ifndef CONVEXHULL_HPP_
#define CONVEXHULL_HPP_

#include "Structs/Vector3.hpp"
#include "Structs/HalfEdgeMesh.hpp"
#include <vector>
#include <map>
#include <fstream>

namespace quickhull {

	template<typename T>
	class ConvexHull {
		std::vector<Vector3<T>> m_vertices;
		std::vector<size_t> m_indices;
	public:
		ConvexHull() {}

		// Construct vertex and index buffers from half edge mesh and pointcloud
		ConvexHull(const Mesh<T>& mesh, const std::vector<Vector3<T>>& pointCloud, bool CCW) {
			std::vector<bool> faceProcessed(mesh.m_faces.size(),false);
			std::vector<size_t> faceStack;
			std::map<size_t,size_t> vertexIndexMapping; // Map vertex indices from original point cloud to the new mesh vertex indices
			for (size_t i = 0;i<mesh.m_faces.size();i++) {
				if (!mesh.m_faces[i].isDisabled()) {
					faceStack.push_back(i);
					break;
				}
			}

			if (faceStack.size()==0) {
				return;
			}

			while (faceStack.size()) {
				auto it = faceStack.end()-1;
				size_t top = *it;
				assert(!mesh.m_faces[top].isDisabled());
				faceStack.erase(it);
				if (faceProcessed[top]) {
					continue;
				}
				else {
					auto halfEdges = mesh.getHalfEdgeIndicesOfFace(mesh.m_faces[top]);
					size_t adjacent[] = {mesh.m_halfEdges[mesh.m_halfEdges[halfEdges[0]].m_opp].m_face,mesh.m_halfEdges[mesh.m_halfEdges[halfEdges[1]].m_opp].m_face,mesh.m_halfEdges[mesh.m_halfEdges[halfEdges[2]].m_opp].m_face};
					for (auto a : adjacent) {
						if (!faceProcessed[a] && !mesh.m_faces[a].isDisabled()) {
							faceStack.push_back(a);
						}
					}
					auto vertices = mesh.getVertexIndicesOfFace(mesh.m_faces[top]);
					for (auto& v : vertices) {
						auto it = vertexIndexMapping.find(v);
						if (it == vertexIndexMapping.end()) {
							m_vertices.push_back(pointCloud[v]);
							v = m_vertices.size()-1;
						}
						else {
							v = it->second;
						}
					}
					m_indices.push_back(vertices[0]);
					if (CCW) {
						m_indices.push_back(vertices[2]);
						m_indices.push_back(vertices[1]);
					}
					else {
						m_indices.push_back(vertices[1]);
						m_indices.push_back(vertices[2]);
					}
					faceProcessed[top]=true;
				}
			}
		}

		std::vector<size_t>& getIndexBuffer() {
			return m_indices;
		}

		std::vector<Vector3<T>>& getVertexBuffer() {
			return m_vertices;
		}
		
		// Export the mesh to a Waveform OBJ file
		void writeWaveformOBJ(const std::string& filename)
		{
			std::ofstream objFile;
			objFile.open (filename);
			objFile << "o quickhull\n";
			for (const auto& v : getVertexBuffer()) {
				objFile << "v " << v.x << " " << v.y << " " << v.z << "\n";
			}
			const auto& indBuf = getIndexBuffer();
			size_t triangleCount = indBuf.size()/3;
			for (int i=0;i<triangleCount;i++) {
				objFile << "f " << indBuf[i*3]+1 << " " << indBuf[i*3+1]+1 << " " << indBuf[i*3+2]+1 << "\n";
			}
			objFile.close();
		}

	};

}

#endif /* CONVEXHULL_HPP_ */
