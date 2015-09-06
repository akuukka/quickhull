#ifndef QuickHull_Face_hpp
#define QuickHull_Face_hpp

#include "Vector3.hpp"
#include "Plane.hpp"
#include "../Types.hpp"

namespace quickhull {
	
	template<typename T>
	struct Face {
		IndexType a,b,c; // Indices of the triangle vertices
		
		std::vector<IndexType> pointsOnPositiveSide; // Points on the "positive side" of the plane defined by this triangle
		Plane<T> P;
		
		Face(IndexType a, IndexType b, IndexType c) : a(a),b(b),c(c) {
		}
		
		Face(IndexType a, IndexType b, IndexType c, const Plane<T>& P) : a(a),b(b),c(c),P(P) {
		}
		
		Face() {
			a=-1;
		}
		
		void disable() {
			a=-1;
		}
		
		bool isDisabled() const {
			return a==-1;
		}
		
		bool containsPoint(IndexType index) const {
			if (index==a) return true;
			if (index==b) return true;
			if (index==c) return true;
			return false;
		}
	};

	
}

#endif
