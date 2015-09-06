
#ifndef QuickHull_MathUtils_hpp
#define QuickHull_MathUtils_hpp

#include "Structs/Vector3.hpp"
#include "Structs/Ray.hpp"

namespace quickhull {
	
	namespace mathutils {
		
		template <typename T>
		inline T getSquaredDistanceBetweenPointAndRay(const Vector3<T>& p, const Ray<T>& r) {
			T t = (p-r.m_S).dotProduct(r.m_V);
			return (p-r.m_S).getLengthSquared() - t*t;
		}
		
		template <typename T>
		inline T getSignedDistanceToPlane(const Vector3<T>& v, const Plane<T>& p) {
			return p.m_N.dotProduct(v) + p.m_D;
		}
		
	}
	
}


#endif
