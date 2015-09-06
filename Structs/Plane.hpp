/*
 * Plane.hpp
 *
 *  Created on: Dec 7, 2013
 *      Author: anttiku
 */

#ifndef QHPLANE_HPP_
#define QHPLANE_HPP_

#include "Vector3.hpp"

namespace quickhull {

	template<typename T>
	class Plane {
	public:
		// Unit normal of the plane
		Vector3<T> m_N;
		
		// Signed distance to the origin from the plane
		T m_D;

		bool isPointOnPositiveSide(const Vector3<T>& Q) const {
			T d = m_N.dotProduct(Q)+m_D;
			if (d>=0) return true;
			return false;
		}

		void setNormal(const Vector3<T>& N) {
			m_N = N.getNormalized();
		}

		const Vector3<T>& getNormal() const {
			return m_N;
		}
		
		T getSignedDistance() const {
			return m_D;
		}

		Plane() {

		}

		// Construct a plane with unit normal N and signed distance D to the origin
		Plane(const Vector3<T>& N, T D) : m_N(N), m_D(D) {
			
		}

		// Construct a plane from unit normal N and any point P on the plane
		Plane(const Vector3<T>& N, const Vector3<T>& P) : m_N(N), m_D(-N.dotProduct(P)) {
			
		}
	};

}


#endif /* PLANE_HPP_ */
