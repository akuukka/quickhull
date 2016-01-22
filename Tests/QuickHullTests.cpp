//
//  QuickHullTests.cpp
//  QuickHull
//
//  Created by Antti Kuukka on 1/5/16.
//  Copyright © 2016 Antti Kuukka. All rights reserved.
//

#include "../QuickHull.hpp"
#include "../MathUtils.hpp"
#include <iostream>
#include <random>
#include <cassert>

namespace quickhull {
	
	namespace tests {
		
		using FloatType = float;
		
		void assertSameValue(FloatType a, FloatType b) {
			assert(std::abs(a-b)<0.0001f);
		}
		
		void testVector3() {
			typedef Vector3<FloatType> vec3;
			vec3 a(1,0,0);
			vec3 b(1,0,0);
			
			vec3 c = a.projection(b);
			assert( (c-a).getLength()<0.00001f);
			
			a = vec3(1,1,0);
			b = vec3(1,3,0);
			c = b.projection(a);
			assert( (c-vec3(2,2,0)).getLength()<0.00001f);
		}
		
		template <typename T>
		std::vector<Vector3<T>> createSphere(T radius, size_t M, Vector3<T> offset = Vector3<T>(0,0,0)) {
			std::vector<Vector3<T>> pc;
			const T pi = 3.14159f;
			for (int i=0;i<=M;i++) {
				FloatType y = std::sin(pi/2 + (FloatType)i/(M)*pi);
				FloatType r = std::cos(pi/2 + (FloatType)i/(M)*pi);
				FloatType K = FloatType(1)-std::abs((FloatType)((FloatType)i-M/2.0f))/(FloatType)(M/2.0f);
				const size_t pcount = (size_t)(1 + K*M + FloatType(1)/2);
				for (size_t j=0;j<pcount;j++) {
					FloatType x = pcount>1 ? r*std::cos((FloatType)j/pcount*pi*2) : 0;
					FloatType z = pcount>1 ? r*std::sin((FloatType)j/pcount*pi*2) : 0;
					pc.emplace_back(x+offset.x,y+offset.y,z+offset.z);
				}
			}
			return pc;
		}
		
		void sphereTest() {
			QuickHull<FloatType> qh;
			FloatType y = 1;
			for (;;) {
				auto pc = createSphere<FloatType>(1, 100, Vector3<FloatType>(0,y,0));
				auto hull = qh.getConvexHull(pc,true,false);
				y *= 15;
				y /= 10;
				if (hull.getVertexBuffer().size()==3) {
					break;
				}
			}
		}
		
		void testPlanes() {
			Vector3<FloatType> N(1,0,0);
			Vector3<FloatType> p(2,0,0);
			Plane<FloatType> P(N,p);
			auto dist = mathutils::getSignedDistanceToPlane(Vector3<FloatType>(3,0,0), P);
			assertSameValue(dist,1);
			dist = mathutils::getSignedDistanceToPlane(Vector3<FloatType>(1,0,0), P);
			assertSameValue(dist,-1);
			dist = mathutils::getSignedDistanceToPlane(Vector3<FloatType>(1,0,0), P);
			assertSameValue(dist,-1);
			N = Vector3<FloatType>(2,0,0);
			P = Plane<FloatType>(N,p);
			dist = mathutils::getSignedDistanceToPlane(Vector3<FloatType>(6,0,0), P);
			assertSameValue(dist,8);
		}
		
		void run() {
			// Setup test env
			const size_t N = 1000;
			std::vector<Vector3<FloatType>> pc;
			QuickHull<FloatType> qh;
			
			// Setup RNG
			std::mt19937 rng;
			std::uniform_real_distribution<> dist(0,1);
			auto rnd = [&](FloatType from, FloatType to) {
				return from + (FloatType)dist(rng)*(to-from);
			};
			
			// Test 1 : Put N points inside unit cube. Result mesh must have exactly 8 vertices because the convex hull is the unit cube.
			for (int i=0;i<8;i++) {
				pc.emplace_back(i&1 ? -1 : 1,i&2 ? -1 : 1,i&4 ? -1 : 1);
			}
			for (size_t i=0;i<N;i++)
			{
				pc.emplace_back(rnd(-1,1),rnd(-1,1),rnd(-1,1));
			}
			auto hull = qh.getConvexHull(pc,true,false);
			assert(hull.getVertexBuffer().size()==8);
			assert(hull.getIndexBuffer().size()==3*2*6); // 6 cube faces, 2 triangles per face, 3 indices per triangle
			assert(&(hull.getVertexBuffer()[0])!=&(pc[0]));
			auto hull2 = hull;
			assert(hull2.getVertexBuffer().size()==hull.getVertexBuffer().size());
			assert(hull2.getVertexBuffer()[0].x==hull.getVertexBuffer()[0].x);
			assert(hull2.getIndexBuffer().size()==hull.getIndexBuffer().size());
			auto hull3 = std::move(hull);
			assert(hull.getIndexBuffer().size()==0);
			
			// Test 1.1 : Same test, but using the original indices.
			hull = qh.getConvexHull(pc,true,true);
			assert(hull.getIndexBuffer().size()==3*2*6);
			assert(hull.getVertexBuffer().size()==pc.size());
			assert(&(hull.getVertexBuffer()[0])==&(pc[0]));
			
			// Test 2 : random N points from the boundary of unit sphere. Result mesh must have exactly N points.
			pc = createSphere<FloatType>(1, 450);
			hull = qh.getConvexHull(pc,true,false);
			assert(pc.size() == hull.getVertexBuffer().size());
			hull = qh.getConvexHull(pc,true,true);
			// Add every vertex twice. This should not affect final mesh
			auto s = pc.size();
			for (size_t i=0;i<s;i++) {
				const auto& p = pc[i];
				pc.push_back(p);
			}
			hull = qh.getConvexHull(pc,true,false);
			assert(pc.size()/2 == hull.getVertexBuffer().size());
			
			// Test 2.1 : Multiply x components of the unit sphere vectors by a huge number => essentially we get a line
			const FloatType mul = 2*2*2;
			while (true) {
				for (auto& p : pc) {
					p.x *= mul;
				}
				hull = qh.getConvexHull(pc,true,false);
				if (hull.getVertexBuffer().size() == 3) {
					break;
				}
			}
			
			// Test 3: 0D and 1D degenerate cases
			pc.clear();
			const Vector3<FloatType> a = Vector3<FloatType>(1,1,1).getNormalized();
			const Vector3<FloatType> b = Vector3<FloatType>(-2,4,9).getNormalized();
			const Vector3<FloatType> c(0,0,0);
			pc.push_back(a);
			hull = std::move(qh.getConvexHull(pc,true,false));
			assert(hull.getVertexBuffer().size() == 3);
			pc.push_back(b);
			hull = qh.getConvexHull(pc,true,false);
			assert(hull.getVertexBuffer().size() == 3);
			for (int i=0;i<N;i++) {
				auto t = rnd(0.001f,0.999f);
				auto d = a + (b-a)*t;
				pc.push_back(d);
			}
			hull = qh.getConvexHull(pc,true,false);
			assert(hull.getVertexBuffer().size() == 3);
			
			// Test 4: 2d degenerate cases
			pc.clear();
			pc.push_back(a);
			pc.push_back(b);
			pc.push_back(c);
			hull = qh.getConvexHull(pc,true,false);
			assert(hull.getVertexBuffer().size() == 3);
			for (int i=0;i<N;i++) {
				auto t = rnd(0.001f,0.999f);
				auto d = a*t;
				auto e = b*(1-t);
				auto f = c + d + e;
				pc.push_back(f);
			}
			hull = qh.getConvexHull(pc,true,false);
			assert(hull.getVertexBuffer().size() == 3);
			hull = qh.getConvexHull(pc,true,true);
			assert(hull.getVertexBuffer().size() == pc.size());
			assert(hull.getIndexBuffer().size() == 3);
			assert(&(hull.getVertexBuffer()[0])==&(pc[0]));
			
			// Other tests
			testPlanes();
			sphereTest();
			testVector3();
		}
		

		
	}
}
