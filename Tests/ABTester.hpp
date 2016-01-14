//
//  ABTesting.hpp
//  QuickHull
//
//  Created by Antti Kuukka on 1/7/16.
//  Copyright © 2016 Antti Kuukka. All rights reserved.
//

#ifndef ABTesting_h
#define ABTesting_h

#include <iostream>
#include <chrono>

namespace quickhull {

	class ABTester {
		public:
		static bool a;
		
		size_t n[2];
		double t[2];
		
		double average(bool A) {
			const size_t ind = A?0:1;
			return t[ind]/n[ind];
		}
		
		void addResult(bool A, double time) {
			const size_t ind = A?0:1;
			n[ind]++;
			t[ind]+=time;
		}
		
		void printResults()
		{
			for (size_t i=0;i<2;i++)
			{
				std::cout << (i==0?"A":"B");
				std::cout << ": ";
				std::cout << "N=" << n[i] << " Average time: " << average(i==0);
				std::cout << std::endl;
			}
		}
		
		template <typename Callable>
		void run(size_t N, Callable C) {
			for (size_t i=0;i<N*2;i++) {
				a = (i%2==0);
				auto start = std::chrono::system_clock::now();
				C();
				auto end = std::chrono::system_clock::now();
				std::chrono::duration<double> elapsed_seconds = end-start;
				addResult(a,elapsed_seconds.count());
			}
			printResults();
		}
		
		ABTester() {
			clear();
		}
		
		void clear() {
			n[0]=0;
			n[1]=0;
			t[0]=0;
			t[1]=0;
		}
	};
	
}


#endif /* ABTesting_h */
