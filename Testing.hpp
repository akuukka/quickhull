/*
 * Testing.hpp
 *
 *  Created on: Aug 16, 2014
 *      Author: anttiku
 */

#ifndef TESTING_HPP_
#define TESTING_HPP_

#include <chrono>

namespace quickhull {
	
	namespace testing {
	
		template<typename T>
		double getExecutionTime(T callable) {
			auto start = std::chrono::system_clock::now();
			callable();
			auto end = std::chrono::system_clock::now();
			std::chrono::duration<double> elapsed_seconds = end-start;
			return elapsed_seconds.count();
		}
		
	}

}



#endif /* TESTING_HPP_ */
