//
//  Pool.hpp
//  QuickHull
//
//  Created by Antti Kuukka on 1/28/16.
//  Copyright © 2016 Antti Kuukka. All rights reserved.
//

#ifndef Pool_h
#define Pool_h

#include <vector>
#include <memory>

namespace quickhull {
	
	template<typename T>
	class Pool {
		std::vector<std::unique_ptr<T>> m_data;
	public:
		
	};
	
}

#endif /* Pool_h */
