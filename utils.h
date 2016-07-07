//
//  utils.h
//  naf
//
//  Created by David CHAMPREDON on 2016-07-04.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#ifndef naf_utils_h
#define naf_utils_h

#include <iostream>
#include <vector>


using namespace std;

// Customized types

using ID = unsigned int;


// Global constants

const ID __UNDEFINED_ID = 999999999;

# define TINY		1e-7
# define SUPERTINY	1e-10


// Useful functions


void stopif(bool condition, string error_msg,
			int error_code=1, const char ff[]=__FUNCTION__);


template <class T> T sumElements(vector<T> x)
{
	/// Returns the sum of all elements from a vector
	
	T s = 0;
	for (int i=0; i<x.size(); i++)	{
		s += x[i];
	}
	return s;
}

#endif
