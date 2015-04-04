/*
 * helper.hpp
 *
 *  Created on: Apr 4, 2015
 *      Author: karthik
 */

#ifndef HELPER_HPP_
#define HELPER_HPP_

template <typename T>
class compareFunction{
private:
	int scaling;
public:
	compareFunction(int scale):scaling(scale){}
	bool operator()(const T & v1, const T& v2){
		return (int)(v1*scaling) < (int) (v2*scaling);
	}
};


#endif /* HELPER_HPP_ */
