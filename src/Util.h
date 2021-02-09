/*
 * Util.h
 *
 * funct
 *
 *  Created on: Oct. 16, 2020
 *      Author: cjustin
 */

#ifndef SRC_UTIL_H_
#define SRC_UTIL_H_
#include <string>
#include <fstream>

using namespace std;

namespace Util {
/*
 * checks if file exists
 */
static bool fexists(const string &filename) {
	ifstream ifile(filename.c_str());
	bool good = ifile.good();
	ifile.close();
	return good;
}

/*
 * Upper triangular matrix to array conversion
 */
static unsigned matToIndex(unsigned i, unsigned j, unsigned n) {
	return ((n * (n - 1) / 2) - (n - i) * ((n - i) - 1) / 2 + j - i - 1);
}


}

#endif /* SRC_UTIL_H_ */
