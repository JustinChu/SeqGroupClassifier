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
#include <cstring>

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

/*
 * Get resident set size of memory
 */
static size_t getRSS(){ //Note: this value is in KB!
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];

    while (fgets(line, 128, file) != NULL){
        if (strncmp(line, "VmRSS:", 6) == 0){
            int i = strlen(line);
            const char* p = line;
            while (*p <'0' || *p > '9') p++;
            line[i-3] = '\0';
            result = atoi(p);
            break;
        }
    }
    fclose(file);
    return result;
}

///*
// * generic hash function
// */
//static long murmur64(long h) {
//	h ^= h >> 33;
//	h *= 0xff51afd7ed558ccdL;
//	h ^= h >> 33;
//	h *= 0xc4ceb9fe1a85ec53L;
//	h ^= h >> 33;
//	return h;
//}

}

#endif /* SRC_UTIL_H_ */
