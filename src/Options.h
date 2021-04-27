/*
 * Options.h
 *
 *  Created on: Oct 14, 2020
 *      Author: cjustin
 */

#ifndef OPTIONS_H
#define OPTIONS_H 1

#include <stdint.h>
#include <string>

/**
 * Global variables that are mostly constant for the duration of the
 * execution of the program.
 */
namespace opt {

typedef uint16_t Count;

int verbose = 0;
unsigned k = 31;
int hashNum = 1;
uint64_t bfSize;
int threads = 1;
std::string bf = "";
}
#endif
