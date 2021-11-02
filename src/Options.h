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

using namespace std;

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
string bf = "";
string outputPrefix = "";
string groupingsFile = "";
string readInput = "";
double pseudoCount = 0.001; //for dealing with kl distance zero values
bool haploid = false;
bool debug = false;
//bool checkDupes = true;
}
#endif
