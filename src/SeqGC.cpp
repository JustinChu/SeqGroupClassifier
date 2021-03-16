/*
 * SeqGC.cpp
 *
 *  Created on: Dec 14, 2020
 *      Author: cjustin
 */

/*
 * PanGenomeRefBuild.cpp
 *
 *  Created on: Oct. 13, 2020
 *      Author: cjustin
 */

#include <sstream>
#include <string>
#include <vector>
#include <getopt.h>
#include <iostream>
#include <stdlib.h>
#include "config.h"
#include "src/Options.h"
#include "src/Util.h"
#include "SeqGroupBuilder.hpp"
#if _OPENMP
# include <omp.h>
#endif

using namespace std;

#define PROGRAM "seqgc"

void printVersion()
{
	const char VERSION_MESSAGE[] = PROGRAM " (" PACKAGE_NAME ") " GIT_REVISION "\n"
	"Written by Justin Chu <cjustin@ds.dfci.harvard.edu>\n"
	"\n"
	"Copyright 2020 Dana-Farber Cancer Institute\n";
	cerr << VERSION_MESSAGE << endl;
	exit(EXIT_SUCCESS);
}

void printHelpDialog()
{
	const char dialog[] =
	"Usage: " PROGRAM " [FUNCTION]\n"
	"Functions:\n"
	"  build\n";

	cerr << dialog << endl;
	exit(EXIT_SUCCESS);
}

void printBuildDialog(){
	const char dialog[] =
	"Usage: " PROGRAM " build [OPTION]... [FILES]...\n"
	"The input is expected to be a set of FASTA files\n\n"
	"  -t, --threads          Number of threads to run.[1]\n"
	"  -b, --bg_bloom         Bloom filter of background sequences.[required]"
	"  -h, --help             Display this dialog.\n"
	"  -v, --verbose          Display verbose output.\n"
	"      --version          Print version information.\n";

	cerr << dialog << endl;
	exit(EXIT_SUCCESS);
}

int main(int argc, char *argv[])
{
	//switch statement variable
	int c;

	//control variables
	bool die = false;
	int OPT_VERSION = 0;

	//long form arguments
	static struct option long_options[] = { {
		"threads", required_argument, NULL, 't' }, {
		"bg_bloom", required_argument, NULL, 'b' }, {
		"help", no_argument, NULL, 'h' }, {
		"version", no_argument, &OPT_VERSION, 1 }, {
		"verbose", no_argument, NULL, 'v' }, {
		NULL, 0, NULL, 0 } };

	int option_index = 0;
	while ((c = getopt_long(argc, argv, "t:vhb:", long_options,
			&option_index)) != -1)
	{
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
		case 't': {
			stringstream convert(optarg);
			if (!(convert >> opt::threads)) {
				cerr << "Error - Invalid parameter t: "
						<< optarg << endl;
				return 0;
			}
			break;
		}
		case 'b': {
			stringstream convert(optarg);
			if (!(convert >> opt::bf)) {
				cerr << "Error - Invalid parameter b: "
						<< optarg << endl;
				return 0;
			}
			break;
		}
		case 'h': {
			printHelpDialog();
			break;
		}
		case 'v': {
			opt::verbose++;
			break;
		}
		case '?': {
			die = true;
			break;
		}
		}
	}

#if defined(_OPENMP)
	if (opt::threads > 0)
	omp_set_num_threads(opt::threads);
#endif

	if (OPT_VERSION) {
		printVersion();
	}

	vector<string> inputFiles;
	while (optind < argc) {
		inputFiles.emplace_back(argv[optind]);
		Util::fexists(inputFiles.back());
		optind++;
		//check if file exists
	}
	//Check needed options
	if (inputFiles.size() == 0) {
		cerr << "Error: Need Input File" << endl;
		die = true;
	}

	//check for background bf
	if(!Util::fexists(opt::bf)){
		cerr << "Error: Need Input background sequence Bloom filter" << endl;
		die = true;
	}

	if (die) {
		cerr << "Try '--help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	//load sequence
	SeqGroupBuilder builder(inputFiles);
	SeqGroupBuilder::GroupHash groups = builder.loadFiles();
	const vector<uint16_t>& sampleCounts = builder.getSampleCount();
	vector<uint64_t> groupCounts(inputFiles.size(),0);
	vector<uint64_t> groupCountInAll(inputFiles.size(),0);
	vector<uint64_t> maxCount(inputFiles.size(),0);

	//mat length
	unsigned matSize = (inputFiles.size() * (inputFiles.size() - 1)) / 2;

	//create similarity matrix
	int *countMat = new int[matSize];
	//init array
	for(unsigned i = 0; i < matSize; ++i){
		countMat[i] = 0;
	}

//	unsigned commonKmers = 0;
	for (SeqGroupBuilder::GroupHash::iterator itr = groups.begin(); itr != groups.end();
			++itr) {
		for (unsigned i = 0; i < inputFiles.size(); ++i) {
			if (itr->second->at(i).m_count != 0) {
				if (opt::verbose > 1) {
					cout << inputFiles[i] << "\t"
							<< itr->second->at(i).m_uniqueCount << "\t"
							<< itr->second->at(i).m_count << "\t"
							<< sampleCounts.at(i) << endl;
				}
				if (sampleCounts.at(i) == itr->second->at(i).m_uniqueCount) {
					++groupCountInAll[i];
				}
				if (maxCount[i] < itr->second->at(i).m_uniqueCount) {
					maxCount[i] = itr->second->at(i).m_uniqueCount;
				}
				for (unsigned j = i + 1; j < inputFiles.size(); ++j) {
					if (itr->second->at(j).m_count != 0) {
						++countMat[Util::matToIndex(i, j, inputFiles.size())];
						if (opt::verbose > 1) {
							cout << Util::matToIndex(i, j, inputFiles.size()) << "\t"
									<< inputFiles[i] << "\t" << inputFiles[j]
									<< "\t" << itr->second->at(i).m_uniqueCount
									<< "\t" << itr->second->at(j).m_uniqueCount
									<< "\t"
									<< countMat[Util::matToIndex(i, j,
											itr->second->size())] << endl;
						}
					}
				}
				++groupCounts[i];
			}
		}
	}
	for (SeqGroupBuilder::GroupID i = 0; i != inputFiles.size(); ++i) {
		cout << inputFiles[i] << "\t" << groupCounts[i] << "\t"
				<< groupCountInAll[i] << "\t" << maxCount[i]
				<< endl;
	}
//	cout << "Common K-mers: " <<  commonKmers << endl;
	for (unsigned i = 0; i < inputFiles.size(); ++i) {
		for (unsigned j = i + 1; j < inputFiles.size(); ++j) {
			cout << Util::matToIndex(i, j, inputFiles.size()) << "\t" << inputFiles[i]
					<< "\t" << inputFiles[j] << "\t"
					<< double(countMat[Util::matToIndex(i, j, inputFiles.size())])
							/ double(
									groupCounts.at(i) + groupCounts.at(j)
											- countMat[Util::matToIndex(i, j,
													inputFiles.size())]) << "\t"
					<< groupCounts.at(i) + groupCounts.at(j)
							- countMat[Util::matToIndex(i, j, inputFiles.size())]
					<< endl;
		}
	}
	//serialize datastructures

	return 0;
}




