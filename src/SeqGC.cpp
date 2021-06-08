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
#include <limits.h>
#include "config.h"
#include "src/Options.h"
#include "src/Util.h"
#include "SeqGroupBuilder.hpp"
#include "SeqGroupClassifier.hpp"
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
//	"Usage: " PROGRAM " build [OPTION]... [FASTA]...\n"
	"Usage: " PROGRAM " [OPTION]... [FASTA]...\n"
	"  -p, --prefix           Output name prefix. [required]\n"
	"  -g, --groups           File outlining groupings.\n"
	"  -l, --haploid          Run in haploid mode.\n"
	"  -i, --input            Input sequences to classify.\n"
	"  -c, --pseudo_count     KL distance value for smoothing.[0.001]\n"
//	"  -t, --threads          Number of threads to run.[1]\n"
	"  -b, --bg_bloom         Bloom filter of background sequences. [required]\n"
	"  -d, --debug            Debug file is produced.\n"
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
		"prefix", required_argument, NULL, 'p' }, {
		"groups", required_argument, NULL, 'g' }, {
		"haploid", no_argument, NULL, 'l' }, {
		"input", required_argument, NULL, 'i' }, {
		"pseudo_count", required_argument, NULL, 'c' }, {
		"threads", required_argument, NULL, 't' }, {
		"bg_bloom", required_argument, NULL, 'b' }, {
		"debug", required_argument, NULL, 'd' }, {
		"help", no_argument, NULL, 'h' }, {
		"version", no_argument, &OPT_VERSION, 1 }, {
		"verbose", no_argument, NULL, 'v' }, {
		NULL, 0, NULL, 0 } };

	int option_index = 0;
	while ((c = getopt_long(argc, argv, "p:g:i:lc:t:vdhb:", long_options,
			&option_index)) != -1)
	{
		istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
		case 'p': {
			stringstream convert(optarg);
			if (!(convert >> opt::outputPrefix)) {
				cerr << "Error - Invalid parameter p: "
						<< optarg << endl;
				return 0;
			}
			break;
		}
		case 'g': {
			stringstream convert(optarg);
			if (!(convert >> opt::groupingsFile)) {
				cerr << "Error - Invalid parameter g: "
						<< optarg << endl;
				return 0;
			}
			break;
		}
		case 'i': {
			stringstream convert(optarg);
			if (!(convert >> opt::readInput)) {
				cerr << "Error - Invalid parameter i: "
						<< optarg << endl;
				return 0;
			}
			break;
		}
		case 'l': {
			opt::haploid = true;
			break;
		}
		case 'c': {
			stringstream convert(optarg);
			if (!(convert >> opt::pseudoCount)) {
				cerr << "Error - Invalid parameter c: "
						<< optarg << endl;
				return 0;
			}
			break;
		}
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
		case 'd': {
			opt::debug = true;
			break;
		}
		case 'h': {
//			printHelpDialog();
			printBuildDialog();
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
		assert(Util::fexists(inputFiles.back()));
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
		cerr << "Error: Need Input background sequence Bloom filter (-b)" << endl;
		die = true;
	}

	//check for background bf
	if(opt::outputPrefix.empty()){
		cerr << "Error: Need output prefex (-p)" << endl;
		die = true;
	}

	if (die) {
		cerr << "Try '--help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	//load sequence
	SeqGroupBuilder builder(inputFiles);
	SeqGroupBuilder::CountsHash counts = builder.loadFiles();

	if(opt::groupingsFile.empty()){
		builder.printMatrixes(counts);
	}
	else{
		SeqGroupClassifier classifier(counts, builder.getSampleIDs());
		if(!opt::readInput.empty()){
			if (Util::fexists(opt::readInput)) {
				if (opt::haploid) {
//					tsl::robin_map<SeqGroupClassifier::SampleID, double> results =
//							classifier.computeAllKLDist(opt::readInput);
//					classifier.printResults(results);
				} else {
//					classifier.computeDiploidKLDist(opt::readInput);
//					classifier.printResults(results);
					classifier.computeDiploid(opt::readInput);
				}
			}
		}
	}

	return 0;
}




