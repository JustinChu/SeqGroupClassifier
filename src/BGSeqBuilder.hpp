/*
 * BGSeqBuilder.hpp
 *
 * Intend to create an object for storing background sequences
 *
 *  Created on: Jan 12, 2021
 *      Author: cjustin
 */

#ifndef SRC_BGSEQBUILDER_HPP_
#define SRC_BGSEQBUILDER_HPP_

#include "vendor/btl_bloomfilter/BloomFilter.hpp"
#include "vendor/ntHash/ntHashIterator.hpp"
#include "Options.h"
#include <vector>
#include <string>
#include <zlib.h>
#include <omp.h>
#include <stdio.h>
#include "Options.h"
#include "vendor/ntHash/ntHashIterator.hpp"
#ifndef KSEQ_INIT_NEW
#define KSEQ_INIT_NEW
#include "vendor/kseq.h"
KSEQ_INIT(gzFile, gzread)
#endif /*KSEQ_INIT_NEW*/

class BGSeqBuilder {
public:
	BGSeqBuilder::BGSeqBuilder() {
		// TODO Auto-generated constructor stub

	}

	BloomFilter storeBG(const vector<string> &filenames) {
		BloomFilter bf(opt::bfSize,opt::hashNum,opt::k);
#pragma omp parallel for
		for (unsigned i = 0; i < filenames.size(); ++i) {
			gzFile fp;
			fp = gzopen(filenames[i].c_str(), "r");
			if (fp == Z_NULL) {
				std::cerr << "file " << filenames[i] << " cannot be opened"
						<< std::endl;
				exit(1);
			}
			kseq_t *seq = kseq_init(fp);
			int l;
			while (l >= 0) {
				l = kseq_read(seq);
				if (l >= 0) {
					//k-merize and insert
					for (ntHashIterator itr(seq->seq.s, opt::hashNum, opt::k);
							itr != itr.end(); ++itr) {
						bf.insert(*itr);
					}
				} else {
					break;
				}
			}
			kseq_destroy(seq);
			gzclose(fp);
		}
		return tbl;
	}
	//load in file
	//populate bloom filter

	BGSeqBuilder::~BGSeqBuilder() {
		// TODO Auto-generated destructor stub
	}


};
