/*
 * SeqGroupBuilder.hpp
 *
 * 	Loads in fast files and breaks them into subsequences that are hashed
 * 	These hased values are then inserted to hash table (may replace datastructure in the future)
 *
 * 	Non-mutable hashtable can then be return for use in classifier
 *
 *  Created on: Jan 4, 2021
 *      Author: cjustin
 */

#ifndef SRC_SEQGROUPBUILDER_HPP_
#define SRC_SEQGROUPBUILDER_HPP_

#include "vendor/tsl/robin_map.h"
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

using namespace std;

class SeqGroupBuilder {
public:

	typedef uint16_t groupID;

	SeqGroupBuilder() {
	}

	tsl::robin_map<uint64_t, groupID> loadFiles(
			const vector<string> &filenames) {
		tsl::robin_map<uint64_t, groupID> tbl;
		//load in file 1 thread per file

#pragma omp parallel for
		for (groupID i = 0; i < filenames.size(); ++i) {
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
#pragma omp critical(tbl)
						tbl[(*itr)[0]] = i;
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

//	const vector<string>& getGroupList() {
//		return m_groupList;
//	}

	virtual ~SeqGroupBuilder() {
	}

private:
//	vector<string> m_groupList;

};

#endif /* SRC_SEQGROUPBUILDER_HPP_ */

