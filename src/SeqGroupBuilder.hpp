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

#include <vector>
#include <string>
#include <omp.h>
#include <stdio.h>
#include <limits>
#include <memory>
//#include <tuple>
#include "Options.h"
#include "vendor/btl_bloomfilter/BloomFilter.hpp"
#include "vendor/ntHash/ntHashIterator.hpp"
#include "vendor/tsl/robin_map.h"
#include "vendor/tsl/robin_set.h"
#include <zlib.h>
#include "vendor/kseq.h"
#ifndef KSEQ_INIT_NEW
#define KSEQ_INIT_NEW
KSEQ_INIT(gzFile, gzread)
#endif /*KSEQ_INIT_NEW*/

using namespace std;

class SeqGroupBuilder {
public:
	typedef tsl::robin_map<uint64_t, shared_ptr<vector<opt::Count>>> CountsHash;

//	static const GroupID s_repeatID = std::numeric_limits<GroupID>::max();
//	static const GroupID s_bgID = std::numeric_limits<GroupID>::max() - 1;

	SeqGroupBuilder(const vector<string> &filenames) :
			m_filenames(filenames) {
	}

	CountsHash loadFiles() {
		CountsHash counts;
		//load in BF
		BloomFilter bf(opt::bf);
		assert(bf.getKmerSize() == opt::k);
//		assert(bf.getHashNum() == opt::hashNum);
		opt::hashNum = bf.getHashNum();

		unsigned sampleCounts = 0;

#pragma omp parallel for
		for (unsigned i = 0; i < m_filenames.size(); ++i) {
			gzFile fp;
			fp = gzopen(m_filenames[i].c_str(), "r");
			if (fp == Z_NULL) {
				std::cerr << "file " << m_filenames[i] << " cannot be opened"
						<< std::endl;
				exit(1);
			} else if (opt::verbose) {
				std::cerr << "Opening " << m_filenames[i] << std::endl;
			}
			kseq_t *seq = kseq_init(fp);
			int l = kseq_read(seq);
			while (l >= 0) {
#pragma omp atomic
				sampleCounts++;
				l = kseq_read(seq);
			}
			gzclose(fp);
		}

		m_sampleIDs.reserve(sampleCounts);

#pragma omp parallel for
		for (unsigned i = 0; i < m_filenames.size(); ++i) {
			gzFile fp;
			fp = gzopen(m_filenames[i].c_str(), "r");
			if (fp == Z_NULL) {
				std::cerr << "file " << m_filenames[i] << " cannot be opened"
						<< std::endl;
				exit(1);
			} else if (opt::verbose) {
				std::cerr << "Opening " << m_filenames[i] << std::endl;
			}
			kseq_t *seq = kseq_init(fp);
			int l = kseq_read(seq);
			unsigned index = 0;
			while (l >= 0) {
				m_sampleIDs.push_back(seq->name.s);
				//k-merize and insert
				for (ntHashIterator itr(seq->seq.s, opt::hashNum, opt::k);
						itr != itr.end(); ++itr) {
					//if kmer exists inside bg set
					if (!bf.contains(*itr)) {
						if (!counts.contains((*itr)[0])) {
#pragma omp critical(counts)
							counts[(*itr)[0]] = shared_ptr<vector<opt::Count>>(
									new vector<opt::Count>(sampleCounts));
						}
#pragma omp atomic
						++(*counts[(*itr)[0]])[index];
					}
				}
				l = kseq_read(seq);
				index++;
			}
			kseq_destroy(seq);
			gzclose(fp);
		}
		return counts;
	}

	void printMatrixes(const CountsHash &counts) {
		unsigned matSize = (m_sampleIDs.size() * (m_sampleIDs.size() - 1)) / 2;

		//create similarity matrix (jaccard)
		int *countMat = new int[matSize];
		//init array
		for (unsigned i = 0; i < matSize; ++i) {

			countMat[i] = 0;
		}

		//populate count matrix
		vector<uint64_t> groupCounts(m_sampleIDs.size(), 0);

		for (CountsHash::const_iterator itr = counts.begin(); itr != counts.end();
				++itr) {
			uint16_t lastCount = itr->second->at(0);
			bool allSame = true;
			for (unsigned i = 1; i < m_sampleIDs.size(); ++i) {
				if (itr->second->at(i) != lastCount) {
					allSame = false;
					break;
				}
			}
			if (!allSame) {
				for (unsigned i = 0; i < m_sampleIDs.size(); ++i) {
					if (itr->second->at(i) != 0) {
						for (unsigned j = i + 1; j < m_sampleIDs.size(); ++j) {
							if (itr->second->at(j) == itr->second->at(i)) {
								++countMat[Util::matToIndex(i, j, m_sampleIDs.size())];
							}
						}
						++groupCounts[i];
					}
				}
			}
		}

		//create similarity matrix (jaccard)
		double *distMat = new double[matSize];
		//init array
		for (unsigned i = 0; i < matSize; ++i) {
			distMat[i] = 0;
		}

		for (unsigned i = 0; i < m_sampleIDs.size(); ++i) {
			for (unsigned j = i + 1; j < m_sampleIDs.size(); ++j) {

				unsigned setSize =
						groupCounts.at(i) + groupCounts.at(j)
								- countMat[Util::matToIndex(i, j,
										m_sampleIDs.size())];
				distMat[Util::matToIndex(i, j, m_sampleIDs.size())] = double(
						countMat[Util::matToIndex(i, j, m_sampleIDs.size())])
						/ double(setSize);

			}

			//output matrix
			printMatrix(distMat, m_sampleIDs, opt::outputPrefix + ".mat.tsv");

			//output sample to matrix table
			printTable(counts, m_sampleIDs, opt::outputPrefix + ".count.tsv");
		}
	}

	const vector<string> &getSampleIDs() const{
		return m_sampleIDs;
	}

	virtual ~SeqGroupBuilder() {
//		delete m_countMat;
	}

private:
	const vector<string> &m_filenames;
	vector<string> m_sampleIDs;

	/*
	 * prints to file a tab separated distance matrix
	 */
	void printMatrix(const double *distMat, const vector<string> &ids,
			string outputFilename) const {
		ofstream out(outputFilename.c_str());
		out << "sample";
		//print out headerIDs
		for (vector<string>::const_iterator itr = ids.begin(); itr != ids.end();
				++itr) {
			out << "," << *itr;
		}
		out << "\n";

//		unsigned matSize = (ids.size() * (ids.size() - 1)) / 2;

		for (unsigned i = 0; i < ids.size(); ++i) {
			out << ids.at(i);
			for (unsigned j = 0; j < ids.size(); ++j) {
				if (i == j) {
					out << ",1";
				} else if (i > j) {
					out << "," << distMat[Util::matToIndex(j, i, ids.size())];
				} else {
					out << "," << distMat[Util::matToIndex(i, j, ids.size())];
				}
			}
			out << "\n";
		}
	}

	/*
	 * prints to file a tab separated matrix of k-mer counts
	 */
	void printTable(const CountsHash &uniqueEntries,
			const vector<string> &ids, string outputFilename) const {
		ofstream out(outputFilename.c_str());
		out << "sample";
		size_t count = 0;
		for (CountsHash::const_iterator itr =
				uniqueEntries.begin(); itr != uniqueEntries.end(); ++itr) {
			opt::Count lastCount = itr->second->at(0);
			bool allSame = true;
			for (unsigned i = 1; i < itr->second->size(); ++i) {
				if (itr->second->at(i) != lastCount) {
					allSame = false;
					break;
				}
			}
			if (!allSame) {
				out << "," << count++;
			}
		}
		out << "\n";

		for (unsigned i = 0; i < ids.size(); ++i) {
			out << ids.at(i);
			for (CountsHash::const_iterator itr =
					uniqueEntries.begin(); itr != uniqueEntries.end(); ++itr) {
				opt::Count lastCount = itr->second->at(0);
				bool allSame = true;
				for (unsigned i = 1; i < itr->second->size(); ++i) {
					if (itr->second->at(i) != lastCount) {
						allSame = false;
						break;
					}
				}
				if (!allSame) {
					out << "," << itr->second->at(i);
				}
			}
			out << "\n";
		}
	}

}
;

#endif /* SRC_SEQGROUPBUILDER_HPP_ */

