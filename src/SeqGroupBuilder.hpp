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
#ifndef KSEQ_INIT_NEW
#define KSEQ_INIT_NEW
#include <zlib.h>
#include "vendor/kseq.h"
KSEQ_INIT(gzFile, gzread)
#endif /*KSEQ_INIT_NEW*/

using namespace std;

class SeqGroupBuilder {
public:

	typedef uint16_t GroupID;
//	typedef std::tuple<GroupID, uint16_t, uint16_t> Entry; //currently only store unique entries
	struct Entry {
		Entry() :
				m_uniqueCount(0), m_count(0) {

		}
//		Entry(uint16_t uniqueCount, uint16_t count) :
//				m_uniqueCount(uniqueCount), m_count(count) {
//		}
		uint16_t m_uniqueCount;
		uint16_t m_count;

	};
	typedef tsl::robin_map<uint64_t, shared_ptr<vector<Entry>>> GroupHash;

//	static const GroupID s_repeatID = std::numeric_limits<GroupID>::max();
//	static const GroupID s_bgID = std::numeric_limits<GroupID>::max() - 1;

	SeqGroupBuilder(const vector<string> &filenames) :
			m_filenames(filenames), m_sampleCounts(filenames.size(), 0) {
	}

	GroupHash loadFiles() {
		GroupHash tbl;
		//load in file 1 thread per file
		GroupID maxCount = m_filenames.size();

		//load in BF
		BloomFilter bf(opt::bf);
		assert(bf.getKmerSize() == opt::k);
//		assert(bf.getHashNum() == opt::hashNum);
		opt::hashNum = bf.getHashNum();

#pragma omp parallel for
		for (GroupID id = 0; id < maxCount; ++id) {
			gzFile fp;
			fp = gzopen(m_filenames[id].c_str(), "r");
			if (fp == Z_NULL) {
				std::cerr << "file " << m_filenames[id] << " cannot be opened"
						<< std::endl;
				exit(1);
			} else if (opt::verbose) {
				std::cerr << "Opening " << m_filenames[id] << std::endl;
			}
			kseq_t *seq = kseq_init(fp);
			int l = kseq_read(seq);
			while (l >= 0) {
#pragma omp atomic
				m_sampleCounts[id]++;
				l = kseq_read(seq);
			}
			gzclose(fp);

			//reopen
			fp = gzopen(m_filenames[id].c_str(), "r");
			seq = kseq_init(fp);
			l = kseq_read(seq);

			GroupHash uniqueEntries;
			vector<string> ids;

			unsigned index = 0;
			while (l >= 0) {
				ids.push_back(seq->name.s);
				//k-merize and insert
				for (ntHashIterator itr(seq->seq.s, opt::hashNum, opt::k);
						itr != itr.end(); ++itr) {
					//if kmer exists inside bg set
					if (!bf.contains(*itr)) {
						if (!tbl.contains((*itr)[0])) {
#pragma omp critical(tbl)
							tbl[(*itr)[0]] = shared_ptr<vector<Entry>>(
									new vector<Entry>(m_filenames.size()));
						}
						if (!uniqueEntries.contains((*itr)[0])) {
							uniqueEntries[(*itr)[0]] =
									shared_ptr<vector<Entry>>(
											new vector<Entry>(
													m_sampleCounts[id]));
						}
						if((*uniqueEntries[(*itr)[0]])[index].m_uniqueCount == 0){
							++(*uniqueEntries[(*itr)[0]])[index].m_uniqueCount;
#pragma omp atomic
							++(*tbl[(*itr)[0]])[id].m_uniqueCount;
						}
						++(*uniqueEntries[(*itr)[0]])[index].m_count;
#pragma omp atomic
						++(*tbl[(*itr)[0]])[id].m_count;
					}
				}
				l = kseq_read(seq);
				index++;
			}
			kseq_destroy(seq);
			gzclose(fp);

			if (m_sampleCounts.at(id) > 1) {
				unsigned matSize = (m_sampleCounts.at(id)
						* (m_sampleCounts.at(id) - 1)) / 2;

				//create similarity matrix (jaccard)
				int *countMat = new int[matSize];
				//init array
				for (unsigned i = 0; i < matSize; ++i) {
					countMat[i] = 0;
				}

				//populate count matrix
				vector<uint64_t> groupCounts(m_sampleCounts.at(id), 0);

				for (SeqGroupBuilder::GroupHash::iterator itr =
						uniqueEntries.begin(); itr != uniqueEntries.end();
						++itr) {
					uint16_t lastCount = itr->second->at(0).m_count;
					bool allSame = true;
					for (unsigned i = 1; i < m_sampleCounts.at(id); ++i) {
						if(itr->second->at(i).m_count != lastCount){
							allSame = false;
							break;
						}
					}
					if (!allSame) {
						for (unsigned i = 0; i < m_sampleCounts.at(id); ++i) {
							if (itr->second->at(i).m_count != 0) {
								for (unsigned j = i + 1;
										j < m_sampleCounts.at(id); ++j) {
									if (itr->second->at(j).m_count
											== itr->second->at(i).m_count) {
										++countMat[Util::matToIndex(i, j,
												m_sampleCounts.at(id))];
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

				for (unsigned i = 0; i < m_sampleCounts.at(id); ++i) {
					for (unsigned j = i + 1; j < m_sampleCounts.at(id); ++j) {

						unsigned setSize = groupCounts.at(i) + groupCounts.at(j)
								- countMat[Util::matToIndex(i, j,
										m_sampleCounts.at(id))];
						distMat[Util::matToIndex(i, j, m_sampleCounts.at(id))] =
								double(
										countMat[Util::matToIndex(i, j,
												m_sampleCounts.at(id))])
										/ double(setSize);

//						cout << Util::matToIndex(i, j, m_sampleCounts.at(id))
//								<< "\t" << i << "\t" << j << "\t"
//								<< double(
//										countMat[Util::matToIndex(i, j,
//												m_sampleCounts.at(id))])
//										/ double(setSize) << "\t" << setSize
//								<< endl;
					}
				}

				//output matrix
				printMatrix(distMat,ids, m_filenames[id] + ".mat.tsv");

				//output sample to matrix table
				printTable(uniqueEntries,ids, m_filenames[id] + ".count.tsv");
			}

		}
		return tbl;
	}

	const vector<uint16_t>& getSampleCount() {
		return m_sampleCounts;
	}

	virtual ~SeqGroupBuilder() {
//		delete m_countMat;
	}

private:
	const vector<string> &m_filenames;
	vector<uint16_t> m_sampleCounts;

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
	void printTable(const SeqGroupBuilder::GroupHash &uniqueEntries,
			const vector<string> &ids, string outputFilename) const {
		ofstream out(outputFilename.c_str());
		out << "sample";
		size_t count = 0;
		for (SeqGroupBuilder::GroupHash::const_iterator itr =
				uniqueEntries.begin(); itr != uniqueEntries.end(); ++itr) {
			opt::Count lastCount = itr->second->at(0).m_count;
			bool allSame = true;
			for (unsigned i = 1; i < itr->second->size(); ++i) {
				if (itr->second->at(i).m_count != lastCount) {
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
			for (SeqGroupBuilder::GroupHash::const_iterator itr =
					uniqueEntries.begin(); itr != uniqueEntries.end(); ++itr) {
				opt::Count lastCount = itr->second->at(0).m_count;
				bool allSame = true;
				for (unsigned i = 1; i < itr->second->size(); ++i) {
					if (itr->second->at(i).m_count != lastCount) {
						allSame = false;
						break;
					}
				}
				if (!allSame) {
					out << "," << itr->second->at(i).m_count;
				}
			}
			out << "\n";
		}
	}



};

#endif /* SRC_SEQGROUPBUILDER_HPP_ */

