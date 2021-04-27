/*
 * SeqGroupClassifier.hpp
 *
 *  Created on: Apr 21, 2021
 *      Author: cjustin
 */

//General plan
//load in group file
//load in k-mer profiles tsl::robin_map<uint64_t, shared_ptr<vector<Entry>>> GroupHash
//create profile
//vector of doubles with k-mer frequency/ number of samples
//for each entry in group sum frequency/
//1) compute total number of k-mers
//2) compute frequencyMatrix
//3) for each k-mer take sum of counts and divide by total number
//4) load into matrix tsl::robin_map<uint64_t, shared_ptr<vector<double>>>
//5) compute running KL distance (w/ Kahan summation algorithm?)
#ifndef SRC_SEQGROUPCLASSIFIER_HPP_
#define SRC_SEQGROUPCLASSIFIER_HPP_

#include <vector>
#include <string>
#include <omp.h>
#include <stdio.h>
#include <math.h>
#include "vendor/tsl/robin_map.h"
#include "Options.h"

//#include "vendor/tsl/robin_set.h"

using namespace std;

class SeqGroupClassifier {
public:

	typedef uint16_t GroupID; //index of vector described by m_groupIDs
	typedef uint16_t SampleID; //index of vector described by m_groupIDs
	typedef tsl::robin_map<uint64_t, uint32_t> indexHash; //indexes the matrix used
	typedef tsl::robin_map<uint64_t, shared_ptr<vector<opt::Count>>> CountHash; //input vector converted to matrix

	SeqGroupClassifier(const CountHash &counts, const vector<string> &sampleIDs,
			const vector<string> &groupIDs,
			const vector<vector<GroupID>> &groupings) :
			m_groupIDs(groupIDs), m_sampleIDs(sampleIDs), m_groupings(
					groupings), m_hashToIndex() {
		computeGroupFreq(counts);
	}

	/*
	 * Compute frequency matrix for each of group ID
	 * Creates a hash value (derived from k-mer) to matrix
	 */
	void computeGroupFreq(const CountHash &counts) {
		//figure out counts of each grouping to calculate frequency
		vector<uint64_t> totalCounts(m_groupings.size(), 0);
		uint64_t freqMatSize = 0;
		//for each k-mer
		for (CountHash::const_iterator i = counts.begin(); i != counts.end();
				i++) {
			//determine if common to all samples
			uint16_t lastCount = i->second->at(0);
			bool allSame = true;
			for (unsigned j = 1; j < m_sampleIDs.size(); ++j) {
				if(i->second->at(j) != lastCount){
					allSame = false;
					break;
				}
			}
			if (!allSame) {
				//if not common to all samples count total number to group and increment size
				//add to index
				m_hashToIndex[i->first] = freqMatSize++;
				//assign count of all k-mers to group
				for (size_t j = 0; j != m_groupings.size(); ++j) {
					for (vector<GroupID>::const_iterator k = m_groupings[j].begin();
							k != m_groupings[j].end(); ++j) {
						totalCounts[j] += i->second->at(*k);
					}
				}
			}
		}
		//set size of vectors
		for (unsigned i = 0; i < m_sampleIDs.size(); ++i) {
			m_groupFreq.push_back(
					shared_ptr<vector<opt::Count>>(
							new vector<opt::Count>(freqMatSize, 0)));
		}

		//populate matrix
		for (CountHash::const_iterator i = counts.begin(); i != counts.end();
				i++) {
			indexHash::iterator itr = m_hashToIndex.find(i->first);
			//k-mer is not the sample in all samples
			if (itr != m_hashToIndex.end()) {
				for (size_t j = 0; j != m_groupings.size(); ++j) {
					uint64_t countOfGroup = 0;
					for (vector<GroupID>::const_iterator k =
							m_groupings[j].begin(); k != m_groupings[j].end();
							++j) {
						countOfGroup += i->second->at(*k);
					}
					double freq = double(countOfGroup) / double(totalCounts[j]);
					(*m_groupFreq[j])[j] = freq;
				}
			}
		}
	}

//	/*
//	 * Creates combination of each group frequencies
//	 * Intended to blend multiple profiles together for diploid genomes
//	 */
//	void createCombFreq() {
//
//	}

//	/*
//	 * Combines the frequency of 2 vectors of counts
//	 */
//	vector<double> combineFrequency(const vector<double> &mat1,
//			const vector<double> &mat2) {
//
//	}

	/*
	 * Assuming diploid genome compute KL distances metric
	 */
	double computeKLDist(const vector<double> &parent1,
			const vector<double> &parent2,
			const vector<double> &sampleFreq) {
		double dist = 0.0;
		for (size_t i = 0; i < sampleFreq.size(); ++i) {
			double blendedFreq = (parent1[i] + parent2[i]) / 2.0;
			//using base 2 log, so units use are "bits" rather than "nats" (base e)
			double subDist = sampleFreq[i] * log2(sampleFreq[i] / blendedFreq);
			dist += subDist;
		}
		return(dist);
	}

//	double computeML(){
//
//	}

//	double computeBaysianEstimate(){
//
//	}

	virtual ~SeqGroupClassifier();
private:
	const vector<string> &m_groupIDs;
	const vector<string> &m_sampleIDs;
	const vector<vector<GroupID>> &m_groupings;
//	tsl::robin_map<SampleID, GroupID> &m_sampleToGroup;
	indexHash &m_hashToIndex; //hashed k-mer value to index
	vector<shared_ptr<vector<double>>> m_groupFreq;
};

#endif /* SRC_SEQGROUPCLASSIFIER_HPP_ */
