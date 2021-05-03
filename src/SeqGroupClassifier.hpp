/*
 * SeqGroupClassifier.hpp
 *
 *  Created on: Apr 21, 2021
 *      Author: cjustin
 */
#ifndef SRC_SEQGROUPCLASSIFIER_HPP_
#define SRC_SEQGROUPCLASSIFIER_HPP_

#include <vector>
#include <string>
#include <omp.h>
#include <stdio.h>
#include <math.h>
#include <fstream>

#include "Options.h"

#include "vendor/tsl/robin_map.h"
#include "vendor/btl_bloomfilter/BloomFilter.hpp"
#include "vendor/ntHash/ntHashIterator.hpp"
#include "vendor/kseq.h"

#ifndef KSEQ_INIT_NEW
#define KSEQ_INIT_NEW
#include <zlib.h>
#include "vendor/kseq.h"
KSEQ_INIT(gzFile, gzread)
#endif /*KSEQ_INIT_NEW*/

using namespace std;

class SeqGroupClassifier {
public:

	// A hash function used to hash a pair of any kind
	struct HashPair {
	    template <class T1, class T2>
	    size_t operator()(const pair<T1, T2>& p) const
	    {
	        auto hash1 = hash<T1>{}(p.first);
	        auto hash2 = hash<T2>{}(p.second);
	        return hash1 ^ hash2;
	    }
	};

	typedef uint16_t GroupID; //index of vector described by m_groupIDs
	typedef uint16_t SampleID; //index of vector described by m_groupIDs
	typedef uint32_t IndexPos;
	typedef tsl::robin_map<uint64_t, IndexPos> indexHash; //indexes the matrix used
	typedef tsl::robin_map<uint64_t, shared_ptr<vector<opt::Count>>> CountHash; //input vector converted to matrix
	typedef tsl::robin_map<pair<GroupID, GroupID>, double, HashPair> ResultsHash;


	SeqGroupClassifier(const CountHash &counts, const vector<string> &sampleIDs) :
			m_sampleIDs(sampleIDs) {
		//create reverse hashtable for sampleIDs
		for(SampleID i = 0; i < m_sampleIDs.size(); ++i){
			m_idToSample[m_sampleIDs[i]] = i;
		}

		//parse and load groups file
		string line;
		ifstream gfh(opt::groupingsFile);
		if (gfh.is_open()) {
			while (getline(gfh, line)) {
				std::string delimiter = "\t";
				size_t pos = line.find(delimiter);
				std::string token;

				string sampleName = line.substr(0, pos);

				//move to last part of string (where group ID lives)
				while ((pos = line.find(delimiter)) != std::string::npos) {
				    token = line.substr(0, pos);
//				    std::cout << token << std::endl;
					line.erase(0, pos + delimiter.length());
				}
				string groupName = line;
				tsl::robin_map<string, GroupID>::iterator itr = m_idToGroup.find(groupName);
				if (itr ==  m_idToGroup.end()) {
					m_idToGroup[groupName] = m_groupings.size();
					m_groupings.push_back(shared_ptr<vector<GroupID>>(
							new vector<GroupID>()));
					m_groupIDs.push_back(groupName);
				}
				m_groupings.back()->push_back(m_idToSample[sampleName]);
			}
			gfh.close();
		}
		else {
			cout << "Unable to open file";
		}
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
				if (i->second->at(j) != lastCount) {
					allSame = false;
					break;
				}
			}
			if (!allSame) {
				//if not common to all samples count total number to group and increment size
				//add to index
				m_hashToIndex[i->first] = freqMatSize++;
				//assign count of all k-mers to group
				for (size_t j = 0; j < m_groupings.size(); ++j) {
					for (vector<GroupID>::const_iterator k =
							m_groupings[j]->begin(); k != m_groupings[j]->end();
							++k) {
						totalCounts[j] += i->second->at(*k);
					}
				}
			}
		}
		//set size of vectors
		for (unsigned i = 0; i < m_sampleIDs.size(); ++i) {
			m_groupFreq.push_back(
					shared_ptr<vector<double>>(
							new vector<double>(freqMatSize, 0)));
		}

		//populate matrix
		for (CountHash::const_iterator i = counts.begin(); i != counts.end();
				i++) {
			indexHash::iterator itr = m_hashToIndex.find(i->first);
			//k-mer is not the sample in all samples
			if (itr != m_hashToIndex.end()) {
				for (size_t j = 0; j < m_groupings.size(); ++j) {
					uint64_t countOfGroup = 0;
					for (vector<GroupID>::const_iterator k =
							m_groupings[j]->begin(); k != m_groupings[j]->end();
							++k) {
						countOfGroup += i->second->at(*k);
					}
					double totalCount = double(totalCounts[j])
							+ double(
									m_hashToIndex.size()
											* m_groupings.at(j)->size())
									* opt::pseudoCount;
					double freq = (double(countOfGroup) + opt::pseudoCount)
							/ totalCount;
					(*m_groupFreq[j])[itr->second] = freq;
				}
			}
		}
	}

	tsl::robin_map<GroupID, double> computeAllKLDist(const string &filename){
		tsl::robin_map<GroupID, double> results;
		assert(!filename.empty());
		vector<double> sampleFreq = loadSeqsToFreq(filename);
		//iterate through all combinations
		for(GroupID i= 0; i < m_groupIDs.size(); ++i){
			double result = computeKLDist(
					*m_groupFreq[i], *m_groupFreq[i],
					sampleFreq);
			cerr << m_groupIDs[i] << "\t" << result << endl;
			results[i] = result;
		}
		//for each compute kl distance
		//return results
		return results;
	}

	/*
	 * Runs all combinations of diploid genotypes
	 * Results pairs are cannonically pair1 <= pair2
	 */
	ResultsHash computeDiploidKLDist(
			const string &filename) {
		ResultsHash results;
		assert(!filename.empty());
		vector<double> sampleFreq = loadSeqsToFreq(filename);
		//iterate through all combinations
		for(GroupID i= 0; i < m_groupIDs.size(); ++i){
			for (GroupID j = i; j < m_groupIDs.size(); ++j) {
				double result = computeKLDist(
						*m_groupFreq[i], *m_groupFreq[j],
						sampleFreq);
//				cerr << m_groupIDs[i] << "," << m_groupIDs[j] << "\t"
//						<< m_groupIDs.size() << "\t" << result << endl;
				results[std::make_pair(i, j)] = result;
			}
		}
		//for each compute kl distance
		//return results
		return results;
	}

	void printResults(ResultsHash results){
		pair<GroupID, GroupID> minGroups;
		double minValue = numeric_limits<double>::max();
		for (SeqGroupClassifier::ResultsHash::iterator itr =
				results.begin(); itr != results.end(); ++itr) {
			if(itr->second < minValue){
				minValue = itr->second;
				minGroups = itr->first;
			}
		}
		cout << minValue << "\t" << m_groupIDs[minGroups.first] << "\t"
				<< m_groupIDs[minGroups.second] << endl;
	}

	void printResults(tsl::robin_map<SeqGroupClassifier::GroupID, double> results){
		GroupID minGroup = 0;
		double minValue = numeric_limits<double>::max();
		for (tsl::robin_map<SeqGroupClassifier::GroupID, double>::iterator itr =
				results.begin(); itr != results.end(); ++itr) {
			if(itr->second < minValue){
				minValue = itr->second;
				minGroup = itr->first;
			}
		}
		cout << minValue << "\t" << m_groupIDs[minGroup] << endl;
	}

//	double computeML(const vector<double> &parent1,
//			const vector<double> &parent2,
//			const vector<double> &sampleFreq){
//		double prob = 0.0;
//		for (size_t i = 0; i < sampleFreq.size(); ++i) {
//			double blendedFreq = (parent1[i] + parent2[i]) / 2.0;
//
//		}
//		return(prob);
//
//	}

//	double computeBaysianEstimate(){
//
//	}

	 virtual ~SeqGroupClassifier(){
	 }

private:
	const vector<string> &m_sampleIDs;
	vector<string> m_groupIDs;
	tsl::robin_map<string, SampleID> m_idToSample;
	tsl::robin_map<string, GroupID> m_idToGroup;
	vector<shared_ptr<vector<GroupID>>> m_groupings;
//	tsl::robin_map<SampleID, GroupID> &m_sampleToGroup;
	indexHash m_hashToIndex; //hashed k-mer value to index
	vector<shared_ptr<vector<double>>> m_groupFreq;

	/*
	 * TODO: possible pitfall - currently k-mers are stored as hashvalues rather
	 * than element itself ad hash collisions could artifically add to counts
	 */
	vector<double> loadSeqsToFreq(const string &filename){
		uint64_t rawCoverage = 0; //total number of k-mers
		vector<opt::Count> counts(m_groupFreq[0]->size(),0);
		vector<double> freqs(m_groupFreq[0]->size());

		BloomFilter bf(opt::bf);
		//read in file
		gzFile fp;
		fp = gzopen(filename.c_str(), "r");
		if (fp == Z_NULL) {
			std::cerr << "file " << filename << " cannot be opened"
					<< std::endl;
			exit(1);
		} else if (opt::verbose) {
			std::cerr << "Opening " << filename << std::endl;
		}
		kseq_t *seq = kseq_init(fp);
		int l = kseq_read(seq);
		unsigned index = 0;
		while (l >= 0) {
			//k-merize
			for (ntHashIterator itr(seq->seq.s, opt::hashNum, opt::k);
					itr != itr.end(); ++itr) {
				//remove background k-mers (optional)
				//if kmer exists inside bg set
				if (!bf.contains(*itr)) {
					indexHash::iterator seqIndex = m_hashToIndex.find((*itr)[0]);
					if(seqIndex != m_hashToIndex.end()){
						++counts[seqIndex->second];
						++rawCoverage;
					}
				}
			}
			l = kseq_read(seq);
			index++;
		}
		kseq_destroy(seq);
		gzclose(fp);

		rawCoverage += double(counts.size()) * opt::pseudoCount;

		//divide counts by coverage
		for (IndexPos i = 0; i < counts.size(); ++i) {
			freqs[i] = (double(counts[i]) + opt::pseudoCount)
					/ double(rawCoverage);
		}

		//return vector
		return(freqs);
	}

	/*
	 * Assuming diploid genome compute KL distances metric
	 */
	double computeKLDist(const vector<double> &parent, const vector<double> &sampleFreq) const{
		double dist = 0.0;
		for (size_t i = 0; i < sampleFreq.size(); ++i) {
//			//check for zero values, if present, remove (TODO: Determine if correct way to handle)
			if (parent.at(i) && sampleFreq.at(i)) {
				//using base 2 log, so units use are "bits" rather than "nats" (base e)
				double subDist = sampleFreq.at(i)
						* log2(sampleFreq.at(i) / parent.at(i));
				dist += subDist;
			}
		}
		return (dist);
	}

	/*
	 * Assuming diploid genome compute KL distances metric
	 */
	double computeKLDist(const vector<double> &parent1,
			const vector<double> &parent2, const vector<double> &sampleFreq) const{
		double dist = 0.0;
		for (size_t i = 0; i < sampleFreq.size(); ++i) {
			double blendedFreq = (parent1.at(i) + parent2.at(i)) / 2.0;
//			//check for zero values, if present, remove (TODO: Determine if correct way to handle)
			if (blendedFreq && sampleFreq.at(i)) {
				//using base 2 log, so units use are "bits" rather than "nats" (base e)
				double subDist = sampleFreq.at(i)
						* log2(sampleFreq.at(i) / blendedFreq);
//				cerr << blendedFreq << "\t" << sampleFreq[i] << "\t" << subDist << endl;
				dist += subDist;
			}
		}
		return (dist);
	}
};

#endif /* SRC_SEQGROUPCLASSIFIER_HPP_ */
