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
#include "Util.h"

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
//	typedef tsl::robin_map<pair<SampleID, SampleID>, double, HashPair> ResultsHash;


	SeqGroupClassifier(const CountHash &counts, const vector<string> &sampleIDs) :
			m_sampleIDs(sampleIDs), m_sampleToGroup(vector<GroupID>(m_sampleIDs.size()))  {
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
					m_groupings.push_back(shared_ptr<vector<SampleID>>(
							new vector<SampleID>()));
					m_groupIDs.push_back(groupName);
				}
				m_groupings[m_idToGroup[groupName]]->push_back(m_idToSample[sampleName]);
				m_sampleToGroup[m_idToSample[sampleName]] = m_idToGroup[groupName];
			}
			gfh.close();
		}
		else {
			cout << "Unable to open file";
		}
//		computeGroupFreq(counts);
		computeSampleFreq(counts);
	}

//	/*
//	 * Compute frequency matrix for each of group ID
//	 * Creates a hash value (derived from k-mer) to matrix
//	 */
//	void computeGroupFreq(const CountHash &counts) {
//		//figure out counts of each grouping to calculate frequency
//		vector<uint64_t> totalCounts(m_groupings.size(), 0);
//		uint64_t freqMatSize = 0;
//		//for each k-mer
//		for (CountHash::const_iterator i = counts.begin(); i != counts.end();
//				i++) {
//			//determine if common to all samples
//			uint16_t lastCount = i->second->at(0);
//			bool allSame = true;
//			for (unsigned j = 1; j < m_sampleIDs.size(); ++j) {
//				if (i->second->at(j) != lastCount) {
//					allSame = false;
//					break;
//				}
//			}
//			if (!allSame) {
//				//if not common to all samples count total number to group and increment size
//				//add to index
//				m_hashToIndex[i->first] = freqMatSize++;
//				//assign count of all k-mers to group
//				for (size_t j = 0; j < m_groupings.size(); ++j) {
//					for (vector<GroupID>::const_iterator k =
//							m_groupings[j]->begin(); k != m_groupings[j]->end();
//							++k) {
//						totalCounts[j] += i->second->at(*k);
//					}
//				}
//			}
//		}
//		//set size of vectors
//		for (unsigned i = 0; i < m_groupings.size(); ++i) {
//			m_groupFreq.push_back(
//					shared_ptr<vector<double>>(
//							new vector<double>(freqMatSize, 0)));
//		}
//
//		//populate matrix
//		for (CountHash::const_iterator i = counts.begin(); i != counts.end();
//				i++) {
//			indexHash::iterator itr = m_hashToIndex.find(i->first);
//			//k-mer is not the sample in all samples
//			if (itr != m_hashToIndex.end()) {
//				for (size_t j = 0; j < m_groupings.size(); ++j) {
//					uint64_t countOfGroup = 0;
//					for (vector<GroupID>::const_iterator k =
//							m_groupings[j]->begin(); k != m_groupings[j]->end();
//							++k) {
//						countOfGroup += i->second->at(*k);
//					}
//					double totalCount = double(totalCounts[j])
//							+ double(
//									m_hashToIndex.size()
//											* m_groupings.at(j)->size())
//									* opt::pseudoCount;
//					double freq = (double(countOfGroup) + opt::pseudoCount)
//							/ totalCount;
//					(*m_groupFreq[j])[itr->second] = freq;
//				}
//			}
//		}
//	}
//
	/*
	 * Compute frequency matrix for each of group ID
	 * Creates a hash value (derived from k-mer) to matrix
	 */
	void computeSampleFreq(const CountHash &counts) {
		//figure out counts of each grouping to calculate frequency
		vector<uint64_t> sampleCount(m_sampleIDs.size(), 0);
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
				for (unsigned j = 0; j < m_sampleIDs.size(); ++j) {
					sampleCount[j] += i->second->at(j);
				}
			}
		}
		//set size of vectors
		for (unsigned i = 0; i < m_sampleIDs.size(); ++i) {
			m_sampleFreq.push_back(
					shared_ptr<vector<double>>(
							new vector<double>(freqMatSize, 0)));
			m_sampleCount.push_back(
					shared_ptr<vector<opt::Count>>(
							new vector<opt::Count>(freqMatSize, 0)));
		}

		//populate matrix
		for (CountHash::const_iterator i = counts.begin(); i != counts.end();
				i++) {
			indexHash::iterator itr = m_hashToIndex.find(i->first);
			if (itr != m_hashToIndex.end()) {
				for (unsigned j = 0; j < m_sampleIDs.size(); ++j) {

					double totalCount = double(sampleCount[j])
							+ double(m_hashToIndex.size()) * opt::pseudoCount;
					double freq = (double(i->second->at(j)) + opt::pseudoCount)
							/ totalCount;
					(*m_sampleCount[j])[itr->second] = i->second->at(j);
					(*m_sampleFreq[j])[itr->second] = freq;
				}
			}
		}
	}

//	tsl::robin_map<SampleID, double> computeAllKLDist(const string &filename){
//		tsl::robin_map<SampleID, double> results;
//		assert(!filename.empty());
//		vector<double> sampleFreq = loadSeqsToFreq(filename);
//		//iterate through all combinations
//		for (SampleID i = 0; i < m_sampleIDs.size(); ++i) {
//			double result = computeKLDist(*m_sampleFreq[i], sampleFreq);
//			results[i] = result;
//		}
//		return results;
//	}

//	/*
//	 * Runs all combinations of diploid genotypes
//	 */
//	void computeDiploidKLDist(
//			const string &filename) {
//		assert(Util::fexists(filename));
//		pair<SampleID, SampleID> minGroups;
//		double minValue = numeric_limits<double>::max();
//		vector<double> sampleFreq = loadSeqsToFreq(filename).first;
//		//iterate through all combinations
//		for (SampleID i = 0; i < m_sampleIDs.size(); ++i) {
//			for (SampleID j = i; j < m_sampleIDs.size(); ++j) {
//				double result = computeKLDist(*m_sampleFreq[i],
//						*m_sampleFreq[j], sampleFreq);
//				if (result < minValue) {
//					minValue = result;
//					minGroups = pair<SampleID, SampleID> (i, j);
//				}
////				cout << result << "\t" << m_sampleIDs[i] << "\t"
////						<< m_sampleIDs[j] << "\t" << m_groupIDs[m_sampleToGroup[i]] << "\t"
////						<< m_groupIDs[m_sampleToGroup[j]] << endl;
//			}
//		}
//		cout << minValue << "\t" << m_sampleIDs[minGroups.first] << "\t"
//				<< m_sampleIDs[minGroups.second] << "\t"
//				<< m_groupIDs[m_sampleToGroup[minGroups.first]] << "\t"
//				<< m_groupIDs[m_sampleToGroup[minGroups.second]] << endl;
//	}

	/*
	 * Runs all combinations of diploid genotypes
	 */
	void computeDiploid(const string &filename) {
		assert(Util::fexists(filename));
		vector<opt::Count> sampleCount = loadSeqsToCount(filename);
		vector<double> mleValues(m_sampleIDs.size()*((m_sampleIDs.size()+1))/2);
		size_t matIndex = 0;
		double maxLogValue = -1.0 * numeric_limits<double>::max();

		for (SampleID i = 0; i < m_sampleIDs.size(); ++i) {
			for (SampleID j = i; j < m_sampleIDs.size(); ++j) {
				double logMLE = computeLogMLE(*m_sampleFreq[i],
						*m_sampleFreq[j], sampleCount);
				mleValues[matIndex++] = logMLE;
				if(maxLogValue < logMLE){
					maxLogValue = logMLE;
				}
			}
		}

		double bayesDeom = computeBayesDeomDiploid(mleValues, maxLogValue);
//		cerr << bayesDeom << endl;
		assert(bayesDeom > 0);

		matIndex = 0;

		double minValueKL = numeric_limits<double>::max();
		double maxBayesProb = 0;
		double maxLogMLE = -1.0 * numeric_limits<double>::max();
		pair<SampleID, SampleID> minIDs;

//		double bayesSum = 0;

		//iterate through all combinations
		for (SampleID i = 0; i < m_sampleIDs.size(); ++i) {
			for (SampleID j = i; j < m_sampleIDs.size(); ++j) {
				double klDist = computeKLDist(*m_sampleFreq[i],
						*m_sampleFreq[j], sampleCount);

				double logMLE = mleValues[matIndex++];
				double mle = pow(2, logMLE - maxLogValue);
				double bayesProb = computeBayesianProb(
						1.0 / double(m_sampleIDs.size()),
						1.0 / double(m_sampleIDs.size()), mle, bayesDeom);
//				bayesSum += bayesProb;
//				cout << m_sampleIDs[i] << "\t" << m_sampleIDs[j] << "\t"
//						<< m_groupIDs[m_sampleToGroup[i]] << "\t"
//						<< m_groupIDs[m_sampleToGroup[j]] << "\t" << klDist
//						<< "\t" << minValueKL << "\t" << logMLE << "\t"
//						<< maxLogMLE << "\t" << bayesProb << "\t"
//						<< maxBayesProb << "\t" << mle << "\t" << bayesSum <<  endl;

				if (klDist < minValueKL) {
//					cerr << klDist << "\t" << minValueKL << "\t" << logMLE
//							<< "\t" << maxLogMLE << "\t" << bayesProb << "\t"
//							<< maxBayesProb << "\n";
					minValueKL = klDist;
					minIDs = pair<SampleID, SampleID>(i, j);
					assert(logMLE >= maxLogMLE);
					maxLogMLE = logMLE;
					assert(bayesProb >= maxBayesProb);
					maxBayesProb = bayesProb;
				}

				if(bayesProb > 0.0) {
					cout << m_sampleIDs[i] << "\t" << m_sampleIDs[j] << "\t"
							<< m_groupIDs[m_sampleToGroup[i]] << "\t"
							<< m_groupIDs[m_sampleToGroup[j]] << "\t" << klDist
							<< "\t" <<  logMLE << "\t" << bayesProb <<  endl;
				}
			}
		}

		if(opt::debug){
			//debugfilehandle
			std::ofstream countFH (opt::outputPrefix + ".debug.counts.tsv", std::ofstream::out);
			std::ofstream infoFH (opt::outputPrefix + ".debug.tsv", std::ofstream::out);

			countFH << "sample" << "\t";
			for (vector<opt::Count>::const_iterator itr =
					sampleCount.begin(); itr != sampleCount.end();
					itr++) {
				countFH << *itr << "\t";
			}

			matIndex = 0;
			for (SampleID i = 0; i < m_sampleIDs.size(); ++i) {
				for (SampleID j = i; j < m_sampleIDs.size(); ++j) {
					double klDist = computeKLDist(*m_sampleFreq[i],
							*m_sampleFreq[j], sampleCount);

					double logMLE = mleValues[matIndex++];
					double mle = pow(2, logMLE - maxLogValue);
					double bayesProb = computeBayesianProb(
							1.0 / double(m_sampleIDs.size()),
							1.0 / double(m_sampleIDs.size()), mle, bayesDeom);

					if (bayesProb == maxBayesProb) {
						infoFH << m_sampleIDs[i] << "_" << m_sampleIDs[j]
								<< "\t" << m_sampleIDs[i] << "\t"
								<< m_sampleIDs[j] << "\t"
								<< m_groupIDs[m_sampleToGroup[i]] << "\t"
								<< m_groupIDs[m_sampleToGroup[j]] << "\t"
								<< klDist << "\t" << logMLE << "\t" << bayesProb
								<< endl;
						countFH << m_sampleIDs[i] << "\t";
						printDebug(i, countFH);
						countFH << m_sampleIDs[j] << "\t";
						printDebug(j, countFH);
						countFH << m_sampleIDs[i] << "_" << m_sampleIDs[j] << "\t";
						printDebugBlend(i, j, countFH);
					}
				}
			}
			countFH.close();
			infoFH.close();
		}
	}

//	void printResults(ResultsHash results){
//		pair<GroupID, GroupID> minGroups;
//		double minValue = numeric_limits<double>::max();
//		for (SeqGroupClassifier::ResultsHash::iterator itr =
//				results.begin(); itr != results.end(); ++itr) {
//			if(itr->second < minValue){
//				minValue = itr->second;
//				minGroups = itr->first;
//			}
//		}
//		cout << minValue << "\t" << m_groupIDs[minGroups.first] << "\t"
//				<< m_groupIDs[minGroups.second] << endl;
//	}

//	void printResults(ResultsHash results) {
//		pair<GroupID, GroupID> minGroups;
////		double minValue = numeric_limits<double>::max();
//		for (SeqGroupClassifier::ResultsHash::iterator itr = results.begin();
//				itr != results.end(); ++itr) {
//			cout << itr->second << "\t" << m_sampleIDs[itr->first.first] << "\t"
//					<< m_sampleIDs[itr->first.second] << endl;
////			if(itr->second < minValue){
////				minValue = itr->second;
////				minGroups = itr->first;
////			}
//		}
////		cout << minValue << "\t" << m_groupIDs[minGroups.first] << "\t"
////				<< m_groupIDs[minGroups.second] << endl;
//	}

	void printResults(tsl::robin_map<SeqGroupClassifier::GroupID, double> results){
//		SampleID minGroup = 0;
//		double minValue = numeric_limits<double>::max();
		for (tsl::robin_map<SeqGroupClassifier::GroupID, double>::iterator itr =
				results.begin(); itr != results.end(); ++itr) {
			cout << itr->second << "\t" << m_sampleIDs[itr->first] << endl;
//			if(itr->second < minValue){
//				minValue = itr->second;
//				minGroup = itr->first;
//			}
		}
//		cout << minValue << "\t" << m_sampleIDs[minGroup] << endl;
	}

	/*
	 * Print out k-mer count profile of specific sample
	 */
	void printDebug(SampleID id, ofstream &debugStream) {
		for (size_t i = 0; i < m_sampleCount.at(id)->size(); ++i) {
			debugStream << m_sampleCount.at(id)->at(i) << "\t";
		}
		debugStream << endl;
	}

	/*
	 * Print out k-mer count profile of specific sample
	 */
	void printDebugBlend(SampleID id1, SampleID id2, ofstream &debugStream) {
		for (size_t i = 0; i < m_sampleCount.at(id1)->size(); ++i) {
			debugStream
					<< m_sampleCount.at(id1)->at(i)
							+ m_sampleCount.at(id2)->at(i) << "\t";
		}
		debugStream << endl;
	}

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
	vector<shared_ptr<vector<SampleID>>> m_groupings;
	vector<GroupID> m_sampleToGroup;
//	tsl::robin_map<SampleID, GroupID> &m_sampleToGroup;
	indexHash m_hashToIndex; //hashed k-mer value to index
//	vector<shared_ptr<vector<double>>> m_groupFreq;
	vector<shared_ptr<vector<double>>> m_sampleFreq;
	vector<shared_ptr<vector<opt::Count>>> m_sampleCount; //for debug purposes

//	/*
//	 * TODO: possible pitfall - currently k-mers are stored as hashvalues rather
//	 * than element itself ad hash collisions could artifically add to counts
//	 */
//	pair<vector<double>, vector<opt::Count>> loadSeqsToFreq(const string &filename) {
//		uint64_t rawCoverage = 0; //total number of k-mers
//		vector<opt::Count> counts(m_hashToIndex.size(), 0);
//		vector<double> freqs(m_hashToIndex.size());
//
//		BloomFilter bf(opt::bf);
//		//read in file
//		gzFile fp;
//		fp = gzopen(filename.c_str(), "r");
//		if (fp == Z_NULL) {
//			std::cerr << "file " << filename << " cannot be opened"
//					<< std::endl;
//			exit(1);
//		} else if (opt::verbose) {
//			std::cerr << "Opening " << filename << std::endl;
//		}
//		kseq_t *seq = kseq_init(fp);
//		int l = kseq_read(seq);
//		unsigned index = 0;
//		while (l >= 0) {
//			//k-merize
//			for (ntHashIterator itr(seq->seq.s, opt::hashNum, opt::k);
//					itr != itr.end(); ++itr) {
//				//remove background k-mers (optional)
//				//if kmer exists inside bg set
//				if (!bf.contains(*itr)) {
//					indexHash::iterator seqIndex = m_hashToIndex.find((*itr)[0]);
//					if(seqIndex != m_hashToIndex.end()){
//						++counts[seqIndex->second];
//						++rawCoverage;
//					}
//				}
//			}
//			l = kseq_read(seq);
//			index++;
//		}
//		kseq_destroy(seq);
//		gzclose(fp);
//
//		double totalCount = double(rawCoverage) + double(counts.size()) * opt::pseudoCount;
//
//		//divide counts by coverage
//		for (IndexPos i = 0; i < counts.size(); ++i) {
//			freqs[i] = (double(counts[i]) + opt::pseudoCount)
//					/ totalCount;
//		}
//
//		//return vector
//		return(pair<vector<double>, vector<opt::Count>>(freqs, counts));
//	}

	/*
	 * TODO: possible pitfall - currently k-mers are stored as hashvalues rather
	 * than element itself ad hash collisions could artifically add to counts
	 */
	vector<opt::Count> loadSeqsToCount(const string &filename) {
		vector<opt::Count> counts(m_hashToIndex.size(), 0);

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
					}
				}
			}
			l = kseq_read(seq);
			index++;
		}
		kseq_destroy(seq);
		gzclose(fp);
		return(counts);
	}



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
			const vector<double> &parent2,
			const vector<opt::Count> &sampleCount) const {
		double dist = 0.0;
		for (size_t i = 0; i < sampleCount.size(); ++i) {
			double blendedFreq = (parent1.at(i) + parent2.at(i)) / 2.0;
			double sampleFreq = double(sampleCount.at(i)) / double(sampleCount.size());
			if(sampleFreq != 0){
				double subDist = sampleFreq	* log2(sampleFreq / blendedFreq);
				dist += subDist;
			}
		}
		return (dist);
	}

	double computeLogMLE(const vector<double> &parent1,
			const vector<double> &parent2,
			const vector<opt::Count> &sampleCount) const{
		double logProb = 0;
		for (size_t i = 0; i < sampleCount.size(); ++i) {
			double blendedFreq = (parent1[i] + parent2[i]) / 2.0;
			logProb += log2(blendedFreq) * double(sampleCount[i]);
		}
		return(logProb);
	}

	/*
	 * Assuming each group a sample
	 */
	double computeBayesDeomDiploid(vector<double> mleValues,
			double scaling) const {
		double probAll = 0.0;
		//iterate through all combinations
		for (size_t matIndex = 0; matIndex < mleValues.size(); ++matIndex) {
			double logMLE = mleValues[matIndex];
			double freqGroup1 = 1.0 / double(m_sampleIDs.size());
			double freqGroup2 = 1.0 / double(m_sampleIDs.size());
			double prob = pow(2, logMLE - scaling) * freqGroup1 * freqGroup2;
			probAll += prob;
		}
		return (probAll);
	}

//	double computeBayesianProb(double freqGroup, double mle,
//			double probAll) const {
//		double prob = (mle * freqGroup) / probAll;
//		return (prob);
//	}

	double computeBayesianProb(double freqGroup1, double freqGroup2,
			double mle, double probAll) const {
		double prob = (mle*freqGroup1*freqGroup2) / probAll;
		return (prob);
	}
};

#endif /* SRC_SEQGROUPCLASSIFIER_HPP_ */
