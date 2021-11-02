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
#include <ostream>

#include "Options.h"
#include "Util.h"

#include "vendor/tsl/robin_map.h"
#include "vendor/btl_bloomfilter/BloomFilter.hpp"
#include "vendor/ntHash/ntHashIterator.hpp"
#include <zlib.h>
#include "vendor/kseq.h"

#ifndef KSEQ_INIT_NEW
#define KSEQ_INIT_NEW
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

//	SeqGroupClassifier(const CountHash &counts, const vector<string> &refIDs) :
//			m_refIDs(refIDs), m_refToGroup(vector<GroupID>(m_refIDs.size())), m_refTotalCounts(
//					m_refIDs.size(), 0) {
//		//create reverse hashtable for sampleIDs
//		for(SampleID i = 0; i < m_refIDs.size(); ++i){
//			m_idToRef[m_refIDs[i]] = i;
//		}
//
//		//parse and load groups file
//		string line;
//		ifstream gfh(opt::groupingsFile);
//		if (gfh.is_open()) {
//			while (getline(gfh, line)) {
//				std::string delimiter = "\t";
//				size_t pos = line.find(delimiter);
//				std::string token;
//
//				string sampleName = line.substr(0, pos);
//
//				//move to last part of string (where group ID lives)
//				while ((pos = line.find(delimiter)) != std::string::npos) {
//				    token = line.substr(0, pos);
////				    std::cout << token << std::endl;
//					line.erase(0, pos + delimiter.length());
//				}
//				string groupName = line;
//				tsl::robin_map<string, GroupID>::iterator itr = m_idToGroup.find(groupName);
//				if (itr ==  m_idToGroup.end()) {
//					m_idToGroup[groupName] = m_groupings.size();
//					m_groupings.push_back(shared_ptr<vector<SampleID>>(
//							new vector<SampleID>()));
//					m_groupIDs.push_back(groupName);
//				}
//				m_groupings[m_idToGroup[groupName]]->push_back(m_idToRef[sampleName]);
//				m_refToGroup[m_idToRef[sampleName]] = m_idToGroup[groupName];
//			}
//			gfh.close();
//		}
//		else {
//			cout << "Unable to open file";
//		}
////		computeGroupFreq(counts);
//		computeSampleFreq(counts);
//	}

	SeqGroupClassifier(const CountHash &counts, const vector<string> &refIDs) :
			m_refIDs(refIDs) {
		//create reverse hashtable for sampleIDs
//		for(SampleID i = 0; i < m_refIDs.size(); ++i){
//			m_idToRef[m_refIDs[i]] = i;
//		}
		computeSampleFreq(counts);

		m_refTotalCounts = vector<size_t>(m_collapsedIDs.size(), 0);
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
//			for (unsigned j = 1; j < m_refIDs.size(); ++j) {
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
		//determine if duplicated
		tsl::robin_map<SampleID, string> keepSet;

		for (SampleID i = 0; i < m_refIDs.size(); ++i) {
			keepSet[i] = m_refIDs.at(i);
		}

		for (SampleID i = 0; i < m_refIDs.size(); ++i) {
			for (SampleID j = i + 1; j < m_refIDs.size(); ++j) {

				if (isIdentical(counts, i, j)) {
					cerr << "Collapsing duplicate: " << m_refIDs[i] << "\t" << m_refIDs[j]
							<< endl;
					keepSet[i] += "," + m_refIDs.at(j);
					if(keepSet.find(j) != keepSet.end()){
						keepSet.erase(j);
					}
				}
			}
		}
		if (keepSet.size() != m_refIDs.size()) {
			cerr << "merged:" << endl;
			for (tsl::robin_map<SampleID, string>::iterator itr = keepSet.begin();
					itr != keepSet.end(); ++itr) {
					cerr << itr->second << endl;
			}
		}

		//set size of vectors
		for (unsigned i = 0; i < m_refIDs.size(); ++i) {
			if (keepSet.find(i) != keepSet.end()) {
				m_refToCollapsedID[i] = m_collapsedIDs.size();
				m_refTotalCounts.push_back(0);
				m_collapsedIDs.push_back(keepSet.at(i));
			}
		}

		//figure out counts of each grouping to calculate frequency
		uint64_t freqMatSize = 0;
		//for each k-mer
		for (CountHash::const_iterator i = counts.begin(); i != counts.end();
				i++) {
			//determine if common to all samples
			uint16_t lastCount = i->second->at(0);
			bool allSame = true;
			for (unsigned j = 1; j < m_refIDs.size(); ++j) {
				if (i->second->at(j) != lastCount) {
					allSame = false;
					break;
				}
			}
			if (!allSame) {
				m_hashToIndex[i->first] = freqMatSize++;
				//if not common to all samples count total number to group and increment size
				//add to index
				//assign count of all k-mers to group
				for (tsl::robin_map<SampleID, string>::iterator j =
						keepSet.begin(); j != keepSet.end(); ++j) {
					m_refTotalCounts[m_refToCollapsedID[j->first]] +=
							i->second->at(j->first);
				}
			}
		}

		for (tsl::robin_map<SampleID, SampleID>::iterator itr =
				m_refToCollapsedID.begin(); itr != m_refToCollapsedID.end(); ++itr) {
			m_kmerCounts.push_back(
					shared_ptr<vector<opt::Count>>(
							new vector<opt::Count>(freqMatSize, 0)));
		}

		//populate matrix
		for (CountHash::const_iterator i = counts.begin(); i != counts.end();
				i++) {
			indexHash::iterator itr = m_hashToIndex.find(i->first);
			if (itr != m_hashToIndex.end()) {
				for (unsigned j = 0; j < m_kmerCounts.size(); ++j) {
					(*m_kmerCounts[j])[itr->second] = i->second->at(j);
				}
			}
		}
	}

//	tsl::robin_map<SampleID, double> computeAllKLDist(const string &filename){
//		tsl::robin_map<SampleID, double> results;
//		assert(!filename.empty());
//		vector<double> sampleFreq = loadSeqsToFreq(filename);
//		//iterate through all combinations
//		for (SampleID i = 0; i < m_refIDs.size(); ++i) {
//			double result = computeKLDist(*m_refFreq[i], sampleFreq);
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
//		for (SampleID i = 0; i < m_refIDs.size(); ++i) {
//			for (SampleID j = i; j < m_refIDs.size(); ++j) {
//				double result = computeKLDist(*m_refFreq[i],
//						*m_refFreq[j], sampleFreq);
//				if (result < minValue) {
//					minValue = result;
//					minGroups = pair<SampleID, SampleID> (i, j);
//				}
////				cout << result << "\t" << m_refIDs[i] << "\t"
////						<< m_refIDs[j] << "\t" << m_groupIDs[m_refToGroup[i]] << "\t"
////						<< m_groupIDs[m_refToGroup[j]] << endl;
//			}
//		}
//		cout << minValue << "\t" << m_refIDs[minGroups.first] << "\t"
//				<< m_refIDs[minGroups.second] << "\t"
//				<< m_groupIDs[m_refToGroup[minGroups.first]] << "\t"
//				<< m_groupIDs[m_refToGroup[minGroups.second]] << endl;
//	}

	/*
	 * Outputs k-mer counts for sample being evaluated
	 */
	void outputSeqCounts(const string &filename){
		assert(Util::fexists(filename));
		vector<opt::Count> sampleCount = loadSeqsToCount(filename);
	}

	/*
	 * Runs all combinations of diploid genotypes
	 */
	void computeDiploid(const string &filename) {
		assert(Util::fexists(filename));

		vector<opt::Count> sampleCount = loadSeqsToCount(filename);

		vector<double> mleValues(m_collapsedIDs.size()*((m_collapsedIDs.size()+1))/2);
		size_t matIndex = 0;
		double maxLogValue = -1.0 * numeric_limits<double>::max();

		size_t sampleCountTotal = 0;

		for (size_t i = 0; i < sampleCount.size(); ++i) {
			sampleCountTotal += sampleCount.at(i);
		}

		for (SampleID i = 0; i < m_collapsedIDs.size(); ++i) {
			for (SampleID j = i; j < m_collapsedIDs.size(); ++j) {
				double logMLE = computeLogMLE(i,j, sampleCount);
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
		for (SampleID i = 0; i < m_collapsedIDs.size(); ++i) {
			for (SampleID j = i; j < m_collapsedIDs.size(); ++j) {
				double klDist = computeKLDist(i, j, sampleCount,
						sampleCountTotal);

				double logMLE = mleValues[matIndex++];
				double mle = pow(2, logMLE - maxLogValue);
				double bayesProb = computeBayesianProb(
						1.0 / double(m_collapsedIDs.size()),
						1.0 / double(m_collapsedIDs.size()), mle, bayesDeom);
//				bayesSum += bayesProb;
//				cout << m_refIDs[i] << "\t" << m_refIDs[j] << "\t"
//						<< m_groupIDs[m_refToGroup[i]] << "\t"
//						<< m_groupIDs[m_refToGroup[j]] << "\t" << klDist
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
//					cout << m_refIDs[i] << "\t" << m_refIDs[j] << "\t"
//							<< m_groupIDs[m_refToGroup[i]] << "\t"
//							<< m_groupIDs[m_refToGroup[j]] << "\t" << klDist
//							<< "\t" <<  logMLE << "\t" << bayesProb <<  endl;
					cout << m_collapsedIDs[i] << "\t" << m_collapsedIDs[j] << "\t" << klDist
							<< "\t" << logMLE << "\t" << bayesProb << endl;
				}
			}
		}

		std::ofstream countFH (opt::outputPrefix + ".debug.counts.tsv", std::ofstream::out);
		countFH << opt::outputPrefix << "sample";

		if(opt::debug){
			//debugfilehandle
			std::ofstream countFH (opt::outputPrefix + ".debug.counts.tsv", std::ofstream::out);
			std::ofstream infoFH (opt::outputPrefix + ".debug.tsv", std::ofstream::out);

			countFH << opt::outputPrefix << "sample";
			for (vector<opt::Count>::const_iterator itr =
					sampleCount.begin(); itr != sampleCount.end();
					itr++) {
				countFH << "\t" << *itr;
			}
			countFH << endl;


			matIndex = 0;
			for (SampleID i = 0; i < m_collapsedIDs.size(); ++i) {
				countFH << m_collapsedIDs[i];
				printDebug(i, countFH);
//				for (SampleID j = i; j < m_collapsedIDs.size(); ++j) {
////					if (opt::outputPrefix == m_refIDs[i] + "_" + m_refIDs[j]) {
////						computeKLDistDebug(i, j, sampleCount, sampleCountTotal);
////					}
//
//					double klDist = computeKLDist(i, j, sampleCount, sampleCountTotal);
//
//					double logMLE = mleValues[matIndex++];
//					double mle = pow(2, logMLE - maxLogValue);
//					double bayesProb = computeBayesianProb(
//							1.0 / double(m_collapsedIDs.size()),
//							1.0 / double(m_collapsedIDs.size()), mle, bayesDeom);
//
//					if (bayesProb > 0) {
////						infoFH << m_refIDs[i] << "_" << m_refIDs[j]
////								<< "\t" << m_refIDs[i] << "\t"
////								<< m_refIDs[j] << "\t"
////								<< m_groupIDs[m_refToGroup[i]] << "\t"
////								<< m_groupIDs[m_refToGroup[j]] << "\t"
////								<< klDist << "\t" << logMLE << "\t" << bayesProb
////								<< endl;
//						infoFH << m_collapsedIDs[i] << "_" << m_collapsedIDs[j]
//								<< "\t" << m_collapsedIDs[i] << "\t"
//								<< m_collapsedIDs[j] << "\t"
//								<< klDist << "\t" << logMLE << "\t" << bayesProb
//								<< endl;
////						cerr << m_refCount.at(i)->size() << endl;
//						countFH << m_collapsedIDs[i] << "_" << m_collapsedIDs[j];
//						printDebugBlend(i, j, countFH);
//					}
//				}
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
//			cout << itr->second << "\t" << m_refIDs[itr->first.first] << "\t"
//					<< m_refIDs[itr->first.second] << endl;
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
			cout << itr->second << "\t" << m_collapsedIDs[itr->first] << endl;
//			if(itr->second < minValue){
//				minValue = itr->second;
//				minGroup = itr->first;
//			}
		}
//		cout << minValue << "\t" << m_refIDs[minGroup] << endl;
	}

	/*
	 * Print out k-mer count profile of specific sample
	 */
	bool isIdentical(const CountHash &counts, SampleID i1, SampleID i2) const {
		bool isIdent = true;
		for (CountHash::const_iterator i = counts.begin(); i != counts.end();
				i++) {
			if (i->second->at(i1) != i->second->at(i2)) {
				isIdent = false;
				break;
			}
		}
		return isIdent;
	}

	/*
	 * Print out k-mer count profile of specific sample
	 */
	void printDebug(SampleID id, ofstream &debugStream) {
		for (size_t i = 0; i < m_kmerCounts.at(id)->size(); ++i) {
			debugStream << "\t" << m_kmerCounts.at(id)->at(i) ;
		}
		debugStream << endl;
	}

//	/*
//	 * Print out k-mer count profile of specific sample
//	 */
//	void printDebugBlend(SampleID id1, SampleID id2, ofstream &debugStream) {
//		for (size_t i = 0; i < m_refCount.at(id1)->size(); ++i) {
//			debugStream << "\t"
//					<< m_refCount.at(id1)->at(i)
//							+ m_refCount.at(id2)->at(i);
//		}
//		debugStream << endl;
//	}

	 virtual ~SeqGroupClassifier(){
	 }

private:
	const vector<string> &m_refIDs;
	vector<string> m_collapsedIDs;
	tsl::robin_map<SampleID, SampleID> m_refToCollapsedID;
//	vector<string> m_groupIDs;
//	tsl::robin_map<string, SampleID> m_idToRef;
//	tsl::robin_map<string, GroupID> m_idToGroup;
//	vector<shared_ptr<vector<SampleID>>> m_groupings;
//	vector<GroupID> m_refToGroup;
//	tsl::robin_map<SampleID, GroupID> &m_refToGroup;
	indexHash m_hashToIndex; //hashed k-mer value to index
//	vector<shared_ptr<vector<double>>> m_groupFreq;
//	vector<shared_ptr<vector<double>>> m_refFreq;
	vector<shared_ptr<vector<opt::Count>>> m_kmerCounts;
	vector<size_t> m_refTotalCounts;

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
#pragma omp parallel
		for (int l;;) {
#pragma omp critical(kseq_read)
			{
				l = kseq_read(seq);
			}
			if (l >= 0) {
				for (ntHashIterator itr(seq->seq.s, opt::hashNum, opt::k);
						itr != itr.end(); ++itr) {
					//if kmer exists inside bg set
					if (!bf.contains(*itr)) {
						indexHash::iterator seqIndex = m_hashToIndex.find((*itr)[0]);
						if(seqIndex != m_hashToIndex.end()){
#pragma omp atomic
							++counts[seqIndex->second];
						}
					}
				}
			} else
				break;
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
	double computeKLDist(SampleID parent1, SampleID parent2,
			const vector<opt::Count> &sampleCount, size_t countTotal) const {
		double dist = 0.0;
		double adjCountTotal = double(countTotal) + double(sampleCount.size()) * opt::pseudoCount;
		double refCountTotal = computeRefTotal(parent1, parent2);
		for (size_t i = 0; i < sampleCount.size(); ++i) {
			double blendedFreq = (double(m_kmerCounts.at(parent1)->at(i))
					+ double(m_kmerCounts.at(parent2)->at(i)) + opt::pseudoCount)
					/ refCountTotal;
			double sampleFreq = double(sampleCount.at(i)) / adjCountTotal;
			if(sampleFreq != 0){
				double subDist = sampleFreq	* log2(sampleFreq / blendedFreq);
				dist += subDist;
			}
		}
		return (dist);
	}

//	double computeKLDistDebug(SampleID ref1, SampleID ref2,
//			const vector<opt::Count> &sampleCount, size_t countTotal) const {
//		double dist = 0.0;
//		double adjCountTotal = double(countTotal) + double(sampleCount.size()) * opt::pseudoCount;
//		double adjCountTotal2 = double(
//				m_refTotalCounts.at(ref1) + m_refTotalCounts.at(ref2))
//				+ double(sampleCount.size()) * opt::pseudoCount;
//		for (size_t i = 0; i < sampleCount.size(); ++i) {
//			double blendedFreq = (m_refFreq.at(ref1)->at(i)
//					+ m_refFreq.at(ref2)->at(i)) / 2.0;
//			double blendedFreq2 = (double(
//					m_refCount.at(ref1)->at(i) + m_refCount.at(ref2)->at(i))
//					+ opt::pseudoCount) / adjCountTotal2;
//			double blendedFreqA = (double(m_refCount.at(ref1)->at(i))
//					+ opt::pseudoCount)
//					/ (double(m_refTotalCounts.at(ref1))
//							+ double(sampleCount.size()) * opt::pseudoCount);
//			double blendedFreqB = (double(m_refCount.at(ref2)->at(i))
//					+ opt::pseudoCount)
//					/ (double(m_refTotalCounts.at(ref2))
//							+ double(sampleCount.size()) * opt::pseudoCount);
//			double sampleFreq = (double(sampleCount.at(i)) + opt::pseudoCount)
//					/ adjCountTotal;
//
//			cout << blendedFreq << "\t" << sampleFreq << "\t" << blendedFreq2
//					<< "\t" << blendedFreqA << "\t" << blendedFreqB << "\t"
//					<< m_refFreq.at(ref1)->at(i) << "\t"
//					<< m_refFreq.at(ref2)->at(i) << "\t" << sampleCount.at(i)
//					<< "\t" << m_refCount.at(ref1)->at(i) << "\t"
//					<< m_refCount.at(ref2)->at(i) << "\t" << countTotal
//					<< "\t" << m_refTotalCounts.at(ref1) << "\t"
//					<< m_refTotalCounts.at(ref2) << endl;
//			if (sampleFreq != 0) {
//				double subDist = sampleFreq	* log2(sampleFreq / blendedFreq);
//				dist += subDist;
//			}
//		}
//		return (dist);
//	}

	double computeRefTotal(SampleID parent1, SampleID parent2) const{
		double refCountTotal = double(
				m_refTotalCounts.at(parent1) + m_refTotalCounts.at(parent2))
				+ double(m_kmerCounts.at(parent1)->size()) * opt::pseudoCount;
		return refCountTotal;
	}

	double computeLogMLE(SampleID parent1, SampleID parent2,
			const vector<opt::Count> &sampleCount) const {
		double logProb = 0;
		double refCountTotal = computeRefTotal(parent1, parent2);
		for (size_t i = 0; i < sampleCount.size(); ++i) {
			double blendedFreq = (double(m_kmerCounts.at(parent1)->at(i))
					+ double(m_kmerCounts.at(parent2)->at(i)) + opt::pseudoCount)
					/ refCountTotal;
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
			double freqGroup1 = 1.0 / double(m_collapsedIDs.size());
			double freqGroup2 = 1.0 / double(m_collapsedIDs.size());
			double prob = pow(2, logMLE - scaling) * freqGroup1 * freqGroup2;
			probAll += prob;
		}
		return (probAll);
	}

	double computeBayesianProb(double freqGroup1, double freqGroup2,
			double mle, double probAll) const {
		double prob = (mle*freqGroup1*freqGroup2) / probAll;
		return (prob);
	}
};

#endif /* SRC_SEQGROUPCLASSIFIER_HPP_ */
