/*
 *  BEANDisco: functions for computing scores
 *  
 *  Copyright 2011-2015 Teppo Niinimäki <teppo.niinimaki(at)helsinki.fi>
 *  
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#include <cmath>

#include "data.hpp"
//#include "stacksubset.hpp"

#ifndef SCORES_HPP
#define SCORES_HPP


class ScoreFun {
public:
	virtual double compute(int nValues, int nParentValues, int* counts) const = 0;
	virtual ~ScoreFun() {};
};


/**
 * BDeu score function.
 */
class BDeuScore : public ScoreFun {
private:
	double ess_;
public:
	BDeuScore(double ess) {
		ess_ = ess;
	}
	double compute(int nValues, int nParentValues, int* counts) const {
		double score = 0;
		double pseudocount = ess_ / (nValues * nParentValues);
		double parentPseudocount = ess_ / nParentValues;
		for (int pv = 0; pv < nParentValues; ++pv) {
			int cumCount = 0;
			for (int v = 0; v < nValues; ++v) {
				int c = counts[pv * nValues + v];
				score += lgamma(c + pseudocount) - lgamma(pseudocount);
				cumCount += c;
			}
			score += lgamma(parentPseudocount) - lgamma(cumCount + parentPseudocount);
		}
		return score;
	}
};

/**
 * BDeu score function with precomputed log gammas.
 */
class BDeuScoreCached : public ScoreFun {
private:
	double ess_;
	double* logGammas_;
	int maxCount_;
public:
	BDeuScoreCached(double ess, int maxCount, int maxNumValues) {
		ess_ = ess;
		maxCount_ =  maxCount;
		logGammas_ = new double[(maxCount_ + 1) * maxNumValues];
		for (int i = 0; i <= maxCount_; ++i)
			for (int j = 0; j < maxNumValues; ++j)
				logGammas_[i + j * (maxCount_ + 1)] = lgamma(i + ess_ / (j+1));
	}
	
	~BDeuScoreCached() {
		delete[] logGammas_;
	}
	
	double compute(int nValues, int nParentValues, int* counts) const {
		double score = 0;
		int pseudoTermIndex = (nParentValues * nValues - 1) * (maxCount_ + 1);
		int parentPseudoTermIndex = (nParentValues - 1) * (maxCount_ + 1);
		
		for (int pv = 0; pv < nParentValues; ++pv) {
			int cumCount = 0;
			for (int v = 0; v < nValues; ++v) {
				int c = counts[pv * nValues + v];
				score += logGammas_[c + pseudoTermIndex] - logGammas_[pseudoTermIndex];
				cumCount += c;
			}
			score += logGammas_[parentPseudoTermIndex] - logGammas_[cumCount + parentPseudoTermIndex];
		}
		return score;
	}
};


/**
 * K2 score function.
 */
class K2Score : public ScoreFun {
public:
	double compute(int nValues, int nParentValues, int* counts) const {
		double score = 0;
		for (int pv = 0; pv < nParentValues; ++pv) {
			int cumCount = 0;
			for (int v = 0; v < nValues; ++v) {
				int c = counts[pv * nValues + v];
				score += lgamma(c + 1);
				cumCount += c;
			}
			score += lgamma(nValues) - lgamma(cumCount + nValues);
		}
		return score;
	}
};

/**
 * K2 score function with precomputed log gammas.
 */
class K2ScoreCached : public ScoreFun {
private:
	double* logGammas_;
public:
	K2ScoreCached(int maxCount) {
		logGammas_ = new double[maxCount + 1];
		for (int i = 0; i <= maxCount; ++i)
			logGammas_[i] = lgamma(i);
	}
	
	~K2ScoreCached() {
		delete[] logGammas_;
	}
	
	double compute(int nValues, int nParentValues, int* counts) const {
		double score = 0;
		for (int pv = 0; pv < nParentValues; ++pv) {
			int cumCount = 0;
			for (int v = 0; v < nValues; ++v) {
				int c = counts[pv * nValues + v];
				score += logGammas_[c + 1];
				cumCount += c;
			}
			score += logGammas_[nValues] - logGammas_[cumCount + nValues];
		}
		return score;
	}
};



double computeLLScore(int nValues, int nParentValues, int* counts) {
	double score = 0;
	for (int pv = 0; pv < nParentValues; ++pv) {
		int cumCount = 0;
		for (int v = 0; v < nValues; ++v) {
			int c = counts[pv * nValues + v];
			if (c > 0)
				score += c * log(c);
			cumCount += c;
		}
		if (cumCount > 0)
			score -= cumCount * log(cumCount);
	}
	return score;
}

/**
 * LL (log-likelihood) score function.
 */
class LLScore : public ScoreFun {
public:
	double compute(int nValues, int nParentValues, int* counts) const {
		return computeLLScore(nValues, nParentValues, counts);
	}
};


/**
 * MDL (minimum description length) score function.
 */
class MDLScore : public ScoreFun {
private:
	int nSamples_;
public:
	MDLScore(int nSamples) {
		nSamples_ = nSamples;
	}
	double compute(int nValues, int nParentValues, int* counts) const {
		return computeLLScore(nValues, nParentValues, counts) - .5 * log(nSamples_) * (nValues - 1) * nParentValues;
	}
};


/**
 * AIC (Akaike Information Criterion) score function.
 */
class AICScore : public ScoreFun {
public:
	double compute(int nValues, int nParentValues, int* counts) const {
		return computeLLScore(nValues, nParentValues, counts) - (nValues - 1) * nParentValues;
	}
};


/*
class CountTable {
private:
	int* counts_;

public:
	int nParentConfs;
	int nNodeConfs;
	
//	CountTable(int _nParentConfs, int _nNodeConfs)
//			: nParentConfs(_nParentConfs), nNodeConfs(_nNodeConfs) {
//		counts_ = new int[nParentConfs * nNodeConfs];
//	}
	
	CountTable(const DataView& dataView, const StackSubset& parents, int node) {
		// compute numbers of configurations
		nParentConfs = 1;
		for (int i = 0; i < parents.size(); ++i)
			nParentConfs *= dataView.getArity(parents[i]);
		nNodeConfs = dataView.getArity(node);
		
		counts_ = new int[nParentConfs * nNodeConfs];
		// get counts
		std::vector<int> vars(parents.size() + 1);
		for (int i = 0; i < parents.size(); ++i)
			vars[i] = parents[i];
		vars[parents.size()] = node;
		int* counts = dataView->getCounts(vars);
	}
	
	~CountTable() {
		delete[] counts_;
	}
	
//	void setCounts(int* counts) {
//		for (int i = 0; i < nParentConfs * nNodeConfs; ++i)
//			counts_[i] = counts[i];
//	}
	
	int operator()(int pa, int i) const {
		return counts_[pa * nNodeConfs + i];
	}

	int& operator()(int pa, int i) {
		return counts_[pa * nNodeConfs + i];
	}
};/**/


template <typename Subset>
double computeScore(const DataView* dataView, const Subset& parents, int node, const ScoreFun* scoreFun) {
	// get counts
	std::vector<int> vars(parents.size() + 1);
	for (int i = 0; i < parents.size(); ++i)
		vars[i] = parents[i];
	vars[parents.size()] = node;
	int* counts = dataView->getCounts(vars);
	
	// compute score
	int nParentValues = 1;
	for (int i = 0; i < parents.size(); ++i)
		nParentValues *= dataView->getArity(parents[i]);
	int nNodeValues = dataView->getArity(node);
	double score = scoreFun->compute(nNodeValues, nParentValues, counts);
	delete[] counts;
	return score;
}/**/



#endif

