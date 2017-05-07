/*
 *  BEANDisco: ADTree implementation
 *  
 *  Copyright 2012-2015 Teppo Niinimäki <teppo.niinimaki(at)helsinki.fi>
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


#include <list>
#include <algorithm>

#include "data.hpp"

#ifndef ADTREE_HPP
#define ADTREE_HPP


class ADTree : public DataView {
private:
	struct ADNode {
		int count;
		union {
			ADNode*** varyNodes;
			int* records;
		};
	};

	ADNode* makeADNode(int depth, int i, const DataColumns& data, const std::list<int>& records) {
		ADNode* adNode = new ADNode();
		adNode->count = records.size();
		if (depth == maxSetSize_) {
			// max depth
			adNode->varyNodes = nullptr;
		} else if (adNode->count < minCount_ || depth >= maxDepth_) {
			// build leaf list
			adNode->records = new int[records.size()];
			int r = 0;
			for (std::list<int>::const_iterator it = records.begin(); it != records.end(); ++it, ++r)
				adNode->records[r] = *it;
		} else {
			// continue recursion
			adNode->varyNodes = new ADNode**[nVariables_ - i];
			for (int j = i; j < nVariables_; ++j)
				adNode->varyNodes[nVariables_-1-j] = makeVaryNode(depth, j, data, records);
		}
		++nAdNodes;
		return adNode;
	}
	
	ADNode** makeVaryNode(int depth, int i, const DataColumns& data, const std::list<int>& records) {
		ADNode** varyNode = new ADNode*[arities_[i]];
		// partition records based on the value of the variable
		std::vector<std::list<int> > childRecords(arities_[i]);
		for (std::list<int>::const_iterator it = records.begin(); it != records.end(); ++it) {
			Datum val = data(i, *it);
			childRecords[val].push_back(*it);
		}
		// continue recursively for each value separately
		for (int val = 0; val < arities_[i]; ++val) {
			if (childRecords[val].empty())
				varyNode[val] = nullptr;
			else
				varyNode[val] = makeADNode(depth+1, i+1, data, childRecords[val]);
		}
		return varyNode;
	}
	
	void freeADNode(ADNode* adNode, int depth, int i) {
		if (depth == maxSetSize_) {
			// do nothing
		} else if (adNode->count < minCount_ || depth >= maxDepth_) {
			// leaf list
			delete[] adNode->records;
		} else {
			// recursion continues
			for (int j = i; j < nVariables_; ++j)
				freeVaryNode(adNode->varyNodes[nVariables_-1-j], depth, j);
			delete[] adNode->varyNodes;
		}
		delete adNode;
	}
	
	void freeVaryNode(ADNode** varyNode, int depth, int i) {
		for (int val = 0; val < arities_[i]; ++val) {
			if (varyNode[val] != nullptr)
				freeADNode(varyNode[val], depth+1, i+1);
		}
		delete[] varyNode;
	}
	
	void fillCounts(const ADNode* adNode, const std::vector<int>& vars,
			int i, const std::vector<int>& cumArities, int* counts) const {
		if (i >= vars.size()) {
			*counts = adNode->count;
			return;
		}
		assert(vars[i] >= 0 && vars[i] < nVariables_);
		if (adNode->count < minCount_ || i >= maxDepth_) {
			// use leaf list
			for (int r = 0; r < adNode->count; ++r) {
				int rec = adNode->records[r];
				int index = 0;
				for (int j = i; j < vars.size(); ++j)
					index += (*data_)(vars[j], rec) * cumArities[j];
				++counts[index];
			}
		} else {
			// continue recursion
			ADNode** varyNode = adNode->varyNodes[nVariables_ - 1 - vars[i]];
			for (int val = 0; val < arities_[vars[i]]; ++val) {
				if (varyNode[val]) {
					fillCounts(varyNode[val], vars, i+1, cumArities, counts+val*cumArities[i]);
				} /*else {
					// filled with zero at initialization
				}*/
			}
		}
	}
	
	struct Cmp {
		const std::vector<int>& vars_;
		Cmp(const std::vector<int>& vars) : vars_(vars) {}
		bool operator()(int i, int j) {return vars_[i] < vars_[j];}
	};
	
	ADNode* root_;
	int nVariables_;
	int* arities_;
	
	const DataColumns* data_;
	int minCount_;
	int maxDepth_;
	int maxSetSize_;
	
	int nAdNodes;
	
public:
	
	ADTree(const DataColumns& data, int minCount = 0, int maxDepth = 0, int maxSetSize = 0) {
		nVariables_ = data.getNumVariables();
		minCount_ = minCount;
		maxDepth_ = (maxDepth <= 0 ? nVariables_ : maxDepth);
		maxSetSize_ = (maxSetSize <= 0 ? nVariables_ : maxSetSize);
		data_ = nullptr;
		if (minCount_ > 0 || maxDepth_ < nVariables_)
			data_ = new DataColumns(data);
		arities_ = new int[nVariables_];
		for (int i = 0; i < nVariables_; ++i)
			arities_[i] = data.getArity(i);
		std::list<int> records;
		for (int i = 0; i < data.getNumSamples(); ++i)
			records.push_back(i);
		nAdNodes = 0;
		root_ = makeADNode(0, 0, data, records);
	}
	
	~ADTree() {
		freeADNode(root_, 0, 0);
		if (data_)
			delete data_;
		delete[] arities_;
	}

	int getNumSamples() const {
		return root_->count;
	}
	
	int getNumVariables() const {
		return nVariables_;
	}
	
	int getArity(int i) const {
		return arities_[i];
	}
	
	int* getCounts(const std::vector<int>& vars) const {
		int n = vars.size();
		assert(n <= maxSetSize_);
		
		// get the order of sorted variables
		std::vector<int> order(n);
		for (int i = 0; i < n; ++i)
			order[i] = i;
		Cmp cmp(vars);
		sort(order.begin(), order.end(), cmp);
		
		// compute cumulative arities
		std::vector<int> cumArities(n);
		cumArities[n-1] = 1;
		for (int i = n-1; i > 0; --i)
			cumArities[i-1] = cumArities[i] * arities_[vars[i]];
		int nValues = cumArities[0] * arities_[vars[0]];
		
		// get the sorted variables and arities
		std::vector<int> sortedVars(n);
		std::vector<int> sortedCumArities(n);
		for (int i = 0; i < n; ++i) {
			sortedVars[i] = vars[order[i]];
			sortedCumArities[i] = cumArities[order[i]];
			assert(i == 0 || sortedVars[i-1] < sortedVars[i]);
		}

		// create the count array and initialize to zero
		int* counts = new int[nValues];
		for (int i = 0; i < nValues; ++i)
			counts[i] = 0;
		
		// fill the array recursively
		fillCounts(root_, sortedVars, 0, sortedCumArities, counts);
		
		return counts;
	}
	
	int getNumADNodes() const {
		return nAdNodes;
	}
	
};

#endif

