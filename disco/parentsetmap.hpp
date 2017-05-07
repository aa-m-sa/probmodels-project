/*
 *  BEANDisco: parentsetmap class
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

#include "common.hpp"
//#include "stacksubset.hpp"
#include "subsetmap.hpp"

#include "dag.hpp"

#ifndef PARENTSETMAP_HPP
#define PARENTSETMAP_HPP

template <class T>
class ParentSetMap {
public:
	const int nNodes;
	const int maxParents;
private:
	SubsetMap<T>** maps_;

	ParentSetMap(const ParentSetMap&) = delete; // disable copying
	ParentSetMap& operator=(const ParentSetMap&) = delete; // disable copying

public:
	ParentSetMap(int n, int k) :
			nNodes(n), maxParents(k)
	{
		maps_ = new SubsetMap<T>*[nNodes];
		//std::vector<int> potParents;
		//potParents.reserve(nNodes - 1);
		SortedArraySubset potParents(nNodes - 1);
		for (int node = 1; node < nNodes; ++node)
			potParents.insertLargest(node);
			//potParents.push_back(node);
		for (int node = 0; node < nNodes; ++node) {
			maps_[node] = new SubsetMap<T>(potParents, maxParents);
			if (node + 1 < nNodes)
				potParents[node] -= 1;
		}
	}
	
	~ParentSetMap() {
		for (int i = 0; i < nNodes; ++i)
			delete maps_[i];
		delete[] maps_;
	}
	
	SubsetMap<T>& operator[](int node) {
		assert(0 <= node && node < nNodes);
		return *maps_[node];
	}

	const SubsetMap<T>& operator[](int node) const {
		assert(0 <= node && node < nNodes);
		return *maps_[node];
	}

	SubsetMap<T>* forNode(int node) {
		assert(0 <= node && node < nNodes);
		return maps_[node];
	}

	const SubsetMap<T>* forNode(int node) const {
		assert(0 <= node && node < nNodes);
		return maps_[node];
	}

	class ParentSets {
	private:
		const ParentSetMap& psMap_;
		std::unique_ptr<typename SubsetMap<T>::SubsetIndex[]> parentSetIndices_;

		std::unique_ptr<typename SubsetMap<T>::SubsetIndex[]> allocParentSetIndices(size_t n) {
			return std::unique_ptr<typename SubsetMap<T>::SubsetIndex[]>(
					static_cast<typename SubsetMap<T>::SubsetIndex*> (
						::operator new[] (n * sizeof(typename SubsetMap<T>::SubsetIndex))));
		}
	public:
		ParentSets(const ParentSetMap& psMap) :
			psMap_(psMap),
			parentSetIndices_(allocParentSetIndices(psMap_.nNodes))
		{
			auto n = psMap_.nNodes;
			for (int node = 0; node < n; ++node)
				parentSetIndices_[node] = psMap_[node].getEmptyIndex();
		}
		ParentSets(const ParentSetMap& psMap, const DagFamily::Instance& dag) :
			psMap_(psMap),
			parentSetIndices_(allocParentSetIndices(psMap_.nNodes))
		{
			assert(psMap_.nNodes == dag.family.n);
			auto n = psMap_.nNodes;
			for (int node = 0; node < n; ++node)
				parentSetIndices_[node] = psMap_[node].getIndex(dag.getParents(node));
		}
		ParentSets(const ParentSets& other) :
			psMap_(other.psMap_),
			parentSetIndices_(allocParentSetIndices(psMap_.nNodes))
		{
			auto n = psMap_.nNodes;
			for (int node = 0; node < n; ++node)
				parentSetIndices_[node] = other.parentSetIndices_[node];
		}
		ParentSets(ParentSets&& other) :
			psMap_(other.psMap_),
			parentSetIndices_(std::move(other.parentSetIndices_))
		{
		}
		~ParentSets() {
		}
		ParentSets& operator=(const ParentSets& other) {
			assert(&psMap_ == &(other.psMap_));
			auto n = psMap_.nNodes;
			for (int node = 0; node < n; ++node)
				parentSetIndices_[node] = other.parentSetIndices_[node];
			return *this;
		}
		ParentSets& operator=(ParentSets&& other) {
			assert(&psMap_ == &(other.psMap_));
			parentSetIndices_ = std::move(other.parentSetIndices_);
			return *this;
		}
		ParentSets& operator=(const DagFamily::Instance& dag) {
			assert(psMap_.nNodes == dag.family.n);
			auto n = psMap_.nNodes;
			for (int node = 0; node < n; ++node)
				parentSetIndices_[node] = psMap_[node].getIndex(dag.getParents(node));
			return *this;
		}
		bool operator==(const ParentSets& other) const {
			if (&psMap_ != &(other.psMap_))
				return false;
			for (int i = 0; i < psMap_.nNodes; ++i) {
				if (!(parentSetIndices_[i] == other.parentSetIndices_[i]))
					return false;
			}
			return true;
		}
		size_t getHash() const {
			const unsigned int m = 0x5bd1e995;
			const int r = 24;

			unsigned int h = 0xc70f6907;

			for (int i = 0; i < psMap_.nNodes * (sizeof(typename SubsetMap<T>::SubsetIndex) / 4); ++i) {
				unsigned int k = ((unsigned int *)(parentSetIndices_.get()))[i];
				k *= m; 
				k ^= k >> r; 
				k *= m; 
				h *= m; 
				h ^= k;
			}

			h ^= h >> 13;
			h *= m;
			h ^= h >> 15;

			return h;
		}

		struct hash {
			size_t operator()(const ParentSets& parentSets) const {
				return parentSets.getHash();
			}
		};
	};
};

/*namespace std {
	template <typename T>
	struct hash<ParentSetMap<T>::ParentSets> {
		size_t operator()(const ParentSetMap<T>::ParentSets& parentSets) const {
			return parentSets.getHash();
		}
	};
}*/




/**
 * Writes local scores to file.
 */
template <typename T>
void writeScores(std::ostream& file, const ParentSetMap<T>* scores) {
	file << scores->nNodes << std::endl;
	file << scores->maxParents << std::endl;
	file.precision(16);
	for (int node = 0; node < scores->nNodes; ++node) {
		for (auto score : scores->forNode(node)->forSubsets(SortedArraySubset(scores->maxParents)))
			file << log(score.second) << " ";
		file << std::endl;
	}
}


/**
 * Reads local scores from file.
 */
template <typename T>
std::unique_ptr<ParentSetMap<T>> readScores(std::istream& file) {
	int nNodes, maxParents;
	file >> nNodes;
	file >> maxParents;
	ParentSetMap<T>* scores = new ParentSetMap<T>(nNodes, maxParents);
	for (int node = 0; node < nNodes; ++node) {
		for (auto score : scores->forNode(node)->forSubsets(SortedArraySubset(scores->maxParents))) {
			double logScore;
			file >> logScore;
			if (file.fail())
				throw Exception("File corrupted; could not read all scores.");
			setLog(score.second, logScore);
		}
	}
	file >> std::ws;
	if (!file.eof())
		throw Exception("File corrupted; contains more data than expected.");
	return std::unique_ptr<ParentSetMap<T>>(scores);
}




#endif



