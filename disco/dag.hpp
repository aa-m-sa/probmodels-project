/*
 *  BEANDisco: DAG family definition
 *  
 *  Copyright 2015 Teppo Niinimäki <teppo.niinimaki(at)helsinki.fi>
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

//#include "stacksubset.hpp"
#include "common.hpp"

#ifndef DAG_HPP
#define DAG_HPP



/**
 * A structure that describes an arc (directed edge).
 */
struct Arc {
	int tail;
	int head;
	
	Arc() {
		tail = 0;
		head = 0;
	}
	
	Arc(int t, int h) {
		tail = t;
		head = h;
	}
	
	bool operator==(const Arc& other) {
		return tail == other.tail && head == other.head;
	}
	bool operator!=(const Arc& other) {
		return tail != other.tail || head != other.head;
	}
};

std::ostream& operator<<(std::ostream& os, const Arc& arc) {
	return os << arc.tail << " -> " << arc.head;
}


/**
 * A range of all possible arcs between n nodes.
 */
class Arcs {
private:
	int n_;
public:
	Arcs(int n) {
		n_ = n;
	}
	
	class Iterator {
	private:
		const Arcs& arcs_;
		int round_;
		Arc arc_;
	public:
		Iterator(const Arcs& arcs, int round, const Arc& arc) :
			arcs_(arcs),
			round_(round),
			arc_(arc)
		{
		}
		
		Iterator& operator++() {
			++arc_.tail;
			if (arc_.head == arc_.tail)
				++arc_.tail;
			if (arc_.tail >= arcs_.n_) {
				arc_.tail = 0;
				++arc_.head;
				if (arc_.head >= arcs_.n_) {
					arc_.head = 0;
					arc_.tail = 1;
					++round_;
				}
			}
			return *this;
		}
		
		bool operator!=(const Iterator& other) {
			return round_ != other.round_ || arc_ != other.arc_;
		}
		
		const Arc& operator*() {
			return arc_;
		}
	};
	
	Iterator begin() {
		return Iterator(*this, 0, Arc(1, 0));
	}
	
	Iterator end() {
		return Iterator(*this, 1, Arc(1, 0));
	}
};

/**
 * Returns a range of all possible arcs between n nodes.
 */
Arcs allArcs(int n) {
	return Arcs(n);
}



class DagFamily {
private:
	using Family = DagFamily;
	
public:
	const int n;

	DagFamily(int _n)
		: n(_n)
	{}
	
	~DagFamily() {
	}

	bool operator==(const DagFamily& other) const {
		return n == other.n;
	}
	
	double logSize() const {
		static double sl[] = {
			0.0000000000000000, 0.0000000000000000, 1.0986122886681097,
			3.2188758248682007, 6.2971093199339354, 10.284694120497512,
			15.145632107613884, 20.853222705439836, 27.387295103238914,
			34.732237719157354, 42.875669738844121, 51.807553542046141,
			61.519594202603898, 72.004824233066829, 83.257309258224103,
			95.271934287254020, 108.04424500329491, 121.57032751723480,
			135.84671562420959, 150.87031813259847, 166.63836111580798,
			183.14834144370819, 200.39798896779145, 218.38523543500190,
			237.10818869718999};
		
		if (n < sizeof(sl) / sizeof(sl[0]))
			return sl[n];
		else
			return lgamma(n + 1) + log(2.0) * 0.5 * (n * (n - 1)) -
					(-0.55449476947084337 + 0.39748572103907781 * n);
	}

	class Instance {
	public:
		const Family& family;
	private:
		//SortedArraySubset* parentSets_;
		int* parentCounts_;
		int* parentSets_;

		bool isDag() {
			std::vector<int> status(family.n, 0);
			std::function<bool(int)> visit = [&] (int node) {
				if (status[node] == 1) {
					return false;
				}
				if (status[node] == 0) {
					status[node] = 1;
					for (int j = 0; j < parentCounts_[node]; ++j)
						if (!visit(parentSets_[node * family.n + j]))
							return false;
					status[node] = 2;
				}
				return true;
			};
			for (int node = 0; node < family.n; ++node) {
				if (status[node] == 0)
					if (!visit(node))
						return false;
			}
			return true;
		}

		friend struct std::hash<Instance>;
	public:
		Instance(const Family&& _family) = delete;
		Instance(const Family& _family) :
			family(_family)
		{
			//parentSets_ = new SortedArraySubset[](family.n, family.n - 1);
			parentCounts_ = new int[family.n];
			parentSets_ = new int[family.n * family.n];
			for (int i = 0; i < family.n; ++i)
				parentCounts_[i] = 0;
			assert(isDag());
		}
		
		~Instance() {
			delete[] parentCounts_;
			delete[] parentSets_;
		}
		
		Instance(const Instance& other) :
			family(other.family)
		{
			parentCounts_ = new int[family.n];
			parentSets_ = new int[family.n * family.n];
			for (int i = 0; i < family.n; ++i) {
				parentCounts_[i] = other.parentCounts_[i];
				for (int j = 0; j < parentCounts_[i]; ++j)
					parentSets_[i * family.n + j] = other.parentSets_[i * family.n + j];
			}
			assert(isDag());
		}
		
		Instance(Instance&& other) :
			family(other.family),
			parentCounts_(std::move(other.parentCounts_)),
			parentSets_(std::move(other.parentSets_))
		{
			assert(isDag());
		}
		
		Instance& operator=(const Instance& other) {
			assert(family == other.family);
			for (int i = 0; i < family.n; ++i) {
				parentCounts_[i] = other.parentCounts_[i];
				for (int j = 0; j < parentCounts_[i]; ++j)
					parentSets_[i * family.n + j] = other.parentSets_[i * family.n + j];
			}
			assert(isDag());
			return *this;
		}
		
		Instance& operator=(Instance&& other) {
			assert(family == other.family);
			delete[] parentSets_;
			delete[] parentCounts_;
			parentSets_ = std::move(other.parentSets_);
			parentCounts_ = std::move(other.parentCounts_);
			assert(isDag());
			return *this;
		}

		bool operator==(const Instance& other) const {
			if (!(family == other.family))
				return false;
			for (int i = 0; i < family.n; ++i) {
				if (parentCounts_[i] != other.parentCounts_[i])
					return false;
				for (int j = 0; j < parentCounts_[i]; ++j)
					if (parentSets_[i * family.n + j] != other.parentSets_[i * family.n + j])
						return false;
			}
			return true;
		}

		void setParents(int node, const SortedArraySubset& parents) {
			assert(0 <= node && node < family.n);
			parentCounts_[node] = parents.size();
			for (int i = 0; i < parents.size(); ++i)
				parentSets_[node * family.n + i] = parents[i];
			assert(isDag());
		}

		SortedSubsetRange getParents(int node) const {
			assert(0 <= node && node < family.n);
			return SortedSubsetRange(parentCounts_[node], parentSets_ + node * family.n);
		}

		/*SortedArraySubset getParents(int node) const {
			assert(0 <= node && node < family.n);
			SortedArraySubset parents(parentCounts_[node]);
			for (int i = 0; i < parentCounts_[node]; ++i)
				parents.insertLargest(parentSets_[node * family.n + i]);
				//parents.insert(parentSets_[node * family.n + i]);
			return parents;
		}*/

		bool contains(const Arc& arc) const {
			// TODO: O(n) -> O(1) ?????
			assert(0 <= arc.head && arc.head < family.n);
			assert(0 <= arc.tail && arc.tail < family.n);
			for (int i = 0; i < parentCounts_[arc.head]; ++i)
				if (parentSets_[arc.head * family.n + i] == arc.tail)
					return true;
			return false;
		}

		void getAdjMat(SquareMatrix<bool>& adjMat) const {
			assert(family.n == adjMat.getN());
			adjMat.setAll(false);
			for (int node = 0; node < family.n; ++node)
				for (int j = 0; j < parentCounts_[node]; ++j)
					adjMat(parentSets_[node * family.n + j], node) = true;
		}

		//void getParents(int i, StackSubset& parentset) const {
		//	parentset.clear();
		//	for (int j = 0; j < family.n; ++j)
		//		if (adjMat_(j, i))
		//			parentset.push(j);
		//}

		void writeTo(std::ostream& os) const {
			for (int i = 0; i < family.n; ++i) {
				if (i > 0)
					os << ", ";
				os << i << " <- {";
				for (int j = 0; j < parentCounts_[i]; ++j) {
					if (j > 0)
						os << ", ";
					os << parentSets_[i * family.n + j];
				}
				os << "}";
			}
		}
	
		//void readFrom(std::istream& is) {
		//	std::string tmp;
		//	for (int b = 0; b < family.nBuckets(); ++b) {
		//		if (b > 0) {
		//			is >> tmp;
		//			assert(tmp == "≺");
		//		}
		//		for (int i = 0; i < family.bucketSize(b); ++i) {
		//			is >> order_[b * family.maxBucketSize + i];
		//		}
		//	assert(isDag());
		//	}
		//}
	};/**/


	/*class RandMadiganYork : public MCProposalDist<Family> {
	private:
		const Family& family_;
	public:
		RandMadiganYork(const Family& family) :
			family_(family)
		{}

		void randStep(Instance& instance) const {
			switch(randuint(3)) {
				case 0:
					break;
				case 1:
					break;
				case 2:
					break;
				default:
					assert(0);
			}
			int i = randuint(family_.n);
			int j = (i + randuint(family_.n - 1)) % family_.n;
			instance.swapNodesAt(i, j);
		}
	};*/
	
/*	class Instance {
	private:
		ArcMap<bool> adjMat_;
		
	public:
		const Family& family;
		
		Instance(const Family& _family) :
			family(_family),
			adjMat_(family.n)
		{
		}
		
		~Instance() {
		}
		
		Instance(const Instance& other) :
			family(other.family),
			adjMat_(other.adjMat_)
		{
		}
		
		Instance(Instance&& other) :
			family(other.family),
			adjMat_(std::move(other.adjMat_))
		{
		}
		
		Instance& operator=(const Instance& other) {
			assert(family == other.family);
			adjMat_ = other.adjMat_;
			return *this;
		}
		
		Instance& operator=(Instance&& other) {
			assert(family == other.family);
			adjMat_ = std::move(other.adjMat_);
			return *this;
		}

		void setParents(int node, const StackSubset& parents) {
			for (int j = 0; j < family.n; ++j)
				adjMat_[Arc(j, node)] = false;
			for (int i = 0; i < parents.size(); ++i)
				adjMat_[Arc(parents[i], node)] = true;
		}
		
		//void swap(size_t i1, size_t i2) {
		//	std::swap(order_[i1], order_[i2]);
		//	std::swap(invOrder_[order_[i1]], invOrder_[order_[i2]]);
		//}

		//void rand() {
		//	for (int i = 0; i < family.n; ++i) {
		//		int j = i + randuint(family.n - i);
		//		swap(i, j);
		//	}
		//}
		
		const ArcMap<bool>& getAdjMat() {
			return adjMat_;
		}

		//void getParents(int i, StackSubset& parentset) const {
		//	parentset.clear();
		//	for (int j = 0; j < family.n; ++j)
		//		if (adjMat_(j, i))
		//			parentset.push(j);
		//}
	};*/
};


std::ostream& operator<<(std::ostream& os, const DagFamily::Instance& instance) {
	instance.writeTo(os);
	return os;
}

/*std::istream& operator>>(std::istream& is, DagFamily::Instance& instance) {
	instance.readFrom(is);
	return is;
}*/

namespace std {

	template <>
	struct hash<DagFamily::Instance> {
		size_t operator()(const DagFamily::Instance& dag) const {
			const unsigned int m = 0x5bd1e995;
			const int r = 24;

			unsigned int h = 0xc70f6907;

			for (int i = 0; i < dag.family.n; ++i) {
				unsigned int k = ((unsigned int *)(dag.parentCounts_))[i];
				k *= m; 
				k ^= k >> r; 
				k *= m; 
				h *= m; 
				h ^= k;

				for (int j = 0; j < dag.parentCounts_[i]; ++j) {
					unsigned int k = ((unsigned int *)(dag.parentSets_))[i * dag.family.n + j];
					k *= m; 
					k ^= k >> r; 
					k *= m; 
					h *= m; 
					h ^= k;
				}
			}

			h ^= h >> 13;
			h *= m;
			h ^= h >> 15;

			return h;
		}
	};

}


#endif
