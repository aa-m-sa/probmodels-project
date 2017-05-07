/*
 *  BEANDisco: bucket order definition
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

//#include "stacksubset.hpp"
#include "common.hpp"

#ifndef BUCKETORDER_HPP
#define BUCKETORDER_HPP


class BucketOrderFamily {
private:
	using Family = BucketOrderFamily;
	int* predSetBucketMasks_;

	BucketOrderFamily(const BucketOrderFamily&) = delete; // disable copy constructor
	BucketOrderFamily& operator=(const BucketOrderFamily&) = delete; // disable assignment
public:
	const int n;
	const int maxBucketSize;
	BucketOrderFamily(int _n, int _bucketSize) : n(_n), maxBucketSize(_bucketSize) {
		assert(1 <= maxBucketSize && maxBucketSize < 32);
		assert(1 <= n);
		
		predSetBucketMasks_ = new int[n * n];
		for (int i = 0; i < n; ++i) {
			int b = getBucket(i);
			for (int j = 0; j < n; ++j) {
				if (i == j || b != getBucket(j))
					predSetBucketMasks_[i * n + j] = 0;
				else {
					int ii = getIndexInBucket(i);
					int ji = getIndexInBucket(j);
					if (ji == bucketSize(b) - 1)
						predSetBucketMasks_[i * n + j] = 1 << ii;
					else
						predSetBucketMasks_[i * n + j] = 1 << ji;
				}
			}
		}
	}
	
	~BucketOrderFamily() {
		delete[] predSetBucketMasks_;
	}
	
	int bucketSize(int b) const {
		assert(0 <= b && b < nBuckets());
		return min(maxBucketSize, n - b * maxBucketSize);
	}
	
	int nBuckets() const {
		return (n + maxBucketSize - 1) / maxBucketSize;
	}
	
	int nDownSets() const {
		int nb = nBuckets();
		return (nb - 1) * (1 << maxBucketSize) + (1 << bucketSize(nb - 1)) - nb + 1;
	}
	
	int getBucket(int v) const {
		return v / maxBucketSize;
	}
	
	int getIndexInBucket(int v) const {
		return v % maxBucketSize;
	}
	
	bool operator==(const Family& bof) const {
		return n == bof.n && maxBucketSize == bof.maxBucketSize;
	}
	
	double logSize() const {
		double s = lgamma(n + 1);
		for (int b = 0; b < nBuckets(); ++b)
			s -= lgamma(bucketSize(b) + 1);
		return s;
	}
	
	double logNumLinExts() const {
		double s = 0.0;
		for (int b = 0; b < nBuckets(); ++b)
			s += lgamma(bucketSize(b) + 1);
		return s;
	}
	
	class Instance;
	class PredSetRange;
	class PredSet;
	class DownSet;
	template <typename T> class PredSetMap;
	template <typename T> class DownSetMap;
	
	class Instance {
	public:
		const Family& family;
	private:
		int* order_;
		int* invOrder_;
		//std::unique_ptr<int[]> order_;
		//std::unique_ptr<int[]> invOrder_;
		friend class BucketOrderFamily;
	public:
		
		Instance(const Family&& _family) = delete;
		Instance(const Family& _family) :
			family(_family)
		{
			order_ = new int[family.n];
			invOrder_ = new int[family.n];
			for (int i = 0; i < family.n; ++i) {
				order_[i] = i;
				invOrder_[i] = i;
			}
		}
		~Instance() {
			delete[] order_;
			delete[] invOrder_;
		}
		
		Instance(const Instance& other) :
			family(other.family)
		{
			order_ = new int[family.n];
			invOrder_ = new int[family.n];
			for (int i = 0; i < family.n; ++i) {
				order_[i] = other.order_[i];
				invOrder_[i] = other.invOrder_[i];
			}
		}
		
		Instance(Instance&& other) :
			family(other.family),
			order_(other.order_),
			invOrder_(other.invOrder_)
		{
			other.order_ = nullptr;
			other.invOrder_ = nullptr;
		}

		Instance& operator=(const Instance& other) {
			assert(family == other.family);
			for (int i = 0; i < family.n; ++i) {
				order_[i] = other.order_[i];
				invOrder_[i] = other.invOrder_[i];
			}
			return *this;
		}

		Instance& operator=(Instance&& other) {
			assert(family == other.family);
			order_ = other.order_;
			invOrder_ = other.invOrder_;
			other.order_ = nullptr;
			other.invOrder_ = nullptr;
			return *this;
		}

		/*Instance(const Family&& _family) = delete;
		Instance(const Family& _family) :
			family(_family),
			order_(new int[family.n]),
			invOrder_(new int[family.n])
		{
			for (int i = 0; i < family.n; ++i) {
				order_[i] = i;
				invOrder_[i] = i;
			}
		}
		
		~Instance() {
		}
		
		Instance(const Instance& other) :
			family(other.family),
			order_(new int[family.n]),
			invOrder_(new int[family.n])
		{
			for (int i = 0; i < family.n; ++i) {
				order_[i] = other.order_[i];
				invOrder_[i] = other.invOrder_[i];
			}
		}
		
		Instance(Instance&& other) :
			family(other.family),
			order_(std::move(other.order_)),
			invOrder_(std::move(other.invOrder_))
		{
		}
		
		Instance& operator=(const Instance& other) {
			assert(family == other.family);
			for (int i = 0; i < family.n; ++i) {
				order_[i] = other.order_[i];
				invOrder_[i] = other.invOrder_[i];
			}
			return *this;
		}
		
		Instance& operator=(Instance&& other) {
			assert(family == other.family);
			order_ = std::move(other.order_);
			invOrder_ = std::move(other.invOrder_);
			return *this;
		}*/
		
		//int operator[](size_t i) const {
		//	assert(0 <= i && i < family.n);
		//	return order_[i];
		//}
		
		int* tail(size_t b) {
			return order_ + b * family.maxBucketSize;
		}
		
		int tailLength(size_t b) const {
			return family.n - b * family.maxBucketSize;
		}
		
		const int* getOrder() const {
			return order_;
		}
		
		int getIndex(int j) const {
			assert(0 <= j && j < family.n);
			return invOrder_[j];
		}
		
		template <typename Subset>
		void getPotPreds(int elem, Subset& potPreds) const {
			int i = getIndex(elem);
			int b = family.getBucket(i);
			for (int j = 0; j < b * family.maxBucketSize + family.bucketSize(b); ++j)
				if (j != i)
					potPreds.insert(order_[j]);
					//potPreds.push(order_[j]);
		}
		
		/*void swap(size_t i1, size_t i2) {
			std::swap(order_[i1], order_[i2]);
			std::swap(invOrder_[order_[i1]], invOrder_[order_[i2]]);
		}*/
		void swapElemsAt(size_t i1, size_t i2) {
			std::swap(order_[i1], order_[i2]);
			std::swap(invOrder_[order_[i1]], invOrder_[order_[i2]]);
		}
		
		void rand() {
			for (int i = 0; i < family.n; ++i) {
				int j = i + randuint(family.n - i);
				swapElemsAt(i, j);
			}
		}
		
		/*template <class ProbComputer, typename Real>
		bool mcmcPOStep(ProbComputer& pc, Real& p) {
			int v1 = randuint(family.n);
			int b1 = v1 / family.maxBucketSize;
			//int bs1 = min(family.maxBucketSize, family.n - b1 * family.maxBucketSize);
			int bs1 = family.bucketSize(b1);
			//int v2 = ((b1 + 1) * family.maxBucketSize + randuint(family.n - bs1)) % family.n;
			int v2 = (b1 * family.maxBucketSize + bs1 + randuint(family.n - bs1)) % family.n;
			swap(v1, v2);
			Real pnew = pc.calcProb(family, *this);
			//Real pnew = pc.updateMargin(*this, v1, v2);
			if (randu() < to<double>(pnew / p)) {
				p = pnew;
				return true;
			} else {
				swap(v1, v2);
				//pc.updateMargin(*this, v1, v2);
				return false;
			}
		}*/
		
		/*void randSwap() {
			int v1 = randuint(family.n);
			int b1 = v1 / family.maxBucketSize;
			int bs1 = family.bucketSize(b1);
			int v2 = (b1 * family.maxBucketSize + bs1 + randuint(family.n - bs1)) % family.n;
			swap(v1, v2);
		}*/
		
		/*void mcmcProposeSwap(int& v1, int& v2) {
			v1 = randuint(family.n);
			int b1 = v1 / family.maxBucketSize;
			//int bs1 = min(family.maxBucketSize, family.n - b1 * family.maxBucketSize);
			int bs1 = family.bucketSize(b1);
			//int v2 = ((b1 + 1) * family.maxBucketSize + randuint(family.n - bs1)) % family.n;
			v2 = (b1 * family.maxBucketSize + bs1 + randuint(family.n - bs1)) % family.n;
			swap(v1, v2);
		}*/
		
		PredSetRange predSets(int elem) const {
			return PredSetRange(*this, elem);
		}
		
		void print(bool newline = true) const {
			for (int b = 0; b < family.nBuckets(); ++b) {
				if (b > 0)
					printf(" ≺  ");
				for (int i = 0; i < family.bucketSize(b); ++i) {
					printf("%d ", order_[b * family.maxBucketSize + i]);
				}
			}
			if (newline)
				printf("\n");
		}

		void writeTo(std::ostream& os) const {
			for (int b = 0; b < family.nBuckets(); ++b) {
				if (b > 0)
					os << " ≺ ";
				for (int i = 0; i < family.bucketSize(b); ++i) {
					if (i > 0)
						os << " ";
					os << order_[b * family.maxBucketSize + i];
				}
			}
		}

		/*void writeTo(std::ostream& os) const {
			for (int b = 0; b < family.nBuckets(); ++b) {
				if (b > 0)
					os << " ≺ ";
				os << "{"
				for (int i = 0; i < family.bucketSize(b); ++i) {
					if (i > 0)
						os << ", ";
					os << order_[b * family.maxBucketSize + i];
				}
				os << "}";
			}
		}*/
	
		void readFrom(std::istream& is) {
			std::string tmp;
			for (int b = 0; b < family.nBuckets(); ++b) {
				if (b > 0) {
					is >> tmp;
					//assert(tmp == "≺");
					if (tmp != "≺")
						throw Exception("Reading a bucket order failed: expected '≺' but encountered '%s'", tmp);
				}
				for (int i = 0; i < family.bucketSize(b); ++i) {
					int ii = b * family.maxBucketSize + i;
					is >> order_[ii];
					invOrder_[order_[ii]] = ii;
				}
			}
		}

		/*void readFrom(std::istream& is) {
			std::string tmp;
			for (int b = 0; b < family.nBuckets(); ++b) {
				if (b > 0) {
					is >> tmp;
					//assert(tmp == "≺");
					if (tmp != "≺")
						throw Exception("Reading a bucket order failed: expected '≺' but encountered '%s'", tmp);
				}
				for (int i = 0; i < family.bucketSize(b); ++i) {
					is >> order_[b * family.maxBucketSize + i];
				}
			}
		}*/
	};

	class UniformRandSwap : public MCProposalDist<Family> {
	private:
		const Family& family_;
		std::vector<std::pair<int,int>> nodePairs_;
	public:
		UniformRandSwap(const Family&& family) = delete;
		UniformRandSwap(const Family& family) :
			family_(family)
		{
			for (int v1 = 0; v1 < family_.maxBucketSize * (family_.nBuckets() - 1); ++v1) {
				int b1 = v1 / family_.maxBucketSize;
				for (int v2 = (b1 + 1) * family_.maxBucketSize; v2 < family_.n; ++v2) {
					nodePairs_.push_back(std::pair<int,int>(v1, v2));
				}
			}
		}

		void randStep(Instance& instance) const {
			auto nodePair = nodePairs_[randuint(nodePairs_.size())];
			instance.swapElemsAt(nodePair.first, nodePair.second);
		}
	};

	/*class UniformRandSwap : public MCProposalDist<Family> {
	private:
		const Family& family_;
	public:
		UniformRandSwap(const Family&& family) = delete;
		UniformRandSwap(const Family& family) :
			family_(family)
		{}

		void randStep(Instance& instance) const {
			int v1 = randuint(family_.n);
			int b1 = v1 / family_.maxBucketSize;
			int bs1 = family_.bucketSize(b1);
			int v2 = (b1 * family_.maxBucketSize + bs1 + randuint(family_.n - bs1)) % family_.n;
			instance.swapElemsAt(v1, v2);
		}
	};*/
	
	// TODO: should also update invOrder_  ->  use instance_.swap(?,?)
	class Enumerator {
	private:
		const Family& family_;
		int** index_;
		Instance instance_;
		int b_;		// current bucket
		int i_;		// current element in current bucket
	public:
		Enumerator(const Family& family) : family_(family), instance_(family_) {
			// allocate index
			index_ = new int*[family_.nBuckets() - 1];
			for (int b = 0; b < family_.nBuckets() - 1; ++b)
				index_[b] = new int[family_.bucketSize(b) + 1];
			
			// init first state
			init();
		}
		~Enumerator() {
			for (int b = 0; b < family_.nBuckets() - 1; ++b)
				delete[] index_[b];
			delete[] index_;
		}
		
		void init() {
			b_ = 0;
			i_ = 0;
			for (; b_ < family_.nBuckets() - 1; ++b_) {
				for (; i_ <= family_.bucketSize(b_); ++i_)
					index_[b_][i_] = -1;
				i_ = 0;
			}
			//printf("  b_=%d   i_=%d\n", b_, i_);
		}
		
		bool next() {
			// find the next elements to change
			while (true) {
				// no elements left in this bucket?
				if (i_ == 0) {
					// if no more buckets left => stop
					if (b_ == 0)
						return false;
					// move to the previous bucket
					--b_;
					i_ = family_.bucketSize(b_);
				}
				// undo the previous swap, if there was such
				if (index_[b_][i_] > index_[b_][i_-1]) {
					std::swap(instance_.tail(b_)[i_-1], instance_.tail(b_+1)[index_[b_][i_]]);
					//printf("  b_=%d   i_=%d  index_[b_][i_]=%d\n", b_, i_, index_[b_][i_]);
					//printf("  "); instance_.print(false); printf("\n");
				}
				// element with unused values left (possible swaps left)? => go on
				if (++index_[b_][i_] < instance_.tailLength(b_+1))
					break;
				// otherwise try previous element
				--i_;
			}
			
			// new swap
			std::swap(instance_.tail(b_)[i_-1], instance_.tail(b_+1)[index_[b_][i_]]);
			//printf("  b_=%d   i_=%d  index_[b_][i_]=%d\n", b_, i_, index_[b_][i_]);
			//printf("  "); instance_.print(false); printf("\n");
			
			// init the tail and move back to the last bucket / element
			for (; b_ < family_.nBuckets() - 1; ++b_) {
				for (++i_; i_ <= family_.bucketSize(b_); ++i_)
					index_[b_][i_] = index_[b_][i_ - 1];
				i_ = 0;
			}
			//printf("  b_=%d   i_=%d\n", b_, i_);
			return true;
		}
		
		const Instance& get() const {
			return instance_;
		}
	};
	
	
	/**
	 * Representation for an element and its predecessor set.
	 */
	class PredSet {
	private:
	//public:
		const Instance& instance_;
		int i_;
		unsigned int setMask_;
	public:
		PredSet(const Instance&& instance, int base) = delete;
		PredSet(const Instance& instance, int base) :
			instance_(instance)
		{
			i_ = instance_.getIndex(base);
			setMask_ = 0;
		}

		PredSet(int elem, const DownSet& downSet) :
			instance_(downSet.instance_)
		{
			i_ = instance_.getIndex(elem);
			int b = instance_.family.getBucket(i_);
			assert(downSet.bucket_ == b);
			int bs = instance_.family.bucketSize(b);
			int j = instance_.family.getIndexInBucket(i_);
			//int vMask = 1 << j;
			int lastMask = 1 << (bs - 1);
			//assert(!(downSet.setMask_ & vMask));
			assert(!(downSet.setMask_ & (1 << j)));
			//setMask_ = downSet.setMask_ | ((downSet.setMask_ >> (bs - 1 - j)) & vMask);
			setMask_ = (downSet.setMask_ & (~lastMask)) | ((downSet.setMask_ & lastMask) >> (bs - 1 - j));
		}
		
		void insert(int elem) {
			int j = instance_.getIndex(elem);
			int elemMask = instance_.family.predSetBucketMasks_[i_ * instance_.family.n + j];
			assert(!elemMask || !(setMask_ & elemMask));
			setMask_ ^= elemMask;
		}
		
		void remove(int elem) {
			int j = instance_.getIndex(elem);
			int elemMask = instance_.family.predSetBucketMasks_[i_ * instance_.family.n + j];
			assert(!elemMask || (setMask_ & elemMask));
			setMask_ ^= elemMask;
		}

		void insertLargest(int elem) {
			insert(elem);
		}

		void removeLargest(int elem) {
			remove(elem);
		}
		
		void getElements(SortedArraySubset& elements) const {
			DownSet(*this).getElements(elements);
		}

		SortedArraySubset getElements() const {
			SortedArraySubset elements(instance_.family.n);
			getElements(elements);
			return elements;
		}

		/*void setElements(StackSubset& ss) const {
			ss.clear();
			int b = instance_.family.getBucket(i_);
			for (int i = 0; i < b * instance_.family.maxBucketSize; ++i)
				ss.push(instance_[i]);
			for (int j = 0; j < instance_.family.bucketSize(b); ++j)
				if (setMask_ & instance_.family.predSetBucketMasks_[
						i_ * instance_.family.n + j])
					ss.push(instance_[j]);
		}*/
		
		/*void setSuperOf(const StackSubset& set) {
			setMask_ = 0;
			int bucket = instance_.family.getBucket(i_);
			for (int v = 0; v < set.size(); ++v) {
				int j = instance_.getIndex(set[v]);
				int b = instance_.family.getBucket(j);
				//int j = instance_.family.getIndexInBucket(i);
				assert(b <= bucket);
				//if (b == bucket) {
				setMask_ |= instance_.family.predSetBucketMasks_[
						i_ * instance_.family.n + j];
				//}
			}
		}*/
		
		friend class PredSetRange;
		template <typename T> friend class PredSetMap;

		friend class DownSet;
	};

	/**
	 * A range/collection of predecessor sets of given element.
	 */
	class PredSetRange {
	private:
		const Instance& instance_;
		int base_;
	public:
		PredSetRange(const Instance&& instance, int base) = delete;
		PredSetRange(const Instance& instance, int base) :
			instance_(instance),
			base_(base)
		{
		}
		
		class Iterator {
		private:
			PredSet predSet_;
		public:
			Iterator(const Instance& instance, int base, unsigned int setMask) :
				predSet_(instance, base)
			{
				predSet_.setMask_ = setMask;
			}
			
			void operator++() {
				++predSet_.setMask_;
			}
			
			PredSet& operator*() {
				return predSet_;
			}
			
			bool operator!=(const Iterator& other) {
				return predSet_.setMask_ != other.predSet_.setMask_;
			}
		};
		
		Iterator begin() {
			return Iterator(instance_, base_, 0);
		}
		
		Iterator end() {
			int b = instance_.family.getBucket(instance_.getIndex(base_));
			int bs = instance_.family.bucketSize(b);
			unsigned int setMask = 1 << (bs - 1);
			return Iterator(instance_, base_, setMask);
		}
	};
	
	
	/**
	 * Representation of downset (lower set, ideal).
	 */
	class DownSet {
	private:
		const Instance& instance_;
		int bucket_;
		int setMask_;
		
	public:
		DownSet(const Instance&& instance) = delete;
		DownSet(const Instance& instance) :
			instance_(instance)
		{
			bucket_ = 0;
			setMask_ = 0;
		}
		
		DownSet(const PredSet& predSet) :
			instance_(predSet.instance_)
		{
			bucket_ = instance_.family.getBucket(predSet.i_);
			int bs = instance_.family.bucketSize(bucket_);
			int i = instance_.family.getIndexInBucket(predSet.i_);
			int vMask = 1 << i;
			setMask_ = (predSet.setMask_ & (~vMask)) |
					((predSet.setMask_ & vMask) << (bs - 1 - i));
		}
		
		void clear() {
			bucket_ = 0;
			setMask_ = 0;
		}

		void insertWithinBucket(int elem) {
			int i = instance_.getIndex(elem);
			assert(instance_.family.getBucket(i) == bucket_);
			int vMask = 1 << instance_.family.getIndexInBucket(i);
			assert(!(setMask_ & vMask));
			setMask_ ^= vMask;
		}

		void insert(int elem) {
			int i = instance_.getIndex(elem);
			assert(instance_.family.getBucket(i) == bucket_);
			int vMask = 1 << instance_.family.getIndexInBucket(i);
			assert(!(setMask_ & vMask));
			setMask_ ^= vMask;
			if (setMask_ == (1 << instance_.family.bucketSize(bucket_)) - 1) {
				setMask_ = 0;
				++bucket_;
			}
		}
		
		void remove(int elem) {
			int i = instance_.getIndex(elem);
			if (setMask_ == 0 ) {
				--bucket_;
				setMask_ = (1 << instance_.family.bucketSize(bucket_)) - 1;
			}
			assert(instance_.family.getBucket(i) == bucket_);
			int vMask = 1 << instance_.family.getIndexInBucket(i);
			assert(setMask_ & vMask);
			setMask_ ^= vMask;
		}

		/*void setSuperOf(const StackSubset& set) {
			clear();
			for (int v = 0; v < set.size(); ++v) {
				int b = set[v] / instance_.family.maxBucketSize;
				int i = set[v] % instance_.family.maxBucketSize;
				if (b > bucket_) {
					bucket_ = b;
					setMask_ = (1 << i);
				} else if (b == bucket_) {
					setMask_ |= (1 << i);
				}
			}
		}*/

		void getElements(SortedArraySubset& elements) const {
			std::vector<int> elems;
			elems.reserve(instance_.family.n);
			int i = 0;
			for (; i < bucket_ * instance_.family.maxBucketSize; ++i)
				elems.push_back(instance_.order_[i]);
			int remainingMask = setMask_;
			while (remainingMask) {
				if (remainingMask & 1)
					elems.push_back(instance_.order_[i]);
				remainingMask >>= 1;
				++i;
			}
			std::sort(elems.begin(), elems.end());
			for (auto elem : elems)
				elements.insertLargest(elem);
		}

		SortedArraySubset getElements() const {
			SortedArraySubset elements(instance_.family.n);
			getElements(elements);
			return elements;
		}
		
		template <typename T> friend class DownSetMap;
		friend class PredSet;
	};
	
	
	/**
	 *  A map for all possible predecessor sets of all elements.
	 */
	template <typename T>
	class PredSetMap {
	private:
		const Family& family_;
		T** data_;
		
		int maxPredSetSize(int i) const {
			int b = family_.getBucket(i);
			return family_.bucketSize(b) - 1;
		}
		
		size_t predSetCount(int i) const {
			return 1 << maxPredSetSize(i);
		}
		
		void init() {
			//printf("init(): family_.n = %d\n", family_.n);
			data_ = new T*[family_.n];
			for (int i = 0; i < family_.n; ++i) {
				data_[i] = new T[predSetCount(i)];
			}
		}

		void copy(const PredSetMap& other) {
			assert(family_ == other.family_);
			for (int i = 0; i < family_.n; ++i) {
				memcpy(data_[i], other.data_[i], predSetCount(i) * sizeof(T));
			}
		}
		
		T operator()(int i, int setMask) const {
			return data_[i][setMask];
		}

		template <typename U> friend class DownSetMap;
	public:
		PredSetMap(const Family&& family) = delete;
		PredSetMap(const Family& family) : family_(family) {
			init();
		}

		PredSetMap(const PredSetMap& other) : family_(other.family_) {
			init();
			copy(other);
		}
		
		~PredSetMap() {	
			//printf("~PredSetMap(): family_.n = %d\n", family_.n);
			for (int i = 0; i < family_.n; ++i)
				delete[] data_[i];
			delete[] data_;
		}
		
		void operator=(const PredSetMap& other) {
			copy(other);
		}
		
		void setAll(T value) {
			for (int i = 0; i < family_.n; ++i)
				for (int ssMask = 0; ssMask < predSetCount(i); ++ssMask)
					data_[i][ssMask] = value;
		}
		
		T& operator[](const PredSet& ps) {
			return data_[ps.i_][ps.setMask_];
		}

		T operator[](const PredSet& ps) const {
			return data_[ps.i_][ps.setMask_];
		}

		T& operator()(const PredSet& ps) {
			return data_[ps.i_][ps.setMask_];
		}

		T operator()(const PredSet& ps) const {
			return data_[ps.i_][ps.setMask_];
		}

		/*void setAll(T value) {
			for (int i = 0; i < family_.nDownSets(); ++i)
				data_[i] = value;
		}
		
		T getEmpty() {
			return data_[0];
		}
		
		T getFull() {
			return data_[family_.nDownSets() - 1];
		}*/
		
		template <typename Op>
		void zetaTransform(const Instance& instance, int elem) {
			int i = instance.getIndex(elem);
			int psSize = maxPredSetSize(i);
			for (int j = 0; j < psSize; ++j) {
				int vMask = 1 << j;
				for (int psMask = 0; psMask < (1 << psSize); ++psMask) {
					if (!(vMask & psMask))
						Op::add(data_[i][psMask | vMask], data_[i][psMask]);
				}
			}
		}

		template <typename Op>
		void upZetaTransform(const Instance& instance, int elem) {
			int i = instance.getIndex(elem);
			int psSize = maxPredSetSize(i);
			for (int j = 0; j < psSize; ++j) {
				int vMask = 1 << j;
				for (int psMask = 0; psMask < (1 << psSize); ++psMask) {
					if (vMask & psMask)
						Op::add(data_[i][psMask ^ vMask], data_[i][psMask]);
				}
			}
		}

		class ZetaTransformData {
		private:
			T*** data_;
			const PredSetMap& psmap_;
		public:
			ZetaTransformData(const PredSetMap&& psmap) = delete;
			ZetaTransformData(const PredSetMap& psmap)
				: psmap_(psmap)
			{
				data_ = new T**[psmap_.family_.n];
				for (int i = 0; i < psmap_.family_.n; ++i) {
					int psSize = psmap_.maxPredSetSize(i);
					data_[i] = new T*[psSize];
					for (int j = 0; j < psSize; ++j) {
						data_[i][j] = new T[psmap_.predSetCount(i)];
					}
				}
			}
			~ZetaTransformData() {
				for (int i = 0; i < psmap_.family_.n; ++i) {
					int psSize = psmap_.maxPredSetSize(i);
					for (int j = 0; j < psSize; ++j)
						delete[] data_[i][j];
					delete[] data_[i];
				}
				delete[] data_;
			}
			//template <typename Op>
			//friend void zetaTransform(const Instance& instance, int elem, ZetaTransformData* partialData);
			friend class PredSetMap;
		};

		template <typename Op>
		void zetaTransform(const Instance& instance, int elem, ZetaTransformData* partialData) {
			int i = instance.getIndex(elem);
			int psSize = maxPredSetSize(i);
			for (int j = 0; j < psSize; ++j) {
				int vMask = 1 << j;
				for (int psMask = 0; psMask < (1 << psSize); ++psMask) {
					partialData->data_[i][j][psMask] = data_[i][psMask];
					if (vMask & psMask)
						Op::add(data_[i][psMask], data_[i][psMask ^ vMask]);
				}
			}
		}

		template <typename Selector>
		void backTrackZetaTransformSet(const Instance& instance, int elem, const ZetaTransformData* partialData, const PredSet& predSet, PredSet& sourceSet) const {
			//int i = instance.getIndex(elem);
			int i = predSet.i_;
			int psSize = maxPredSetSize(i);
			//int psMask = (1 << psSize) - 1;
			int psMask = predSet.setMask_;
			T curVal = data_[i][psMask];
			Selector selector(2);
			for (int j = psSize - 1; j >= 0; --j) {
				int vMask = 1 << j;
				if (psMask & vMask) {
					assert(curVal == partialData->data_[i][j][psMask ^ vMask] + partialData->data_[i][j][psMask]);
					selector.init(curVal);
					selector.add(0, partialData->data_[i][j][psMask]);
					selector.add(vMask, partialData->data_[i][j][psMask ^ vMask]);
					psMask ^= selector.get();
					curVal = partialData->data_[i][j][psMask];
				}
				else {
					assert(curVal == partialData->data_[i][j][psMask]);
				}
			}
			sourceSet.i_ = i;
			sourceSet.setMask_ = psMask;
		}
		
	};


	/**
	 *  A map for all possible downsets (lower sets, ideals).
	 */
	template <typename T>
	class DownSetMap {
	private:
		const Family& family_;
		T* data_;
		
		void init() {
			data_ = new T[family_.nDownSets()];
		}
		
		void copy(const DownSetMap& other) {
			assert(&family_ == other.family_);
			memcpy(data_, other.data_, family_.nDownSets() * sizeof(T));
		}
		
	public:
		DownSetMap(const Family&& family) = delete;
		DownSetMap(const Family& family) : family_(family) {
			init();
		}

		DownSetMap(const DownSetMap& other) : family_(other.family_) {
			init();
			copy(other);
		}

		~DownSetMap() {
			delete[] data_;
		}
		
		void operator=(const DownSetMap& other) {
			copy(other);
		}
		
		T operator[] (const DownSet& downSet) const {
			return data_[downSet.bucket_ * ((1 << family_.maxBucketSize) - 1) + downSet.setMask_];
		}

		void setAll(T value) {
			for (int i = 0; i < family_.nDownSets(); ++i)
				data_[i] = value;
		}
		
		T getEmpty() const {
			return data_[0];
		}
		
		T getFull() const {
			return data_[family_.nDownSets() - 1];
		}
		
		template <typename Op>
		void forwardSum(const PredSetMap<T>& psmap) {
			data_[0] = Op::one();
			// for each bucket
			for (int b = 0; b < family_.nBuckets(); ++b) {
				int bucketSize = family_.bucketSize(b);
				// index increment from bucket number
				int bi = b * ((1 << family_.maxBucketSize) - 1);
				// enumerate all Ŷ:s in the bucket
				for (int yHatMask = 1; yHatMask < (1 << bucketSize); ++yHatMask) {
					data_[bi + yHatMask] = Op::zero();
					// for each direct subideal
					//int psMask = yHatMask >> 1;
					for (int v = 0; v < bucketSize; ++v) {
						// variable mask
						int vMask = (1 << v);
						if (yHatMask & vMask) {
							int sHatMask = yHatMask ^ vMask;
							int lastV = bucketSize - 1;
							int psMask = (sHatMask & ~(1 << lastV)) |
									((sHatMask >> (lastV - v)) & vMask);
							// variable index
							int i = b * family_.maxBucketSize + v;
							Op::add(data_[bi + yHatMask],
									Op::prod(psmap(i, psMask), data_[bi + sHatMask]));
						}
					}
				}
			}
		}
		
		template <class Selector>
		void backTrackForwardSumOrder(const PredSetMap<T>& psmap, const Instance& instance, std::vector<int>& order) const {
			assert(order.size() == family_.n);
			Selector selector(family_.n);
			int i = family_.n - 1;
			for (int b = family_.nBuckets() - 1; b >= 0; --b) {
				int bucketSize = family_.bucketSize(b);
				int bi = b * ((1 << family_.maxBucketSize) - 1);
				int yHatMask = (1 << bucketSize) - 1;
				while (yHatMask) {
					selector.init(data_[bi + yHatMask]);
					for (int v = 0; v < bucketSize; ++v) {
						int vMask = (1 << v);
						if (yHatMask & vMask) {
							int sHatMask = yHatMask ^ vMask;
							int lastV = bucketSize - 1;
							int psMask = (sHatMask & ~(1 << lastV)) |
									((sHatMask >> (lastV - v)) & vMask);
							// variable index
							int j = b * family_.maxBucketSize + v;
							selector.add(v, psmap(j, psMask) * data_[bi + sHatMask]);
						}
					}
					int v = selector.get();
					int j = b * family_.maxBucketSize + v;
					order[i] = instance.order_[j];
					//predSets[i] = selector.get();
					int vMask = (1 << v);
					yHatMask ^= vMask;
					--i;
				}
			}
		}
		
		template <typename Op>
		void backwardSum(const PredSetMap<T>& psmap) {
			data_[family_.nDownSets() - 1] = Op::one();
			// for each bucket
			for (int b = family_.nBuckets() - 1; b >= 0; --b) {
				int bucketSize = family_.bucketSize(b);
				// index increment from bucket number
				int bi = b * ((1 << family_.maxBucketSize) - 1);
				// enumerate all Y̌:s in the bucket
				for (int yHatMask = (1 << bucketSize) - 2; yHatMask >= 0; --yHatMask) {
					data_[bi + yHatMask] = Op::zero();
					// for each direct subideal
					for (int v = 0; v < bucketSize; ++v) {
						// variable mask
						int vMask = (1 << v);
						if (!(yHatMask & vMask)) {
							int tHatMask = yHatMask | vMask;
							int lastV = bucketSize - 1;
							int psMask = (yHatMask & ~(1 << lastV)) |
									((yHatMask >> (lastV - v)) & vMask);
							// variable index
							int i = b * family_.maxBucketSize + v;
							Op::add(data_[bi + yHatMask],
									Op::prod(psmap(i, psMask), data_[bi + tHatMask]));
						}
					}
				}
			}
		}
	};
};


std::ostream& operator<<(std::ostream& os, const BucketOrderFamily::Instance& instance) {
	instance.writeTo(os);
	return os;
}

std::istream& operator>>(std::istream& is, BucketOrderFamily::Instance& instance) {
	instance.readFrom(is);
	return is;
}


#endif


