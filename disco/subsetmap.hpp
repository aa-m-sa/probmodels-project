/*
 *  BEANDisco: SubsetMap class
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

#include "common.hpp"
#include "sortedarraysubset.hpp"

#ifndef SUBSETMAP_HPP
#define SUBSETMAP_HPP

unsigned long numKSubsets(unsigned long n, unsigned long k) {
	unsigned long sum = 1;
	unsigned long perLevel = 1;
	for (unsigned long i = 1; i <= k; ++i) {
		perLevel *= (n-i+1);
		perLevel /= i;
		sum += perLevel;
	}
	return sum;
}


template <typename T, typename I = size_t>
class FullSubsetMap {
private:
	ConstSortedArraySubset groundSet_;
	std::vector<int> elementIndices_;

	T* data_;

	I getSetMask(const SortedArraySubset& subset) {
		I setMask = 0;
		for (int i = 0; i < subset.size(); ++i)
			 setMask |= I(1) << elementIndices_[subset[i]];
		return setMask;
	}

		
	FullSubsetMap(const FullSubsetMap&) = delete; // disable copying
	FullSubsetMap& operator=(const FullSubsetMap&) = delete; // disable copying
public:
	FullSubsetMap(const ConstSortedArraySubset& groundSet) :
		groundSet_(groundSet),
		elementIndices_(groundSet.getLargest() + 1, -1)
	{
		for (int i = 0; i < groundSet_.size(); ++i)
			elementIndices_[groundSet_[i]] = i;
		data_ = new T[I(1) << groundSet_.size()];
	}
	
	~FullSubsetMap() {
		delete[] data_;
	}
	
	class SubsetIndex {
	private:
		I setMask_;
		SubsetIndex(I setMask) : setMask_(setMask) {}
	public:
		friend class FullSubsetMap;
	};

	template <typename SubsetElements>
	void getElements(const SubsetIndex& index, SubsetElements& subset) const {
		backTrackElements(index.datum_, subset);
	}

	class Subset {
	private:
		const FullSubsetMap& map_;
		I setMask_;
	public:
		Subset(const FullSubsetMap& map) :
			map_(map),
			setMask_(0)
		{}

		Subset(const FullSubsetMap& map, const SubsetIndex& index) :
			map_(map),
			setMask_(index.setMask_)
		{}

		void insert(int elem) {
			I eMask = I(1) << map_.elementIndices_[elem];
			assert(!(setMask_ & eMask));
			setMask_ ^= eMask;
		}
		void remove(int elem) {
			I eMask = I(1) << map_.elementIndices_[elem];
			assert(setMask_ & eMask);
			setMask_ ^= eMask;
		}
		void insertLargest(int elem) {
			assert((setMask_ >> map_.elementIndices_[elem]) == 0);
			insert(elem);
		}
		void removeLargest(int elem) {
			assert((setMask_ >> map_.elementIndices_[elem]) == 1);
			remove(elem);
		}

		void insertSmallest(int elem) {
			assert((setMask_ & ((I(1) << elem) - I(1))) == 0);
			insert(elem);
		}
		void removeSmallest(int elem) {
			assert((setMask_ & ((I(1) << elem) - I(1))) == 0);
			remove(elem);
		}

		SubsetIndex getIndex() const {
			return SubsetIndex(setMask_);
		}

		template <typename SubsetElements>
		void getElements(SubsetElements& subset) const {
			int n = map_.groundSet_.size();
			for (int i = 0; i < n; ++i) {
				if ((I(1) << i) & setMask_)
					subset.insertLargest(map_.groundSet_[i]);
			}
		}

		SortedArraySubset getElements() const {
			int n = map_.groundSet_.size();
			SortedArraySubset elems(n);
			getElements(elems);
			return elems;
		}

		friend class FullSubsetMap;
	};

	void setAll(const T& value) {
		int n = groundSet_.size();
		for (I i = 0; i < I(1) << n; ++i)
			data_[i] = value;
	}

	T& operator[] (const Subset& subset) {
		return data_[subset.setMask_];
	}

	const T& operator[] (const Subset& subset) const {
		return data_[subset.setMask_];
	}

	T& operator[] (const SubsetIndex& subsetIndex) {
		return data_[subsetIndex.setMask_];
	}

	const T& operator[] (const SubsetIndex& subsetIndex) const {
		return data_[subsetIndex.setMask_];
	}

	T& operator[] (const SortedArraySubset& subset) {
		I setMask = getSetMask(subset);
		return data_[setMask];
	}

	const T& operator[] (const SortedArraySubset& subset) const {
		I setMask = getSetMask(subset);
		return data_[setMask];
	}

	T getEmpty() const {
		return data_[0];
	}
	
	T getFull() const {
		return data_[(I(1) << groundSet_.size()) - I(1)];
	}

	template <typename Subset>
	struct ForSubsets {
	private:
		FullSubsetMap& map_;
		ConstSortedArraySubset set_;
		Subset subset_;
	
	public:
		ForSubsets(FullSubsetMap& map, const ConstSortedArraySubset& set, const Subset& subset) :
			map_(map),
			set_(set),
			subset_(subset)
		{
		}
		
		using Value = std::pair<const Subset&, T&>;

		struct Iterator {
		private:
			ForSubsets range_;
			int round_;
			I setMask_;
			Subset subset_;
		public:
			Iterator(const ForSubsets& range, int round, I setMask, Subset subset) :
				range_(range), round_(round), setMask_(setMask), subset_(subset)
			{}
			
			Iterator& operator++() {
				int i = 0;
				while (i < range_.set_.size()) {
					int elem = range_.set_[i];
					I eMask = I(1) << range_.map_.elementIndices_[elem];
					setMask_ ^= eMask;
					if (setMask_ & eMask) {
						subset_.insertSmallest(elem);
						return *this;
					}
					subset_.removeSmallest(elem);
					++i;
				}
				++round_;
				return *this;
			}
			
			bool operator!=(const Iterator& other) {
				return round_ != other.round_ || setMask_ != other.setMask_;
			}
			
			Value operator*() {
				return Value(subset_, range_.map_.data_[setMask_]);
			}
		};
		
		Iterator begin() {
			return Iterator(*this, 0, I(0), subset_);
		}
		
		Iterator end() {
			return Iterator(*this, 1, I(0), subset_);
		}
		
	};
	
	template <typename Subset>
	ForSubsets<Subset> forSubsetsOf(const ConstSortedArraySubset& set, const Subset& subsetInitializer) {
		return ForSubsets<Subset>(*this, set, subsetInitializer);
	}

	template <typename Subset>
	ForSubsets<Subset> forSubsets(const Subset& subsetInitializer) {
		return ForSubsets<Subset>(*this, groundSet_, subsetInitializer);
	}

	template <typename Op>
	void zetaTransform() {
		int n = groundSet_.size();
		for (int j = 0; j < n; ++j) {
			I eMask = I(1) << j;
			for (I setMask = 0; setMask < (I(1) << n); ++setMask) {
				if (!(eMask & setMask))
					Op::add(data_[setMask | eMask], data_[setMask]);
			}
		}
	}

	template <typename Op>
	void upZetaTransform() {
		int n = groundSet_.size();
		for (int j = 0; j < n; ++j) {
			I eMask = I(1) << j;
			for (I setMask = 0; setMask < (I(1) << n); ++setMask) {
				if (eMask & setMask)
					Op::add(data_[setMask ^ eMask], data_[setMask]);
			}
		}
	}

		
	template <typename Op>
	void forwardSum(const std::vector<FullSubsetMap<T>*>& predSetMaps) {
		data_[0] = Op::one();
		int n = groundSet_.size();
		for (I setMask = 1; setMask < (I(1) << n); ++setMask) {
			data_[setMask] = Op::zero();
			//I predSetMask = (setMask >> 1);
			for (int i = 0; i < n; ++i) {
				I eMask = (I(1) << i);
				if (setMask & eMask) {
					I subSetMask = setMask ^ eMask;
					I predSetMask = (subSetMask & (eMask - 1)) | ((subSetMask >> 1) & ~(eMask - 1));
					// variable index
					Op::add(data_[setMask],
							Op::prod(predSetMaps[i]->data_[predSetMask], data_[subSetMask]));
				}
				//predSetMask = (predSetMask & ~eMask) | (setMask & eMask);
			}
		}
	}


	template <typename Op>
	void backwardSum(const std::vector<FullSubsetMap<T>*>& predSetMaps) {
		int n = groundSet_.size();
		I fullSetMask = (I(1) << n) - I(1);
		data_[fullSetMask] = Op::one();
		for (I setMask = fullSetMask - 1; setMask != I(0) - I(1); --setMask) {
			data_[setMask] = Op::zero();
			for (int i = 0; i < n; ++i) {
				I eMask = (I(1) << i);
				if (!(setMask & eMask)) {
					I superSetMask = setMask ^ eMask;
					I predSetMask = (setMask & (eMask - 1)) | ((setMask >> 1) & ~(eMask - 1));
					Op::add(data_[setMask],
							Op::prod(predSetMaps[i]->data_[predSetMask], data_[superSetMask]));
				}
			}
		}
	}
};

/*template <typename T, typename I>
std::ostream& operator<<(std::ostream& os, const typename FullSubsetMap<T, I>::Subset& ss) {
	os << ss.getElements;
	return os;
}*/








template <typename T>
class SubsetMap {
private:
	ConstSortedArraySubset groundSet_;
	std::vector<int> elementIndices_;
	const int maxSubsetSize_;

	struct Datum {
		Datum* parent;
		Datum* child0;
		T value;
	};
	Datum* root_;

	void buildRecursive(Datum* x, int depth, int nextIndex, Datum*& freeDatums) {
		//for (int i = 0; i < setSize; ++i)
		//	printf("%d", inSet[i]);
		//printf(" => %d \n", x);

		if (depth == 0) {
			x->child0 = nullptr;
			return;
		}
		else {
			x->child0 = freeDatums - nextIndex;
			freeDatums += groundSet_.size() - nextIndex;
			for (int i = nextIndex; i < groundSet_.size(); ++i) {
				Datum* child = x->child0 + i;
				child->parent = x;
				buildRecursive(child, depth - 1, i + 1, freeDatums);
			}
		}

	}

	Datum* getDatum(const SortedSubsetRange& subset) const {
		assert(subset.size() <= maxSubsetSize_);
		Datum* x = root_;
		for (int i = 0; i < subset.size(); ++i)
			x = x->child0 + elementIndices_[subset[i]];
		return x;
	}

	template <typename SubsetElements>
	void backTrackElementsRecursive(Datum* x, SubsetElements& subset) const {
		if (x != root_) {
			Datum* parent = x->parent;
			backTrackElementsRecursive(parent, subset);
			subset.insertLargest(groundSet_[x - parent->child0]);
		}
	}

	template <typename SubsetElements>
	void backTrackElements(Datum* x, SubsetElements& subset) const {
		subset.clear();
		backTrackElementsRecursive(x, subset);
		/*int size = 0;
		while (x != root_) {
			x = x->parent;
			++size;
		}
		subset.resize(size);
		int i = 0;
		while (x != root_) {
			Datum* parent = x->parent;
			subset[i] = groundSet_[x - parent->child0];
			x = parent;
			++i;
		}*/
	}
		
	SubsetMap(const SubsetMap&) = delete; // disable copying
	SubsetMap& operator=(const SubsetMap&) = delete; // disable copying
public:
	SubsetMap(const ConstSortedArraySubset& groundSet, int maxSubsetSize) :
		groundSet_(groundSet),
		//elementIndices_(max(groundSet), -1),
		elementIndices_(groundSet.getLargest() + 1, -1),
		maxSubsetSize_(maxSubsetSize)
	{
		//std::sort(groundSet_.begin(), groundSet_.end());
		for (int i = 0; i < groundSet_.size(); ++i)
			elementIndices_[groundSet_[i]] = i;
		unsigned long nSubsets = numKSubsets(groundSet_.size(), maxSubsetSize_);
		Datum* freeDatums = new Datum[nSubsets];
		root_ = freeDatums;
		root_->parent = nullptr;
		freeDatums += 1;
		buildRecursive(root_, maxSubsetSize_, 0, freeDatums);
		
	}
	
	~SubsetMap() {
		delete[] root_;
	}
	
	class SubsetIndex {
	private:
		Datum* datum_;
		SubsetIndex(Datum* datum) : datum_(datum) {}
	public:
		bool operator==(const SubsetIndex& other) const {
			return datum_ == other.datum_;
		}
		friend class SubsetMap;
	};

	template <typename SubsetElements>
	void getElements(const SubsetIndex& index, SubsetElements& subset) const {
		backTrackElements(index.datum_, subset);
	}

	SubsetIndex getIndex(const SortedSubsetRange& subset) const {
		return SubsetIndex(getDatum(subset));
	}

	SubsetIndex getEmptyIndex() const {
		return SubsetIndex(root_);
	}

	unsigned long numSubsets() const {
		return numKSubsets(groundSet_.size(), maxSubsetSize_);
	}

	class Subset {
	private:
		const SubsetMap& map_;
		Datum* datum_;
	public:
		Subset(const SubsetMap& map) :
			map_(map),
			datum_(map.root_)
		{}

		Subset(const SubsetMap& map, const SubsetIndex& index) :
			map_(map),
			datum_(index.datum_)
		{}

		void insertLargest(int elem) {
			datum_ = datum_->child0 + map_.elementIndices_[elem];
		}
		void removeLargest(int elem) {
			assert(datum_->parent != nullptr);
			assert(getLargest() == elem);
			datum_ = datum_->parent;
		}
		int getLargest() const {
			assert(datum_->parent != nullptr);
			return map_.groundSet_[datum_ - datum_->parent->child0];
		}

		SubsetIndex getIndex() const {
			return SubsetIndex(datum_);
		}

		template <typename SubsetElements>
		void getElements(SubsetElements& subset) const {
			map_.backTrackElements(datum_, subset);
		}

		friend class SubsetMap;
	};

	T& operator[] (const Subset& subset) {
		return subset->datum_->value;
	}

	const T& operator[] (const Subset& subset) const {
		return subset->datum_->value;
	}

	T& operator[] (const SubsetIndex& subsetIndex) {
		return subsetIndex.datum_->value;
	}

	const T& operator[] (const SubsetIndex& subsetIndex) const {
		return subsetIndex.datum_->value;
	}

	T& operator[] (const SortedSubsetRange& subset) {
		Datum* x = getDatum(subset);
		return x->value;
	}

	const T& operator[] (const SortedSubsetRange& subset) const {
		Datum* x = getDatum(subset);
		return x->value;
	}

	template <typename Subset>
	struct ForSubsets {
	private:
		SubsetMap& map_;
		ConstSortedArraySubset set_;
		Subset subset_;
	
	public:
		ForSubsets(SubsetMap& map, const ConstSortedArraySubset& set, const Subset& subset) :
			map_(map),
			set_(set),
			subset_(subset)
		{
		}
		
		using Value = std::pair<const Subset&, T&>;

		struct Iterator {
		private:
			ForSubsets range_;
			int round_;
			SortedArraySubset ssi_;
			Datum* datum_;
			Subset subset_;
		public:
			Iterator(const ForSubsets& range, int round, const SortedArraySubset& ssi,
					Datum* datum, Subset subset) :
				range_(range), round_(round), ssi_(ssi), datum_(datum), subset_(subset)
			{
				assert(datum_ != nullptr);
			}
			
			Iterator& operator++() {
				int next = ssi_.empty() ? 0 : ssi_.getLargest() + 1;
				if (next < range_.set_.size()) {
					if (ssi_.size() < range_.map_.maxSubsetSize_) {
						ssi_.insertLargest(next);
						datum_ = datum_->child0 + range_.map_.elementIndices_[range_.set_[next]];
						subset_.insertLargest(range_.set_[next]);
						return *this;
					}
				}
				else {
					if (!ssi_.empty()) {
						int prev = ssi_.getLargest();
						ssi_.removeLargest(prev);
						datum_ = datum_->parent;
						assert(datum_ != nullptr);
						subset_.removeLargest(range_.set_[prev]);
					}
				}
				if (ssi_.empty()) {
					assert(datum_ == range_.map_.root_);
					++round_;
					return *this;
				}
				int prev = ssi_.getLargest();
				ssi_.removeLargest(prev);
				datum_ = datum_->parent;
				assert(datum_ != nullptr);
				subset_.removeLargest(range_.set_[prev]);
				next = prev + 1;
				ssi_.insertLargest(next);
				datum_ = datum_->child0 + range_.map_.elementIndices_[range_.set_[next]];
				subset_.insertLargest(range_.set_[next]);
				return *this;
			}
			
			bool operator!=(const Iterator& other) {
				return round_ != other.round_ || ssi_ != other.ssi_;
			}
			
			Value operator*() {
				return Value(subset_, datum_->value);
			}
		};
		
		Iterator begin() {
			return Iterator(*this, 0, SortedArraySubset(set_.size() + 1), map_.root_, subset_);
		}
		
		Iterator end() {
			return Iterator(*this, 1, SortedArraySubset(set_.size() + 1), map_.root_, subset_);
		}
		
	};
	
	template <typename Subset>
	ForSubsets<Subset> forSubsetsOf(const ConstSortedArraySubset& set, const Subset& subsetInitializer) {
		return ForSubsets<Subset>(*this, set, subsetInitializer);
	}

	template <typename Subset>
	ForSubsets<Subset> forSubsets(const Subset& subsetInitializer) {
		return ForSubsets<Subset>(*this, groundSet_, subsetInitializer);
	}

	template <typename Subset>
	struct ForSubsetsConst {
	private:
		const SubsetMap& map_;
		ConstSortedArraySubset set_;
		Subset subset_;
	
	public:
		ForSubsetsConst(const SubsetMap& map, const ConstSortedArraySubset& set, const Subset& subset) :
			map_(map),
			set_(set),
			subset_(subset)
		{
		}
		
		using Value = std::pair<const Subset&, const T&>;

		struct Iterator {
		private:
			ForSubsetsConst range_;
			int round_;
			SortedArraySubset ssi_;
			Datum* datum_;
			Subset subset_;
		public:
			Iterator(const ForSubsetsConst& range, int round, const SortedArraySubset& ssi,
					Datum* datum, Subset subset) :
				range_(range), round_(round), ssi_(ssi), datum_(datum), subset_(subset)
			{
				assert(datum_ != nullptr);
			}
			
			Iterator& operator++() {
				int next = ssi_.empty() ? 0 : ssi_.getLargest() + 1;
				if (next < range_.set_.size()) {
					if (ssi_.size() < range_.map_.maxSubsetSize_) {
						ssi_.insertLargest(next);
						datum_ = datum_->child0 + range_.map_.elementIndices_[range_.set_[next]];
						subset_.insertLargest(range_.set_[next]);
						return *this;
					}
				}
				else {
					if (!ssi_.empty()) {
						int prev = ssi_.getLargest();
						ssi_.removeLargest(prev);
						datum_ = datum_->parent;
						assert(datum_ != nullptr);
						subset_.removeLargest(range_.set_[prev]);
					}
				}
				if (ssi_.empty()) {
					assert(datum_ == range_.map_.root_);
					++round_;
					return *this;
				}
				int prev = ssi_.getLargest();
				ssi_.removeLargest(prev);
				datum_ = datum_->parent;
				assert(datum_ != nullptr);
				subset_.removeLargest(range_.set_[prev]);
				next = prev + 1;
				ssi_.insertLargest(next);
				datum_ = datum_->child0 + range_.map_.elementIndices_[range_.set_[next]];
				subset_.insertLargest(range_.set_[next]);
				return *this;
			}
			
			bool operator!=(const Iterator& other) {
				return round_ != other.round_ || ssi_ != other.ssi_;
			}
			
			Value operator*() {
				return Value(subset_, datum_->value);
			}
		};
		
		Iterator begin() {
			return Iterator(*this, 0, SortedArraySubset(set_.size() + 1), map_.root_, subset_);
		}
		
		Iterator end() {
			SortedArraySubset ssi(set_.size() + 1);
			ssi.insertLargest(0);
			return Iterator(*this, 1, SortedArraySubset(set_.size() + 1), map_.root_, subset_);
		}
		
	};
	
	template <typename Subset>
	ForSubsetsConst<Subset> forSubsetsOf(const ConstSortedArraySubset& set, const Subset& subsetInitializer) const {
		return ForSubsetsConst<Subset>(*this, set, subsetInitializer);
	}

	template <typename Subset>
	ForSubsetsConst<Subset> forSubsets(const Subset& subsetInitializer) const {
		return ForSubsetsConst<Subset>(*this, groundSet_, subsetInitializer);
	}
};


#endif
