/*
 *  BEANDisco: SortedArraySubset class
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


#include <cassert>
#include <iostream>


#ifndef SORTEDARRAYSUBSET_HPP
#define SORTEDARRAYSUBSET_HPP

struct SortedSubsetRange {
protected:
	size_t size_;
	int* elements_;
	
	SortedSubsetRange& operator=(const SortedSubsetRange&) = delete; // disable assignment

	SortedSubsetRange() {
	}
public:
	SortedSubsetRange(size_t size, int* elements) {
		size_ = size;
		elements_ = elements;
	}
	
	/*SortedSubsetRange(const SortedSubsetRange& other) {
		size_ = other.size_;
		elements_ = other.elements_;
	}*/
	
	bool operator==(const SortedSubsetRange& other) const {
		if (size_ != other.size_)
			return false;
		for (int i = 0; i < size_; ++i)
			if (elements_[i] != other.elements_[i])
				return false;
		return true;
	}
	
	bool operator!=(const SortedSubsetRange& other) const {
		return !operator==(other);
	}

	size_t size() const {
		return size_;
	}
	
	bool empty() const {
		return size_ == 0;
	}
	
	int operator[] (int i) const {
		assert(0 <= i && i < size_);
		return elements_[i];
	}
	
	int getLargest() const {
		assert(size_ > 0);
		return elements_[size_ - 1];
	}

	bool contains(int x) const {
		for (int i = 0; i < size_; ++i)
			if (x == elements_[i])
				return true;
		return false;
	}

	//typedef int* Iterator;
	struct ConstIterator {
	private:
		int* elem_;
	public:
		ConstIterator(int* elem) {
			elem_ = elem;
		}
		ConstIterator& operator++() {
			++elem_;
			return *this;
		}
		
		bool operator!=(const ConstIterator& other) {
			return elem_ != other.elem_;
		}
		
		int operator*() {
			return *elem_;
		}
	};

	ConstIterator begin() const {
		return ConstIterator(elements_);
	}
		
	ConstIterator end() const {
		return ConstIterator(elements_ + size_);
	}

	friend std::ostream& operator<<(std::ostream& os, const SortedSubsetRange& ss);
};



struct SortedArraySubset : public SortedSubsetRange {
private:
	size_t maxSize_;
	
	SortedArraySubset& operator=(const SortedArraySubset&) = delete; // disable assignment
public:
	SortedArraySubset() {
		maxSize_ = 0;
		elements_ = NULL;
		size_ = 0;
	}
	
	SortedArraySubset(int maxSize) {
		maxSize_ = maxSize;
		elements_ = new int[maxSize_];
		size_ = 0;
	}
	
	SortedArraySubset(const SortedArraySubset& other) :
		SortedSubsetRange(other)
	{
		maxSize_ = other.maxSize_;
		elements_ = new int[maxSize_];
		size_ = 0;
		for (; size_ < other.size_; ++size_)
			elements_[size_] = other.elements_[size_];
	}
	
	SortedArraySubset(SortedArraySubset&& other) {
		maxSize_ = other.maxSize_;
		elements_ = other.elements_;
		size_ = other.size_;
		other.maxSize_ = 0;
		other.elements_ = nullptr;
		other.size_ = 0;
	}

	SortedArraySubset(int maxSize, const SortedSubsetRange& other) :
		SortedSubsetRange(other)
	{
		maxSize_ = maxSize;
		elements_ = new int[maxSize_];
		size_ = 0;
		for (; size_ < other.size(); ++size_)
			elements_[size_] = other[size_];
	}
	
	SortedArraySubset(int maxSize, std::vector<int> elems) {
		maxSize_ = maxSize;
		assert(elems.size() <= maxSize);
		std::sort(elems.begin(), elems.end());
		elements_ = new int[maxSize_];
		size_ = 0;
		for (int i = 0; i < elems.size(); ++i)
			insertLargest(elems[i]);
	}
	
	~SortedArraySubset() {
		delete[] elements_;
	}
	
	void clear() {
		size_ = 0;
	}
	
	// TODO: ei pitäisi olla sallittu
	int& operator[] (int i) {
		assert(0 <= i && i < size_);
		return elements_[i];
	}

	int operator[] (int i) const {
		assert(0 <= i && i < size_);
		return elements_[i];
	}
	
	void insertLargest(int elem) {
		assert(size_ < maxSize_);
		assert(size_ == 0 || elements_[size_ - 1] < elem);
		elements_[size_++] = elem;
	}

	void removeLargest(int elem) {
		assert(size_ > 0);
		--size_;
		assert(elements_[size_] == elem);
	}
	
	void insert(int elem) {
		assert(size_ < maxSize_);
		int i = size_;
		++size_;
		while(i > 0 && elements_[i - 1] > elem) {
			elements_[i] = elements_[i - 1];
			--i;
		}
		assert(i == 0 || elements_[i - 1] < elem);
		elements_[i] = elem;
	}
	
	void remove(int elem) {
		assert(size_ > 0);
		--size_;
		int i = size_;
		int curr = elements_[i];
		--i;
		while(i >= 0 && curr > elem) {
			int tmp = elements_[i];
			elements_[i] = curr;
			curr = tmp;
			--i;
		}
		assert(curr == elem);
	}

	static SortedArraySubset fullSet(int n) {
		SortedArraySubset s(n);
		for (int i = 0; i < n; ++i)
			s.insertLargest(i);
		return s;
	}
	
	//friend std::ostream& operator<<(std::ostream& os, const SortedArraySubset& ss);
};/**/



/*struct SortedArraySubset {
private:
	size_t maxSize_;
	size_t size_;
	int* elements_;
	
	SortedArraySubset& operator=(const SortedArraySubset&) = delete; // disable assignment
public:
	SortedArraySubset() {
		maxSize_ = 0;
		elements_ = NULL;
		size_ = 0;
	}
	
	SortedArraySubset(int maxSize) {
		maxSize_ = maxSize;
		elements_ = new int[maxSize_];
		size_ = 0;
	}
	
	SortedArraySubset(const SortedArraySubset& other) {
		maxSize_ = other.maxSize_;
		elements_ = new int[maxSize_];
		size_ = 0;
		for (; size_ < other.size_; ++size_)
			elements_[size_] = other.elements_[size_];
	}
	
	SortedArraySubset(SortedArraySubset&& other) {
		maxSize_ = other.maxSize_;
		elements_ = other.elements_;
		size_ = other.size_;
		other.maxSize_ = 0;
		other.elements_ = nullptr;
		other.size_ = 0;
	}

	SortedArraySubset(int maxSize, std::vector<int> elems) {
		maxSize_ = maxSize;
		assert(elems.size() <= maxSize);
		std::sort(elems.begin(), elems.end());
		elements_ = new int[maxSize_];
		size_ = 0;
		for (int i = 0; i < elems.size(); ++i)
			insertLargest(elems[i]);
	}
	
	~SortedArraySubset() {
		delete[] elements_;
	}
	
	bool operator==(const SortedArraySubset& other) const {
		if (size_ != other.size_)
			return false;
		for (int i = 0; i < size_; ++i)
			if (elements_[i] != other.elements_[i])
				return false;
		return true;
	}
	
	bool operator!=(const SortedArraySubset& other) const {
		return !operator==(other);
	}

	size_t size() const {
		return size_;
	}
	
	bool empty() const {
		return size_ == 0;
	}
	
	void clear() {
		size_ = 0;
	}
	
	// TODO: ei pitäisi olla sallittu
	int& operator[] (int i) {
		assert(0 <= i && i < size_);
		return elements_[i];
	}

	int operator[] (int i) const {
		assert(0 <= i && i < size_);
		return elements_[i];
	}
	
	void insertLargest(int elem) {
		assert(size_ < maxSize_);
		assert(size_ == 0 || elements_[size_ - 1] < elem);
		elements_[size_++] = elem;
	}

	void removeLargest(int elem) {
		assert(size_ > 0);
		--size_;
		assert(elements_[size_] == elem);
	}

	int getLargest() const {
		assert(size_ > 0);
		return elements_[size_ - 1];
	}
	
	void insert(int elem) {
		assert(size_ < maxSize_);
		int i = size_;
		++size_;
		while(i > 0 && elements_[i - 1] > elem) {
			elements_[i] = elements_[i - 1];
			--i;
		}
		assert(i == 0 || elements_[i - 1] < elem);
		elements_[i] = elem;
	}
	
	void remove(int elem) {
		assert(size_ > 0);
		--size_;
		int i = size_;
		int curr = elements_[i];
		--i;
		while(i >= 0 && curr > elem) {
			int tmp = elements_[i];
			elements_[i] = curr;
			curr = tmp;
			--i;
		}
		assert(curr == elem);
	}
	
	bool contains(int x) const {
		for (int i = 0; i < size_; ++i)
			if (x == elements_[i])
				return true;
		return false;
	}

	static SortedArraySubset fullSet(int n) {
		SortedArraySubset s(n);
		for (int i = 0; i < n; ++i)
			s.insertLargest(i);
		return s;
	}
	
	friend std::ostream& operator<<(std::ostream& os, const SortedArraySubset& ss);
};/**/



struct ConstSortedArraySubset : public SortedSubsetRange {
private:
	ConstSortedArraySubset& operator=(const ConstSortedArraySubset&) = delete; // disable assignment

public:
	ConstSortedArraySubset() {
		size_ = 0;
		elements_ = nullptr;
	}
	
	ConstSortedArraySubset(const SortedSubsetRange& other) {
		size_ = other.size();
		elements_ = new int[size_];
		for (int i = 0; i < size_; ++i)
			elements_[i] = other[i];
	}

	// should be unnecessary! (the above should suffice)
	ConstSortedArraySubset(const ConstSortedArraySubset& other) :
		SortedSubsetRange(other.size(), new int[other.size()])
	{
		for (int i = 0; i < size_; ++i)
			elements_[i] = other[i];
	}

	// should be unnecessary! (the above should suffice)
	ConstSortedArraySubset(const SortedArraySubset& other) {
		size_ = other.size();
		elements_ = new int[size_];
		for (int i = 0; i < size_; ++i)
			elements_[i] = other[i];
	}
	
	ConstSortedArraySubset(ConstSortedArraySubset&& other) {
		size_ = other.size_;
		elements_ = other.elements_;
		other.size_ = 0;
		other.elements_ = nullptr;
	}

	ConstSortedArraySubset(std::vector<int> elems) {
		std::sort(elems.begin(), elems.end());
		size_ = elems.size();
		elements_ = new int[size_];
		for (int i = 0; i < size_; ++i)
			elements_[i] = elems[i];
	}
	
	~ConstSortedArraySubset() {
		delete[] elements_;
	}
	
	friend std::ostream& operator<<(std::ostream& os, const ConstSortedArraySubset& ss);
};/**/




/*struct ConstSortedArraySubset {
private:
	size_t size_;
	int* elements_;
	
	ConstSortedArraySubset& operator=(const ConstSortedArraySubset&) = delete; // disable assignment
public:
	ConstSortedArraySubset() {
		size_ = 0;
		elements_ = NULL;
	}
	
	ConstSortedArraySubset(const ConstSortedArraySubset& other) {
		size_ = other.size_;
		elements_ = new int[size_];
		for (int i = 0; i < size_; ++i)
			elements_[i] = other.elements_[i];
	}

	ConstSortedArraySubset(const SortedArraySubset& other) {
		size_ = other.size();
		elements_ = new int[size_];
		for (int i = 0; i < size_; ++i)
			elements_[i] = other[i];
	}
	
	ConstSortedArraySubset(ConstSortedArraySubset&& other) {
		size_ = other.size_;
		elements_ = other.elements_;
		other.size_ = 0;
		other.elements_ = nullptr;
	}

	ConstSortedArraySubset(std::vector<int> elems) {
		std::sort(elems.begin(), elems.end());
		size_ = elems.size();
		elements_ = new int[size_];
		for (int i = 0; i < size_; ++i)
			elements_[i] = elems[i];
	}
	
	~ConstSortedArraySubset() {
		delete[] elements_;
	}
	
	bool operator==(const ConstSortedArraySubset& other) const {
		if (size_ != other.size_)
			return false;
		for (int i = 0; i < size_; ++i)
			if (elements_[i] != other.elements_[i])
				return false;
		return true;
	}
	
	bool operator!=(const ConstSortedArraySubset& other) const {
		return !operator==(other);
	}

	size_t size() const {
		return size_;
	}
	
	bool empty() const {
		return size_ == 0;
	}
	
	int operator[] (int i) const {
		assert(0 <= i && i < size_);
		return elements_[i];
	}
	
	int getLargest() const {
		assert(size_ > 0);
		return elements_[size_ - 1];
	}

	bool contains(int x) const {
		for (int i = 0; i < size_; ++i)
			if (x == elements_[i])
				return true;
		return false;
	}

	friend std::ostream& operator<<(std::ostream& os, const ConstSortedArraySubset& ss);
};/**/



bool isSubsetOf(const SortedArraySubset& a, const SortedArraySubset& b) {
	int i = 0;
	int j = 0;
	while (i < a.size()) {
		if (j >= b.size()) 
			return false;
		while (b[j] < a[i]) {
			++j;
			if (j >= b.size()) 
				return false;
		}
		if (a[i] != b[j])
			return false;
		++i;
		++j;
	}
	return true;
}


std::ostream& operator<<(std::ostream& os, const SortedSubsetRange& ss) {
	if (ss.size() == 0)
		os << "∅";
	else {
		os << "{" << ss.elements_[0];
		for (int i = 1; i < ss.size(); ++i)
			os << "," << ss.elements_[i];
		os << "}";
	}
	//for (int i = s.length; i < width; ++i)
	//	os << "  ";
	return os;
}


#endif

