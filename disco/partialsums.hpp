/*
 *  BEANDisco: algorithms for computing sums of subsequences
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
#include <vector>

#ifndef PARTIALSUMS_HPP
#define PARTIALSUMS_HPP


template <typename T>
class PartialSumsNaive {
private:
	T zero_;
	std::vector<T> values_;
public:
	PartialSumsNaive(const T& zero) :
		zero_(zero)
	{
		reset();
	}
	void reset() {
		values_.clear();
	}
	void add(const T& x) {
		values_.push_back(x);
	}

	size_t size() const {
		return values_.size();
	}

	T get(size_t begin, size_t end) const {
		assert(begin <= end && end <= size());
		T sum = 0;
		//printf("naive: ");
		for (size_t i = begin; i < end; ++i) {
			sum += values_[i];
			//printf("%d ", values_[i]);
		}
		//printf("\n");
		return sum;
	}

};/**/

template <typename T>
class PartialSumsSubtract {
private:
	T zero_;
	std::vector<T> cumSum_;
public:
	PartialSumsSubtract(const T& zero) :
		zero_(zero)
	{
		reset();
	}
	void reset() {
		cumSum_.clear();
		cumSum_.push_back(zero_);
	}
	void add(const T& x) {
		cumSum_.push_back(cumSum_.back() + x);
	}

	size_t size() const {
		return cumSum_.size() - 1;
	}

	T get(size_t begin, size_t end) const {
		assert(begin <= end && end <= size());
		return cumSum_[end] - cumSum_[begin];
	}

};/**/


template <typename T>
class PartialSumsTree {
private:
	T zero_;
	std::vector<std::vector<T>> sums_;
	size_t n;
public:
	PartialSumsTree(const T& zero) :
		zero_(zero)
	{
		reset();
	}
	void reset() {
		n = 0;
		sums_.clear();
	}
	void add(const T& x) {
		int i = 0;
		size_t m = n;
		T s = x;
		while (true) {
			if (i == sums_.size())
				sums_.push_back(std::vector<T>());
			sums_[i].push_back(s);
			if ((m & 1) == 0)
				break;
			s += sums_[i][m - 1];
			m /= 2;
			++i;
		}
		++n;
	}

	size_t size() const {
		return n;
	}

	void disp() {
		int i = 0;
		for (auto sums : sums_) {
			printf("%d: ", i);
			for (auto sum : sums) {
				for (int j = 0; j < (1 << i) - 1; ++j)
					printf(" ");
				printf("%d ", sum);
				for (int j = 0; j < (1 << i) - 1; ++j)
					printf(" ");
			}
			printf("\n");
			++i;
		}
	}

	T get(size_t begin, size_t end) const {
		assert(begin <= end && end <= n);
		T sum = zero_;
		//printf("other: ");
		int i = 0;
		size_t m = 1;
		while (begin + m < end) {
			if (begin & m) {
				sum += sums_[i][begin >> i];
				begin += m;
				//printf("%d ", sums_[i][begin >> i]);
			}
			i += 1;
			m *= 2;
		}
		do {
			if (begin + m <= end) {
				sum += sums_[i][begin >> i];
				//printf("%d ", sums_[i][begin >> i]);
				begin += m;
			}
			i -= 1;
			m /= 2;
		} while (i >= 0);
		//printf("\n");
		return sum;
	}

};/**/



/*template <typename T>
class PartialSums {
private:
	T zero_;
	std::vector<std::vector<T>> sums_;
	size_t n;
public:
	PartialSums(const T& zero) :
		zero_(zero)
	{
		reset();
	}
	void reset() {
		n = 0;
		sums_.clear();
	}
	void add(const T& x) {
		int i = 0;
		size_t m = n;
		T s = x;
		while (true) {
			if (i == sums_.size())
				sums_.push_back(std::vector<T>());
			sums_[i].push_back(s);
			if ((m & 1) == 0)
				break;
			s += sums_[i][m];
			m /= 2;
			++i;
		}
		++n;
	}
	void disp() {
		int i = 0;
		for (auto sums : sums_) {
			printf("%d: ", i);
			for (auto sum : sums) {
				for (int j = 0; j < (1 << i) - 1; ++j)
					printf(" ");
				printf("%d ", sum);
				for (int j = 0; j < (1 << i) - 1; ++j)
					printf(" ");
			}
			printf("\n");
			++i;
		}
	}

	class SlidingSum {
	private:
		const PartialSums& partialSums_;
		size_t begin_;
		size_t end_;
		std::vector<T> leftSums_;
		std::vector<T> rightSums_;
	public:
		SlidingSum(const PartialSums& partialSums) :
			partialSums_(partialSums)
		{
			left_ = 0;
			right_ = 0;
			//leftSums_.clear();
			//rightSums_.clear();
		}

		void moveBegin() {

		}

		void moveEnd() {
			int i = 0;
			size_t m = end_;
			while (m & 1) {
				if (i == rightSums_.size()) {
					leftSums_.push_back(std::vector<T>());
					rightSums_.push_back(std::vector<T>());
				}
				m /= 2;
				++i;
			}
			++end_;
		}
			//size_t k = size_t(ceil(n * tailFraction_));
			//size_t s = n - k;
			//return (cumSum_.back() - cumSum_[cumSum_.size() - 1 - n]) / n;
		T get() const {
			
		}
	};

};/**/


#endif

