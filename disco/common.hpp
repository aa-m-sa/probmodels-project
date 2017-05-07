/*
 *  BEANDisco: general stuff
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

#include <memory>
#include <algorithm>
#include <random>
#include <boost/tokenizer.hpp>

#include "format.hpp"

#ifndef COMMON_HPP
#define COMMON_HPP

//#define NDEBUG
#include <cassert>

// make_unique definition
/*namespace std {
	template<typename T, typename ...Args>
	std::unique_ptr<T> make_unique(Args&& ...args) {
	    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
	}
}*/

// min and max utility function
template <typename T>
T min(T x, T y) {
	return x < y ? x : y;
}

template <typename T>
T max(T x, T y) {
	return x > y ? x : y;
}

template <typename R>
const typename R::value_type& max(const R& range) {
	return *(std::max_element(range.begin(), range.end()));
}


// assign the logarithm of a variable
template <typename T, typename F>
void setLog(T& x, F l){
	return x = exp(l);
}

// binomial coefficient
template <typename T>
T binom(int n, int k) {
	T res;
	setLog(res, lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1));
	return res;
}



// create a random number generator
typedef std::mt19937 Rng;
Rng rng;

// random number generation
double randu() {
	std::uniform_real_distribution<double> dist(0, 1);
	return dist(rng);
}

int randuint(int ceiling) {
	std::uniform_int_distribution<int> dist(0, ceiling - 1);
	return dist(rng);
}

/*long randulong(long ceiling) {
	std::uniform_int_distribution<long> dist((long)0, ceiling - 1);
	return dist(rng);
}*/


void randomShuffle(std::vector<int>& perm) {
	int n = perm.size();
	for (int i = 0; i < n; ++i) {
		int j = i + randuint(n - i);
		std::swap(perm[i], perm[j]);
	}
}


/**
 * An implementation of the discrete distribution using the alias method: O(n) initialization time
 * and O(1) per-sample time.
 */
template <typename T>
class DiscreteDist {
private:
	T* prob_;
	size_t* alias_;
	std::uniform_real_distribution<T> uniDist_;
public:
	DiscreteDist(const std::vector<T>& probs)
		: uniDist_(0, probs.size())
	{
		size_t n = probs.size();
		prob_ = new T[n];
		alias_ = new size_t[n];
		
		size_t* list = new size_t[n];
		size_t s = 0;
		size_t l = n;
		for (size_t i = 0; i < n; ++i) {
			alias_[i] = i; // does not matter in theory
			prob_[i] = probs[i] * n;
			if (prob_[i] < 1)
				list[s++] = i;
			else
				list[--l] = i;
		}
		while (s > 0 && l < n) {
			size_t small = list[--s];
			size_t large = list[l];
			alias_[small] = large;
			prob_[large] = (prob_[large] + prob_[small]) - 1;
			if (prob_[large] < 1) {
				++l;
				list[s++] = large;
			}
		}
		while (l < n)
			prob_[list[l++]] = 1;
		while (s > 0)
			prob_[list[--s]] = 1;
		
		delete[] list;
	}
	
	~DiscreteDist() {
		delete[] prob_;
		delete[] alias_;
	}
	
	size_t rand() {
		T i, f;
		f = uniDist_(rng);

		assert(0 <= f && f < uniDist_.max());
		f = modf(f, &i);
		assert(0 <= i && i < uniDist_.max());
		size_t j = (size_t)i;
		if (f <= prob_[j])
			return j;
		else
			return alias_[j];
	}
};


template <typename T>
int randDiscrete(const std::vector<T>& probs) {
	double u = randu();
	int i = 0;
	while (i < probs.size() - 1) {
		u -= static_cast<double>(probs[i]);
		if (u <= 0)
			break;
		++i;
	}
	return i;
}


/*template <typename T>
void normalizeProbs(std::vector<T> probs) {
	T normConst = 0.0;
	for (int i = 0; i < probs.size(); ++i)
		normConst += probs[i];
	T invNormConst = 1.0 / normConst;
	for (int i = 0; i < probs.size(); ++i)
		probs[i] *= invNormConst;
}*/

template <typename U, typename T>
void getNormalizedProbs(const std::vector<U>& unProbs, std::vector<T>& probs) {
	assert(unProbs.size() == probs.size());
	U normConst = 0.0;
	for (int i = 0; i < unProbs.size(); ++i)
		normConst += unProbs[i];
	U invNormConst = U(1.0) / normConst;
	for (int i = 0; i < unProbs.size(); ++i)
		probs[i] = T(unProbs[i] * invNormConst);
}



template <typename SampleSpace>
class MCProposalDist {
public:
	//class Proposal { virtual void apply(SampleSpace::Instance& instance) const = 0; };
	//virtual Proposal randProposal(SampleSpace::Instance& instance) const = 0;
	virtual void randStep(typename SampleSpace::Instance& instance) const = 0;
};



template <class T>
class SumProdOp {
public:
	static T sum(const T& x, const T& y) {
		return x + y;
	}
	static void add(T& x, const T& y) {
		x += y;
	}
	static T prod(const T& x, const T& y) {
		return x * y;
	}
	
	static constexpr T zero() { return T(0); }
	static constexpr T one() { return T(1); }
};

template <class T>
class MaxProdOp {
public:
	static T sum(const T& x, const T& y) {
		return std::max(x, y);
	}
	static void add(T& x, const T& y) {
		if (y > x)
			return y;
		else
			return x;
	}
	static T prod(const T& x, const T& y) {
		return x * y;
	}

	static constexpr T zero() { 
		return std::numeric_limits<T>::is_signed && std::numeric_limits<T>::has_infinity ?
				-std::numeric_limits<T>::infinity() : std::numeric_limits<T>::lowest();
	}
	static constexpr T one() { return T(1); };
};

template <class T>
class RandSelector {
private:
	std::vector<int> items_;
	std::vector<double> probs_;
	T totalUnnormProb_;
public:
	RandSelector(int n) {
		items_.reserve(n);
		probs_.reserve(n);
	}
	void init(T sum) {
		items_.clear();
		probs_.clear();
		totalUnnormProb_ = sum;
	}
	void add(int item, T unnormProb) {
		items_.push_back(item);
		probs_.push_back(double(unnormProb / totalUnnormProb_));
	}
	int get() {
		return items_[randDiscrete(probs_)];
	}
};

template <class T>
class MaxSelector {
private:
	int currentItem_;
	T currentValue_;
public:
	MaxSelector(int n) {
	}
	void init(T max) {
		currentItem_ = -1;
		currentValue_ = std::numeric_limits<T>::is_signed && std::numeric_limits<T>::has_infinity ?
				-std::numeric_limits<T>::infinity() : std::numeric_limits<T>::lowest();
	}
	void add(int item, T value) {
		if (value > currentValue_) {
			currentItem_ = item;
			currentValue_ = value;
		}
	}
	int get() {
		return currentItem_;
	}
};



/**
 * Splits string by separator.
 */
using Tokens = boost::tokenizer<boost::char_separator<char>>;
Tokens split(const std::string& str, const std::string& sep,
		bool keepEmpty = false) {
	return Tokens(str, boost::char_separator<char>(sep.c_str(), "",
			keepEmpty ? boost::keep_empty_tokens : boost::drop_empty_tokens));
}

/**
 * Tests if a range contains a given value.
 */
template <typename Range, typename T>
bool contains(const Range& range, const T& value) {
//std::cout << "contains(" << range << ", " << value << ") -> "  <<
//(std::find(range.begin(), range.end(), value) != range.end()) << std::endl;
//std::cout << "length(range) = " << range.size() << std::endl;
	return std::find(range.begin(), range.end(), value) != range.end();
}



/**
 * General exception.
 */
class Exception : public std::exception {
protected:
	std::string msg_;
	
public:
	template <typename... Args>
	Exception(const char* fmt, Args... args) {
		msg_ = format(fmt, args...).str();
	}

	template <typename... Args>
	Exception(const std::string& fmt, Args... args) {
		msg_ = format(fmt, args...).str();
	}
	
	~Exception() throw() {}
	
	virtual const char* what() const throw() {
		return msg_.c_str();
	}
};




struct DummySet {
	void clear() {}
	void insert(int x) {}
	void remove(int x) {}
	void insertLargest(int x) {}
	void removeLargest(int x) {}
};


template <class Set1, class Set2>
struct SetPair : public std::pair<Set1, Set2> {
	using std::pair<Set1, Set2>::first;
	using std::pair<Set1, Set2>::second;
	
	constexpr SetPair() :
		std::pair<Set1, Set2>()
	{}

	SetPair(const Set1& set1, const Set2& set2) :
		std::pair<Set1, Set2>(set1, set2)
	{}

	SetPair(Set1&& set1, Set2&& set2) :
		std::pair<Set1, Set2>(set1, set2)
	{}
	
	SetPair(const SetPair& other) :
		std::pair<Set1, Set2>(other)
	{}

	SetPair(SetPair&& other) :
		std::pair<Set1, Set2>(other)
	{}
	
	void clear() {
		first.clear;
		second.clear;
	}
	
	void insert(int x) {
		first.insert(x);
		second.insert(x);
	}
	
	void remove(int x) {
		first.remove(x);
		second.remove(x);
	}

	void insertLargest(int x) {
		first.insertLargest(x);
		second.insertLargest(x);
	}
	
	void removeLargest(int x) {
		first.removeLargest(x);
		second.removeLargest(x);
	}
};



/**
 * Square matrix of any type.
 */
template <class T>
class SquareMatrix {
private:
	int n_;
	T* data_;
	
public:
	SquareMatrix(int n) :
		n_(n),
		data_(new T[n * n])
	{}
	
	SquareMatrix(int n, const T& val) :
		//SquareMatrix(n)
		n_(n),
		data_(new T[n * n])
	{
		setAll(val);
	}

	SquareMatrix(const SquareMatrix& other) {
		n_ = other.n_;
		data_ = new T[n_ * n_];
		for (int i = 0; i < n_ * n_; ++i)
			data_[i] = other.data_[i];
	}
	
	SquareMatrix(SquareMatrix&& other) {
		n_ = other.n_;
		data_ = other.data_;
		other.n_ = 0;
		other.data_ = nullptr;
	}
	
	~SquareMatrix() {
		delete[] data_;
	}

	SquareMatrix& operator=(const SquareMatrix& other) {
		assert(other.n_ == n_);
		for (int i = 0; i < n_ * n_; ++i)
			data_[i] = other.data_[i];
		return *this;
	}

	SquareMatrix& operator=(SquareMatrix&& other) {
		assert(this != &other);
		assert(other.n_ == n_);
		delete[] data_;
		data_ = other.data_;
		other.n_ = 0;
		other.data_ = nullptr;
		return *this;
	}

	SquareMatrix& operator+=(const SquareMatrix& other) {
		assert(n_ == other.n_);
		for (int i = 0; i < n_ * n_; ++i)
			data_[i] += other.data_[i];
		return *this;
	}

	SquareMatrix& operator-=(const SquareMatrix& other) {
		assert(n_ == other.n_);
		for (int i = 0; i < n_ * n_; ++i)
			data_[i] -= other.data_[i];
		return *this;
	}

	SquareMatrix operator+(const SquareMatrix& other) const {
		assert(n_ == other.n_);
		SquareMatrix sum(*this);
		sum += other;
		return sum;
	}

	SquareMatrix operator-(const SquareMatrix& other) const {
		assert(n_ == other.n_);
		SquareMatrix sum(*this);
		sum -= other;
		return sum;
	}

	SquareMatrix& operator+=(const T& val) {
		for (int i = 0; i < n_ * n_; ++i)
			data_[i] += val;
		return *this;
	}

	SquareMatrix& operator-=(const T& val) {
		for (int i = 0; i < n_ * n_; ++i)
			data_[i] -= val;
		return *this;
	}

	SquareMatrix& operator*=(const T& val) {
		for (int i = 0; i < n_ * n_; ++i)
			data_[i] *= val;
		return *this;
	}

	SquareMatrix& operator/=(const T& val) {
		for (int i = 0; i < n_ * n_; ++i)
			data_[i] /= val;
		return *this;
	}

	SquareMatrix operator+(const T& val) const {
		SquareMatrix sum(*this);
		sum += val;
		return sum;
	}

	SquareMatrix operator-(const T& val) const {
		SquareMatrix sum(*this);
		sum -= val;
		return sum;
	}

	SquareMatrix operator*(const T& val) const {
		SquareMatrix sum(*this);
		sum *= val;
		return sum;
	}

	SquareMatrix operator/(const T& val) const {
		SquareMatrix sum(*this);
		sum /= val;
		return sum;
	}

	int getSize() const {
		return n_;
	}
	int getN() const {
		return n_;
	}
	
	void setAll(const T& value) {
		for (int i = 0; i < n_ * n_; ++i)
			data_[i] = value;
	}
	
	SquareMatrix& operator=(const T& value) {
		setAll(value);
	}
	
	T& operator() (int i, int j) {
		assert(0 <= i && i < n_);
		assert(0 <= j && j < n_);
		return data_[i + j * n_];
	}

	const T& operator() (int i, int j) const {
		assert(0 <= i && i < n_);
		assert(0 <= j && j < n_);
		return data_[i + j * n_];
	}
};

/**
 * (Floyd-Warshall)
 */
void adjMatToPredMat(const SquareMatrix<bool>& adjMat, SquareMatrix<bool>& predMat) {
	predMat = adjMat;
	for (int k = 0; k < predMat.getN(); ++k)
		for (int i = 0; i < predMat.getN(); ++i)
			for (int j = 0; j < predMat.getN(); ++j)
				predMat(i, j) |= predMat(i, k) & predMat(k, j);
}


#endif


