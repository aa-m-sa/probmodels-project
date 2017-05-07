/*
 *  BEANDisco: algorithms for counting the linear extensions of a DAG
 *  
 *  Copyright 2013-2015 Teppo Niinimäki <teppo.niinimaki(at)helsinki.fi>
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

#include <bitset>
#include <cstdint>
#include <queue>
#include <unordered_map>

#define BOOST_DYNAMIC_BITSET_DONT_USE_FRIENDS
#include <boost/dynamic_bitset.hpp>
#include <boost/functional/hash.hpp>

#include "common.hpp"


#ifndef LINEXT_HPP
#define LINEXT_HPP


template <typename T>
struct IdealInfo {
	T leCount;
	unsigned char* s;
};

//size_t countLinExtsExactFast_maxMapSize = 0;
//size_t countLinExtsExactFast_maxReserved = 0;
//size_t countLinExtsExactFast_numIdeals = 0;

template <typename T, typename Set>
T countLinExtsExactImpl(const SquareMatrix<bool>& predMat, const Set& emptySet) {
	int n = predMat.getSize();
	std::unordered_map<Set, IdealInfo<T> > infos;
	std::queue<Set> ideals;
	
	ideals.push(emptySet);
	infos[emptySet].s = new unsigned char[n];
	for (int j = 0; j < n; ++j) {
		infos[emptySet].s[j] = 0;
		for (int i = 0; i < n; ++i)
			if (predMat(i, j))
				++infos[emptySet].s[j];
	}
	infos[emptySet].leCount = 1;
	
	T count = 0;
	size_t numIdeals = 0;
	size_t maxSize = 0;
	size_t maxReserved = 0;
	Set x, y;
	while (!ideals.empty()) {
		++numIdeals;
		x = ideals.front(); ideals.pop();
		y = x;
		
		IdealInfo<T> xinfo = infos[x];
		count = xinfo.leCount;
		for (int i = 0; i < n; ++i) {
			if (y[i] || xinfo.s[i] > 0)
				continue;
			y.set(i);
			IdealInfo<T> yinfo;
			if (infos.count(y) == 0) {
				yinfo.s = new unsigned char[n];
				for (int j = 0; j < n; ++j)
					if (!y[j])
						yinfo.s[j] = xinfo.s[j] - predMat(i, j);
				yinfo.leCount = 0;
				ideals.push(y);
			} else {
				yinfo = infos[y];
			}
			yinfo.leCount += count;
			infos[y] = yinfo;
			y.reset(i);
		}
		delete[] xinfo.s;
		maxSize = std::max(maxSize, infos.size());
		maxReserved = std::max(maxReserved, infos.bucket_count());
		infos.erase(x);
	}
	//countLinExtsExactFast_maxMapSize = maxSize;
	//countLinExtsExactFast_maxReserved = maxReserved;
	//countLinExtsExactFast_numIdeals = numIdeals;
	return count;
}


template <class T>
class UIntSet {
private:
	T bits_;
	
public:
	UIntSet() {
	}
	
	UIntSet(T u) {
		bits_ = u;
	}
	
	bool operator== (const UIntSet& other) const {
		return bits_ == other.bits_;
	}
	
	bool operator[] (size_t i) const {
		return bits_ & ((T)1 << i);
	}
	
	void set(size_t i) {
		bits_ |= ((T)1 << i);
	}

	void reset(size_t i) {
		bits_ &= ~((T)1 << i);
	}
	
	friend struct std::hash<UIntSet<T> >;
};

namespace std {
	template <class T>
	struct hash<UIntSet<T> > {
	public:
		size_t operator() (const UIntSet<T>& set) const {
			std::hash<T> hash;
			return hash(set.bits_);
		}
	};
}

namespace std {
	template <>
	struct hash<boost::dynamic_bitset<> > {
	public:
		size_t operator() (const boost::dynamic_bitset<>& set) const {
			return boost::hash_value(set.m_bits);
		}
	};
}

/**
 * Counts the number of linear extensions of a partial order. Parameter: the predecessor matrix
 */
template <typename T>
T countLinExtsExact(const SquareMatrix<bool>& predMat) {
	int n = predMat.getSize();
	T count;
	if (n <= 64) {
		count = countLinExtsExactImpl<T, UIntSet<uint64_t> >(predMat, 0);
	}
	else if (n <= 128) {
		count = countLinExtsExactImpl<T, std::bitset<128> >(predMat, 0);
	}
	else if (n <= 256) {
		count = countLinExtsExactImpl<T, std::bitset<256> >(predMat, 0);
	}
	else {
		count = countLinExtsExactImpl<T, boost::dynamic_bitset<> >(predMat,
				boost::dynamic_bitset<>(n));
	}
	return count;
}

#endif

