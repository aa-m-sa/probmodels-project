/*
 *  BEANDisco: structural features
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

#include <iostream>

#include "common.hpp"
//#include "stacksubset.hpp"
#include "dag.hpp"

#ifndef FEATURE_HPP
#define FEATURE_HPP



/**
 * A boolean feature that is always true.
 */
template <typename T> 
class TrueFeature {
public:
	using ValueType = T;
	using Subset = DummySet;

	bool holds(int v, const Subset& pa) const {
		return true;
	}

	//bool holds(int v, const StackSubset& pa) const {
	//	return true;
	//}
};

/**
 * A boolean feature describing the presence of an arc.
 */
template <typename T>
class ArcProbFeature {
private:
	Arc arc_;
	
public:
	using ValueType = T;

	ArcProbFeature(const Arc& arc)
		: arc_(arc) {
	}
	
	T zeroValue() const {
		return T(0);
	}

	struct Subset {
	private:
		const int tail_;
		bool tailPresent_;
	
	public:
		Subset(const ArcProbFeature& feature) :
			tail_(feature.arc_.tail)
		{
			clear();
		}
	
		void clear() {
			tailPresent_ = false;
		}
	
		void insert(int v) {
			assert(!(v == tail_) || !tailPresent_);
			tailPresent_ ^= (v == tail_);
		}
	
		void remove(int v) {
			assert(!(v == tail_) || tailPresent_);
			tailPresent_ ^= (v == tail_);
		}

		void insertLargest(int v) {
			insert(v);
		}
		void removeLargest(int v) {
			remove(v);
		}
		
		friend class ArcProbFeature;
	};

	bool holds(int v, const Subset& pa) const {
		return (arc_.head != v || pa.tailPresent_);
	}

	//bool holds(int v, const StackSubset& pa) const {
	//	return (arc_.head != v || pa.contains(arc_.tail));
	//}
	
	T value(int v, const Subset& pa) const {
		return holds(v, pa) ? 1.0 : 0.0;
	}
	
	//T value(int v, const StackSubset& pa) const {
	//	return holds(v, pa) ? 1.0 : 0.0;
	//}

	T value (const DagFamily::Instance& dag) const {
		return dag.contains(arc_) ? 1.0 : 0.0;
	}
};



/**
 * A templated map data structure with Arc as index type.
 */
/*template <class T>
class ArcMap {
private:
	SquareMatrix<T> data_;
	//ArcMap(const ArcMap&) = delete; // disable copy constructor
	//ArcMap& operator=(const ArcMap&) = delete; // disable copying
public:
	ArcMap(int nNodes) :
		data_(nNodes)
	{
	}

	ArcMap(int nNodes, const T& initValue) :
		data_(nNodes, initValue)
	{
	}
	
	void setAll(T value) {
		data_.setAll(value);
	}
	
	T& operator[] (Arc arc) {
		return data_(arc.head, arc.tail);
	}

	const T& operator[] (Arc arc) const {
		return data_(arc.head, arc.tail);
	}

	const ArcMap& operator+=(const ArcMap& other) {
		data_ += other.data;
	}

	ArcMap operator+(const ArcMap& other) {

		data_ += other.data;
	}
};*/

template <class T>
class ArcMap : public SquareMatrix<T> {
public:
	ArcMap(int nNodes) :
		SquareMatrix<T>(nNodes)
	{}

	ArcMap(int nNodes, const T& initValue) :
		SquareMatrix<T>(nNodes, initValue)
	{}

	/*ArcMap(const ArcMap& other) :
		SquareMatrix<T>(other)
	{}*/

	ArcMap(const SquareMatrix<T>& other) :
		SquareMatrix<T>(other)
	{}

	ArcMap(SquareMatrix<T>&& other) :
		SquareMatrix<T>(other)
	{}

	/*ArcMap& operator=(const ArcMap& other) {
		SquareMatrix<T>(other)
	}*/

	using SquareMatrix<T>::setAll;
	
	T& operator[] (Arc arc) {
		return SquareMatrix<T>::operator()(arc.head, arc.tail);
	}

	const T& operator[] (Arc arc) const {
		return SquareMatrix<T>::operator()(arc.head, arc.tail);
	}
};


/*template <typename T>
class ArcProbsFeature {
private:
	int n_;
public:
	typedef ArcMap<T> ValueType;

	ArcProbsFeature(int n)
		: n_(n)
	{}
	
	T zeroValue() const {
		return ArcMap<T>(n_);
	}

	T value (const DagFamily::Instance& dag) const {
		//return dag.contains(arc_) ? 1.0 : 0.0;
	}
};*/


#endif
