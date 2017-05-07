/*
 *  BEANDisco: linear order definition
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

#ifndef ORDER_HPP
#define ORDER_HPP

class OrderFamily {
private:
	using Family = OrderFamily;
	
public:
	const int n;

	OrderFamily(int _n) :
			n(_n)
	{}
	
	bool operator==(const OrderFamily& other) const {
		return n == other.n;
	}
	
	double logSize() const {
		return lgamma(n + 1);
	}
	
	class Instance {
	public:
		const Family& family;
	private:
		//int* order_;
		//int* invOrder_;
		std::unique_ptr<int[]> order_;
		std::unique_ptr<int[]> invOrder_;
		
	public:
		
		/*Instance(const Family& _family) : family(_family) {
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
		
		Instance(const Instance& other) : family(other.family) {
			order_ = new int[family.n];
			invOrder_ = new int[family.n];
			for (int i = 0; i < family.n; ++i) {
				order_[i] = other.order_[i];
				invOrder_[i] = other.invOrder_[i];
			}
		}
		
		Instance& operator=(const Instance& other) {
			assert(family == other.family);
			for (int i = 0; i < family.n; ++i) {
				order_[i] = other.order_[i];
				invOrder_[i] = other.invOrder_[i];
			}
			return *this;
		}

		Instance(Instance&& other) : family(other.family) {
			order_ = other.order_;
			invOrder_ = other.invOrder_;
			other.order_ = nullptr;
			other.invOrder_ = nullptr;
		}
		
		Instance& operator=(Instance&& other) {
			assert(family == other.family);
			order_ = other.order_;
			invOrder_ = other.invOrder_;
			other.order_ = nullptr;
			other.invOrder_ = nullptr;
			return *this;
		}*/

		Instance(const Family&& _family) = delete;
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
		
		Instance(const Family&& _family, const std::vector<int>& order) = delete;
		Instance(const Family& _family, const std::vector<int>& order) :
			family(_family),
			order_(new int[family.n]),
			invOrder_(new int[family.n])
		{
			assert(family.n == order.size());
			for (int i = 0; i < family.n; ++i) {
				assert(0 <= order[i] && order[i] < family.n);
				order_[i] = order_[i];
				invOrder_[order_[i]] = i;
			}
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
		}
		
		int getIndex(int elem) const {
			assert(0 <= elem && elem < family.n);
			return invOrder_[elem];
		}
		
		template <typename Subset>
		void getPreds(int elem, Subset& preds) const {
			int i = getIndex(elem);
			preds.clear(); // ????
			for (int j = 0; j < i; ++j)
				preds.insert(order_[j]);
		}
		
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

		void writeTo(std::ostream& os) const {
			for (int i = 0; i < family.n; ++i) {
				if (i > 0)
					os << " ";
				os << order_[i];
			}
		}
	
		void readFrom(std::istream& is) {
			for (int i = 0; i < family.n; ++i) {
				is >> order_[i];
				invOrder_[order_[i]] = i;
			}
		}
	};

	class UniformRandSwap : public MCProposalDist<Family> {
	private:
		const Family& family_;
	public:
		UniformRandSwap(const Family& family) :
			family_(family)
		{}

		void randStep(Instance& instance) const {
			int i = randuint(family_.n);
			int j = (i + randuint(family_.n - 1)) % family_.n;
			instance.swapElemsAt(i, j);
		}
	};
};


std::ostream& operator<<(std::ostream& os, const OrderFamily::Instance& instance) {
	instance.writeTo(os);
	return os;
}

std::istream& operator>>(std::istream& is, OrderFamily::Instance& instance) {
	instance.readFrom(is);
	return is;
}


#endif
