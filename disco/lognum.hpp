/*
 *  BEANDisco: class for handling number in logarithmic form
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

#include <limits>

#include <cassert>


#ifndef LOGNUM_HPP
#define LOGNUM_HPP



template <class T>
struct Lognum {
private:
	T logx_;
	
public:
	
	Lognum() {
		logx_ = -std::numeric_limits<T>::infinity();
	}
	
	template <class U>
	Lognum(U x) {
		logx_ = log(x);
	}

	Lognum operator=(Lognum x) {
		logx_ = x.logx_;
		return *this;
	}
	
	template <class U>
	Lognum operator=(U x) {
		logx_ = log(x);
		return *this;
	}
	
	void setLog(T e) {
		logx_ = e;
	}
	
	T getLog() const {
		return logx_;
	}
	
	T getVal() const {
		return exp(logx_);
	}

	template <typename F>
	explicit operator F() const {
		return (F) getVal();
	}
	
	/*template <typename F>
	operator F() const {
		static_assert(std::is_floating_point<F>::value,
			"implicit conversion from Lognum only possible to floating point types");
		return (F) getVal();
	}*/
	
	//template <typename FL>
	//operator Lognum<FL>() {
	//	Lognum<FL> tmp;
	//	tmp.setLog(getLog());
	//	return tmp;
	//}
	
	Lognum operator+(Lognum b) const {
		Lognum res;
		if (logx_ == -std::numeric_limits<T>::infinity())
			res = b;
		else if (b.logx_ == -std::numeric_limits<T>::infinity())
			res = *this;
		else
			res.logx_ = (logx_ > b.logx_) ?
					logx_ + log1p(exp(b.logx_ - logx_)) :
					b.logx_ + log1p(exp(logx_ - b.logx_));
		return res;
	}
	
	Lognum operator+=(Lognum b) {
		*this = *this + b;
		return *this;
	}
	
	Lognum operator-(Lognum b) const {
		//assert(logx_ >= b.logx_);
		Lognum res;
		if (logx_ == -std::numeric_limits<T>::infinity())
			res.logx_ = -std::numeric_limits<T>::infinity();
		else if (logx_ < b.logx_)
			res.logx_ = -std::numeric_limits<T>::infinity();
		else
			res.logx_ = logx_ + log1p(-exp(b.logx_ - logx_));
		return res;
	}
	
	Lognum operator-=(Lognum b) {
		*this = *this - b;
		return *this;
	}
	
	Lognum operator*(Lognum b) const {
		Lognum res;
		res.logx_ = logx_ + b.logx_;
		return res;
	}
	
	Lognum operator*=(Lognum b) {
		logx_ += b.logx_;
		return *this;
	}
	
	Lognum operator/(Lognum b) const {
		Lognum res;
		res.logx_ = logx_ - b.logx_;
		return res;
	}
	
	Lognum operator/=(Lognum b) {
		logx_ -= b.logx_;
		return *this;
	}

	bool operator==(Lognum b) const {
		return logx_ == b.logx_;
	}

	bool operator!=(Lognum b) const {
		return logx_ != b.logx_;
	}

	bool operator<(Lognum b) const {
		return logx_ < b.logx_;
	}

	bool operator>(Lognum b) const {
		return logx_ > b.logx_;
	}

	bool operator<=(Lognum b) const {
		return logx_ <= b.logx_;
	}

	bool operator>=(Lognum b) const {
		return logx_ >= b.logx_;
	}
	
	Lognum pow(T e) {
		Lognum res;
		res.logx_ = logx_ * e;
		return res;
	}
};

template <typename T>
Lognum<T> lognumFromLog(T l) {
	Lognum<T> x; x.setLog(l); return x;
}


template <typename T>
Lognum<T> pow(Lognum<T> x, T e){
	return x.pow(e);
}


template <typename T>
T log(Lognum<T> x){
	return x.getLog();
}

template <typename T, typename F>
void setLog(Lognum<T>& x, F l){
	return x.setLog(l);
}




template <typename T>
Lognum<T> gamma(Lognum<T> x){
	Lognum<T> res;
	res.setLog(lgamma(double(x)));
	return res;
}



template <typename T, typename F>
T to(F x){
	return T(x);
}


template <class T>
std::ostream& operator<<(std::ostream& os, const Lognum<T>& x) {
	os << "e^(" << x.getLog() << ")";
	return os;
}

/*template <class T>
std::istream& operator>>(std::istream& is, Lognum<T>& x) {
	std::string tmp;
	is >> tmp;
	assert(tmp == "≺");
	is >> "e^(" >> x.getLog() >> ")";
	return is;
}*/



namespace std {
	template <typename T> class numeric_limits<Lognum<T>> {
	public:
		static constexpr bool is_specialized = true;
		static constexpr Lognum<T> min() noexcept {
			return lognumFromLog(std::numeric_limits<T>::lowest());
		}
		static constexpr Lognum<T> max() noexcept {
			return lognumFromLog(std::numeric_limits<T>::max());
		}
		static constexpr Lognum<T> lowest() noexcept {
			return Lognum<T>(0);
		}
		static constexpr int  digits = 0;    // ???
		static constexpr int  digits10 = 0;  // ???
		static constexpr bool is_signed = false;
		static constexpr bool is_integer = false;
		static constexpr bool is_exact = false;
		static constexpr int radix = 0;      // ???
		static constexpr Lognum<T> epsilon() noexcept {
			return lognumFromLog(std::numeric_limits<T>::epsilon());
		}
		static constexpr Lognum<T> round_error() noexcept {
			return Lognum<T>(0.5); // ?
		}

		static constexpr int  min_exponent = 0;     // ???
		static constexpr int  min_exponent10 = 0;   // ???
		static constexpr int  max_exponent = 0;     // ???
		static constexpr int  max_exponent10 = 0;   // ???

		static constexpr bool has_infinity = true;
		static constexpr bool has_quiet_NaN = true;
		static constexpr bool has_signaling_NaN = true;
		static constexpr float_denorm_style has_denorm = denorm_absent; // ?
		static constexpr bool has_denorm_loss = false;
		static constexpr Lognum<T> infinity() noexcept {
			return Lognum<T>(std::numeric_limits<T>::infinity());
		}
		static constexpr Lognum<T> quiet_NaN() noexcept {
			return Lognum<T>(std::numeric_limits<T>::quiet_NaN());
		}
		static constexpr Lognum<T> signaling_NaN() noexcept {
			return Lognum<T>(std::numeric_limits<T>::signaling_NaN());
		}
		static constexpr Lognum<T> denorm_min() noexcept {
			return lognumFromLog(std::numeric_limits<T>::lowest());
		}

		static constexpr bool is_iec559 = true;
		static constexpr bool is_bounded = true;
		static constexpr bool is_modulo = false;

		static constexpr bool traps = false;    // ??
		static constexpr bool tinyness_before = false; // ???
		static constexpr float_round_style round_style = round_toward_zero; // ???
	};
}



//template <typename T, typename FL>
//T to(Lognum<FL> x){
//	return to<T>(x.getVal());
//}


//template <typename T>
//T to(Lognum<T> x){
//	return x.getVal();
//}

//template <typename T, typename FL>
//T to(Lognum<FL> x){
//	return T(x.getVal());
//}

//template <typename TL, typename FL>
//Lognum<TL> to<Lognum<TL> >(Lognum<FL> x){
//	return Lognum<TL>(x);
//}



//template <class U>
//U to(double x){
//	return U(x);
//}
//
//template <class U>
//U to(Lognum<double> x){
//	return U(x);
//}
//
//template <>
//double to<double>(Lognum<double> x){
//	return x.getVal();
//}


#endif


