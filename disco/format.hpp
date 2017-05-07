/*
 *  BEANDisco: formatting
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

#include <boost/format.hpp>
#include <string>

#ifndef FORMAT_HPP
#define FORMAT_HPP

boost::format format(boost::format fmt) {
	return fmt;
}

template <typename Arg, typename... Args>
boost::format format(boost::format fmt, Arg arg, Args... args) {
	return format(fmt % arg, args...);
}

template <typename... Args>
boost::format format(const char* fmt, Args... args) {
	return format(boost::format(fmt), args...);
}

template <typename... Args>
boost::format format(const std::string& fmt, Args... args) {
	return format(boost::format(fmt), args...);
}



class PrettyDuration {
private:
	double seconds_;
	friend std::ostream& operator<<(std::ostream& os, const PrettyDuration& prettyDur);
public:
	PrettyDuration(double seconds) : seconds_(seconds) {}
};

PrettyDuration prettyDuration(double seconds) {
	return PrettyDuration(seconds);
}

std::ostream& operator<<(std::ostream& os, const PrettyDuration& prettyDur) {
	double seconds = prettyDur.seconds_;
	double minutes = floor(seconds / 60);
	if (minutes > 0) {
		seconds -= 60 * minutes;
		double hours = floor(minutes / 60);
		if (hours > 0) {
			minutes -= 60 * hours;
			double days = floor(hours / 24);
			if (days > 0) {
				hours -= 24 * days;
				os << format("%.0f days ", days);
			}
			os << format("%.0f h ", hours);
		}
		os << format("%.0f min ", minutes);
	}
	if (prettyDur.seconds_ >= 10)
		os << format("%.1f s", seconds);
	else
		os << format("%.2f s", seconds);
	return os;
}



namespace units {
	struct None {
		static constexpr const char* str = "";
	};
	//const char* None::str = "";

	struct B {
		static constexpr const char* str = "B";
	};
	//const char* B::str = "";
}

template <typename Unit>
class PrettySize {
private:
	size_t size_;
	template <typename Unit2>
	friend std::ostream& operator<<(std::ostream& os, const PrettySize<Unit2>& prettySiz);
public:
	PrettySize(size_t size) : size_(size){}
};

template <typename Unit = units::None>
PrettySize<Unit> prettySize(size_t size) {
	return PrettySize<Unit>(size);
}

template <typename Unit>
std::ostream& operator<<(std::ostream& os, const PrettySize<Unit>& prettySiz) {
	double size = prettySiz.size_;
	if (size < 1024) {
		if (strlen(Unit::str) > 0)
			os << format("%d %s", size, Unit::str);
		else
			os << format("%d", size);
		return os;
	}
	size /= 1024;
	if (size < 1024) {
		os << format("%.3g k%s", size, Unit::str);
		return os;
	}
	size /= 1024;
	if (size < 1024) {
		os << format("%.3g M%s", size, Unit::str);
		return os;
	}
	size /= 1024;
	if (size < 1024) {
		os << format("%.3g G%s", size, Unit::str);
		return os;
	}
	size /= 1024;
	os << format("%.3g T%s", size, Unit::str);
	return os;
}


#endif

