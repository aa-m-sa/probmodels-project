/*
 *  BEANDisco: debug
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


#include <memory>
#include <string>
#include <cxxabi.h>


#ifndef DEBUG_HPP
#define DEBUG_HPP


const char* demangle(const char* name) {
	int status = -4;
	const char* demangledName = abi::__cxa_demangle(name, nullptr, nullptr, &status);
	return (status == 0) ? demangledName : name;
}

template <class T>
const char* typeName(const T& t) {
	return demangle(typeid(t).name());
}
template <class T>
const char* typeName() {
	return demangle(typeid(T).name());
}
const char* name(const std::type_info& i) {
	return demangle(i.name());
}


template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec) {
	//os << "(std::vector<" << typeName<T>() << ">) ";
	os << "[";
	if (!vec.empty()) {
		os << vec[0];
		for (int i = 1; i < vec.size(); ++i)
			os << "," << vec[i];
	}
	os << "]";
	return os;
}


/*std::string demangle(const char* name) {
	int status = -4;
	std::unique_ptr<char, void(*)(void*)> res {
		abi::__cxa_demangle(name, NULL, NULL, &status),
		std::free
	};
	return (status==0) ? res.get() : name ;
}

template <class T>
std::string typeName(const T& t) {
	return demangle(typeid(t).name());
}
std::string name(const std::type_info& i) {
	return demangle(i.name());
}*/

#endif

