/*
 *  BEANDisco: data class for reading and writing data
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

#include <cstddef>
#include <exception>
#include <boost/format.hpp>
#include <iostream>
#include <sstream>
#include <string>
#include <list>

#include "common.hpp"

#ifndef DATA_HPP
#define DATA_HPP

using Datum = unsigned char;


/*class DataReadException : public Exception {
public:
	DataReadException(const string& msg) : Exception(msg) {
	}
};*/


class DataView {
public:
	virtual ~DataView() {};
	virtual int getNumSamples() const = 0;
	virtual int getNumVariables() const = 0;
	virtual int getArity(int i) const = 0;
	virtual int* getCounts(const std::vector<int>& vars) const = 0;
};


class Data : public DataView {
private:
	Data(const Data&); // disable copying
	Data& operator=(const Data&); // disable copying
	
	void computeArities() {
		assert(arities == nullptr);
		arities = (int*)malloc(sizeof(int) * nVariables);
		for (int v = 0; v < nVariables; ++v) {
			int arity = 0;
			for (int i = 0; i < nSamples; ++i) {
				if ((*this)(v,i) >= arity)
					arity = (*this)(v,i) + 1;
			}
			arities[v] = arity;
		}
	}
public:
	int nVariables;
	int nSamples;
	Datum* data;
	int* arities;
	
	Data() {
		nVariables = 0;
		nSamples = 0;
		data = nullptr;
		arities = nullptr;
	}
	
	Data(Data&& other) {
		nVariables = other.nVariables;
		nSamples = other.nSamples;
		data = other.data;
		arities = other.arities;
		other.nVariables = 0;
		other.nSamples = 0;
		other.data = nullptr;
		other.arities = nullptr;
	}
	
	void clear() {
		nVariables = 0;
		nSamples = 0;
		if (data)
			free(data);
		data = nullptr;
		if (arities)
			free(arities);
		arities = nullptr;
	}
	
	~Data() {
		clear();
	}
	
	int getNumSamples() const {
		return nSamples;
	}
	
	int getNumVariables() const {
		return nVariables;
	}
	
	Datum& operator()(int v, int i) {
		//return data[v * nSamples + i];
		return data[i * nVariables + v];
	}
	
	Datum operator()(int v, int i) const {
		//return data[v * nSamples + i];
		return data[i * nVariables + v];
	}
	
	void selectFrom(const Data& other, std::list<int> variables, std::list<int> samples) {
		nVariables = variables.size();
		nSamples = samples.size();
		data = (Datum*)malloc(sizeof(Datum) * nVariables * nSamples);
		int v = 0;
		for (int v2 : variables) {
			int i = 0;
			for (int i2 : samples) {
				(*this)(v, i) = other(v2, i2);
				++i;
			}
			++v;
		}
		computeArities();
	}
	
	void read(std::istream& file, int nVars, int nSamps) {
		nVariables = nVars;
		nSamples = nSamps;
		data = (Datum*)malloc(sizeof(Datum) * nVariables * nSamples);
		
		// load the data
		std::string row;
		for (int i = 0; i < nSamples; ++i) {
			if (file.eof())
				throw Exception("Not enough rows (%d while %d expected).", i, nSamples);
			getline(file, row);
			std::istringstream rowStream(row);
			for (int v = 0; v < nVariables; ++v) {
				int tmp;
				rowStream >> tmp;
				if (rowStream.fail())
					throw Exception("Could not read value on row %d column %d", i + 1, v + 1);
				(*this)(v, i) = (Datum) tmp;
			}
		}
		
		computeArities();
	}

	void read(std::istream& file) {
		nVariables = 1;
		nSamples = 1;
		data = (Datum*)malloc(sizeof(Datum) * nVariables * nSamples);
		
		std::string row;
		
		// read the first row (and get the number of variables)
		getline(file, row);
		std::istringstream rowStream(row);
		int v = 0;
		rowStream >> std::ws;
		while(!rowStream.eof()) {
			if (v >= nVariables) {
				nVariables *= 2;
				data = (Datum*)realloc(data, sizeof(Datum) * nVariables * nSamples);
			}
			int tmp;
			rowStream >> tmp;
			if (rowStream.fail())
				throw Exception("Could not read value on row 1 column %d.", v + 1);
			data[v] = (Datum) tmp;
			rowStream >> std::ws;
			++v;
		}
		file >> std::ws;
		nVariables = v;
		data = (Datum*)realloc(data, sizeof(Datum) * nVariables * nSamples);
		
		// load the rest of the data
		int i = 1;
		while(!file.eof()) {
			if (i >= nSamples) {
				nSamples *= 2;
				data = (Datum*)realloc(data, sizeof(Datum) * nVariables * nSamples);
			}
			getline(file, row);
			std::istringstream rowStream(row);
			for (int v = 0; v < nVariables; ++v) {
				int tmp;
				rowStream >> tmp;
				if (rowStream.fail())
					throw Exception("Could not read %dth value on row %d", v + 1, i + 1);
				data[i * nVariables + v] = (Datum) tmp;
			}
			file >> std::ws;
			++i;
		}
		nSamples = i;
		data = (Datum*)realloc(data, sizeof(Datum) * nVariables * nSamples);
		
		computeArities();
	}
	
	int getArity(int v) const {
		return arities[v];
	}

	int* getCounts(const std::vector<int>& vars) const {
		// create count table and initialize to zero
		int nValues = 1;
		for (size_t i = 0; i < vars.size(); ++i)
			nValues *= getArity(vars[i]);
		int* counts = new int[nValues];
		for (int i = 0; i < nValues; ++i)
			counts[i] = 0;
		
		// fill counts
		for (int j = 0; j < nSamples; ++j) {
			int index = 0;
			for (size_t i = 0; i < vars.size(); ++i)
				index = index * getArity(vars[i]) + (*this)(vars[i], j);
			++counts[index];
		}
		return counts;
	}
};


class DataColumns : public DataView {
private:
	const Data& data_;
	std::vector<int> variables_;
public:
	DataColumns(const Data& data, const std::vector<int>& vars)
			: data_(data), variables_(vars) {}

	DataColumns(const DataColumns& dataCols)
			: data_(dataCols.data_), variables_(dataCols.variables_) {}

	DataColumns(const Data& data)
			: data_(data), variables_(data.getNumVariables()) {
		for (size_t i = 0; i < variables_.size(); ++i)
			variables_[i] = i;
	}
	
	int getNumSamples() const {
		return data_.getNumSamples();
	}
	
	int getNumVariables() const {
		return variables_.size();
	}
	
	int getArity(int v) const {
		return data_.getArity(variables_[v]);
	}
	
	Datum operator()(int v, int i) const {
		return data_(variables_[v], i);
	}
	
	int* getCounts(const std::vector<int>& vars) const {
		std::vector<int> transVars(vars.size());
		for (size_t i = 0; i < vars.size(); ++i)
			transVars[i] = variables_[vars[i]];
		return data_.getCounts(transVars);
	}
};

#endif

