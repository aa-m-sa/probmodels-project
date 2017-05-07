/*
 *  BEANDisco: input and output streams (for files and stdio)
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
#include <fstream>
#include <boost/iostreams/stream_buffer.hpp>
#include <boost/iostreams/device/file_descriptor.hpp> 
//#include <boost/iostreams/device/null.hpp>
#include <unistd.h>

#ifndef IO_H
#define IO_H


class InStream {
public:
	virtual ~InStream() {};
	virtual std::istream& stream() = 0;
	virtual bool isOpen() = 0;
	virtual bool hasInput() = 0;
};


class InFile : public InStream {
private:
	std::ifstream in_;
public: 
	InFile(std::string name) :
		in_(name)
	{
	}
	
	std::istream& stream() {
		return in_;
	}
	
	bool isOpen() {	
		return in_.is_open();
	}
	
	bool hasInput() {
		return true;
	}
};

class StdIn : public InStream {
private:
public: 
	StdIn()
	{
	}
	
	std::istream& stream() {
		return std::cin;
	}
	
	bool isOpen() {	
		return true;
	}
	
	bool hasInput() {
		return true;
	}
};



class OutStream {
public:
	virtual ~OutStream() {};
	virtual std::ostream& stream() = 0;
	virtual void newField() = 0;
	virtual void endField() = 0;
	virtual bool isOpen() = 0;
	virtual bool hasOutput() = 0;

	template <typename T>
	void writeField(const T& value) {
		newField();
		stream() << value;
		endField();
	}

	template <typename T>
	void writelnField(const T& value) {
		newField();
		stream() << value << "\n";
		endField();
	}
};

class NullOutStream : public OutStream {
private:
	std::ostream out_;
public: 
	NullOutStream() :
		out_(0)
	{
	}
	
	std::ostream& stream() {
		return out_;
	}
	
	void newField() {
	}
	
	void endField() {
	}
	
	bool isOpen() {	
		return true;
	}
	
	bool hasOutput() {
		return false;
	}
};


class CumulativeOutFile : public OutStream {
private:
	std::ofstream out_;
public: 
	CumulativeOutFile(std::string name) :
		out_(name)
	{
	}
	
	std::ostream& stream() {
		return out_;
	}
	
	void newField() {
	}
	
	void endField() {
		out_.flush();
	}
	
	bool isOpen() {	
		return out_.is_open();
	}
	
	bool hasOutput() {
		return true;
	}
};

class CumulativeStdOut : public OutStream {
private:
public: 
	CumulativeStdOut()
	{
	}
	
	std::ostream& stream() {
		return std::cout;
	}
	
	void newField() {
	}
	
	void endField() {
		std::cout.flush();
	}
	
	bool isOpen() {	
		return true;
	}
	
	bool hasOutput() {
		return true;
	}
};

class OverwriteOutFile : public OutStream {
private:
	boost::iostreams::file_descriptor_sink fd_;
	boost::iostreams::stream_buffer<boost::iostreams::file_descriptor_sink> buf_;
	std::ostream out_;

public:
	OverwriteOutFile(std::string name) :
		fd_(name), buf_(fd_), out_(&buf_)
	{
	}
	
	std::ostream& stream() {
		return out_;
	}
	
	void newField() {
		out_.seekp(0);
		//out_.flush();
		//ftruncate(fd_.handle(), 0);
		// TODO (onko edes tarpeellinen?)
		int err = ftruncate(fd_.handle(), 0);
		if (err)
			throw std::ostream::failure("Error: failed to truncate file");
	}
	
	void endField() {
		out_.flush();
	}
	
	bool isOpen() {	
		return fd_.is_open();
	}
	
	bool hasOutput() {
		return true;
	}
};



std::unique_ptr<OutStream> openOutStream(const std::string& filename, bool overwrite = false) {
	OutStream* os = nullptr;
	if (filename == "-") {
		//assert(!overwrite);
		if (overwrite)
			throw Exception("Can't open stdout in update mode.");
		os = new CumulativeStdOut();
	}
	else if (filename == "") {
		os = new NullOutStream();
	}
	else {
		if (!overwrite)
			os = new CumulativeOutFile(filename);
		else
			os = new OverwriteOutFile(filename);
	}
	if (!os->isOpen()) {
		delete os;
		throw Exception("Couldn't open file '%s' for writing.", filename);
	}
	return std::unique_ptr<OutStream>(os);
}


std::unique_ptr<InStream> openInStream(const std::string& filename) {
	InStream* is = nullptr;
	if (filename == "-") {
		is = new StdIn();
	}
	else if (filename == "") {
		//is = new NullInStream();
		is = nullptr;
	}
	else {
		is = new InFile(filename);
	}
	if (is != nullptr && !is->isOpen()) {
		if (is) delete is;
		throw Exception("Couldn't open file '%s' for reading.", filename);
	}
	return std::unique_ptr<InStream>(is);
}


#endif

