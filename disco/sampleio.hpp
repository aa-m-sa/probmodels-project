/*
 *  BEANDisco: sample input and output
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

#include "sample.hpp"
#include "io.hpp"
#include "dist/dist.hpp"

#ifndef SAMPLEIO_HPP
#define SAMPLEIO_HPP



template <typename SampleSpace, typename Model, typename T>
class SampleFileSource : public SampleSource<SampleSpace, Model, T> {
private:
	using SampleSource<SampleSpace, Model, T>::sampleSpace_;
	using SampleSource<SampleSpace, Model, T>::model_;
	//std::istream& in_;
	InStream* in_;

	void readSample(typename SampleSpace::Instance& instance, T& weight) {
		in_->stream() >> instance;
		double tmp;
		in_->stream() >> tmp;
		setLog(weight, tmp);
		in_->stream() >> std::ws;
	}
public:
	//SampleFileSource(const SampleSpace& sampleSpace, const Model* model, std::istream& in)
	SampleFileSource(const SampleSpace& sampleSpace, const Model* model, InStream* in)
		: SampleSource<SampleSpace, Model, T>(sampleSpace, model, 0, 0), in_(in)
	{
	}
	
	std::unique_ptr<Sample<SampleSpace, Model, T>> getSample() {
		typename SampleSpace::Instance instance(sampleSpace_);
		T weight;
		readSample(instance, weight);
		auto modelSample = std::make_shared<const ModelSample<SampleSpace, Model, T>>(model_, instance);
		return std::unique_ptr<Sample<SampleSpace, Model, T>>(
				new Sample<SampleSpace, Model, T>(modelSample, weight));
	}

	bool eof() const {
		return in_->stream().eof();
	}

	bool skip(size_t n) {
		typename SampleSpace::Instance instance(sampleSpace_);
		T weight;
		for (size_t i = 0; i < n; ++i) {
			if (eof())
				return false;
			readSample(instance, weight);
		}
		return true;
	}
};

template <typename SampleSpace, typename Model, typename T>
class SampleDistSource : public SampleSource<SampleSpace, Model, T> {
private:
	using SampleSource<SampleSpace, Model, T>::sampleSpace_;
	std::shared_ptr<ModelDist<SampleSpace, Model, T>> dist_;
	mutable std::unique_ptr<Sample<SampleSpace, Model, T>> sample_;
public:
	SampleDistSource(std::shared_ptr<ModelDist<SampleSpace, Model, T>> dist) :
		SampleSource<SampleSpace, Model, T>(dist->getSampleSpace(), dist->getModel(), dist->nLevels(), dist->nHarmonicLevels()),
		dist_(dist),
		sample_(nullptr)
	{
	}

	std::shared_ptr<ModelDist<SampleSpace, Model, T>> getDist() {
		return dist_;
	}
	
	std::unique_ptr<Sample<SampleSpace, Model, T>> getSample() {
		assert(sample_);
		return std::move(sample_);
	}

	bool eof() const {
		sample_ = dist_->rand();
		return !sample_;
	}

	void setInterruptCallback(std::function<bool()> interruptCallback) {
		dist_->setInterruptCallback(interruptCallback);
	}
};


template <typename SubSampleSpace, typename SubModel, typename T>
class SubSampleSource {
public:
	virtual std::unique_ptr<SampleSource<SubSampleSpace, SubModel, T>> newSubSource(const SubModel* subModel) = 0;
	virtual ~SubSampleSource() {}
};

template <typename SubSampleSpace, typename SubModel, typename SubDist, typename T>
class SubSampleDistSource : public SubSampleSource<SubSampleSpace, SubModel, T> {
private:
	const SubSampleSpace& subSampleSpace_;
public:
	SubSampleDistSource(const SubSampleSpace& subSampleSpace) :
		subSampleSpace_(subSampleSpace)
	{}
	
	std::unique_ptr<SampleSource<SubSampleSpace, SubModel, T>> newSubSource(const SubModel* subModel) {
		return std::unique_ptr<SampleSource<SubSampleSpace, SubModel, T>>(
				new SampleDistSource<SubSampleSpace, SubModel, T>(
					std::unique_ptr<SubDist>(new SubDist(subSampleSpace_, subModel))));
	}
};



template <typename SampleSpace, typename Model, typename T>
class SampleFileSink : public SampleSink<SampleSpace, Model, T> {
private:
	OutStream* out_;
public:
	SampleFileSink(OutStream* out)
		: out_(out) {
	}
	
	void addSample(const Sample<SampleSpace, Model, T>* sample) {
		out_->newField();
		//out_->stream() << (*sample) << std::endl;
		out_->stream() << sample->getInstance();
		out_->stream() << " ";
		out_->stream() << std::setiosflags(std::ios::scientific);
		out_->stream() << std::setprecision(10);
		out_->stream() << log(sample->getWeight());
		out_->stream() << std::endl;
		out_->endField();
	}
};

template <typename SampleSpace, typename Model, typename T>
class SubSampleCountFileSink : public SampleSink<SampleSpace, Model, T> {
private:
	OutStream* out_;
public:
	SubSampleCountFileSink(OutStream* out)
		: out_(out) {
	}
	
	void addSample(const Sample<SampleSpace, Model, T>* sample) {
		out_->newField();
		out_->stream() << sample->getNumSubSamples();
		out_->stream() << std::endl;
		out_->endField();
	}
};

template <typename SampleSpace, typename Model, typename T>
class CallbackFileSink : public SampleSink<SampleSpace, Model, T> {
private:
	using Callback = std::function<void(OutStream*, const Sample<SampleSpace, Model, T>*)>;
	OutStream* out_;
	Callback callback_;
	size_t spacing_;
	size_t n_;
public:
	CallbackFileSink(OutStream* out, Callback callback, int spacing)
		: out_(out), callback_(callback), spacing_(spacing), n_(0)
	{
		assert(spacing > 0);
	}
	
	void addSample(const Sample<SampleSpace, Model, T>* sample) {
		++n_;
		if (n_ % spacing_ == 0) {
			out_->newField();
			callback_(out_, sample);
			out_->endField();
		}
	}
};

template <typename SampleSpace, typename Model, typename T>
std::unique_ptr<CallbackFileSink<SampleSpace, Model, T>> makeFileSink(
		OutStream* out,
		std::function<void(OutStream*, const Sample<SampleSpace, Model, T>*)> callback,
		int spacing = 1) {
	return std::unique_ptr<CallbackFileSink<SampleSpace, Model, T>>(new CallbackFileSink<SampleSpace, Model, T>(out, callback, spacing));
}

#endif
