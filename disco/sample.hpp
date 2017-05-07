/*
 *  BEANDisco: sample classes and abstract SampleSource and SampleSink
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

#include "model.hpp"

#ifndef SAMPLE_HPP
#define SAMPLE_HPP


template <typename SampleSpace, typename Model, typename T>
class ModelSample {
private:
	using RestModel = RestrictedModel<Model, SampleSpace>;
	using Instance = typename SampleSpace::Instance;
	Instance instance_;
	RestModel restrictedModel_;

public:
	ModelSample(const Model* model, const Instance& instance) :
		instance_(instance),
		restrictedModel_(model, instance_)
	{}

	const Instance& getInstance() const {
		return instance_;
	}
	
	const RestModel& getRestrictedModel() const {
		return restrictedModel_;
	}

	T getUnnormProb() const {
		return restrictedModel_.getUnnormProb();
	}
};



template <typename SampleSpace, typename Model, typename T>
class Sample {
private:
	using RestModel = RestrictedModel<Model, SampleSpace>;
	using Instance = typename SampleSpace::Instance;
	std::shared_ptr<const ModelSample<SampleSpace, Model, T>> modelSample_;
	T weight_;
	const std::vector<T> normConstRatioWeights_;
	const std::vector<T> invNormConstRatioWeights_;
	long long int nSubSamples_;
public:
	Sample(
			std::shared_ptr<const ModelSample<SampleSpace, Model, T>> modelSample, 
			T weight = T(1.0),
			const std::vector<T>& normConstRatioWeights = std::vector<T>(),
			const std::vector<T>& invNormConstRatioWeights = std::vector<T>()) :
		modelSample_(modelSample),
		weight_(weight),
		normConstRatioWeights_(normConstRatioWeights),
		invNormConstRatioWeights_(invNormConstRatioWeights)
	{}

	const Instance& getInstance() const {
		return modelSample_->getInstance();
	}

	std::shared_ptr<const ModelSample<SampleSpace, Model, T>> getModelSample() const {
		return modelSample_;
	}

	void setWeight(T weight) {
		weight_ = weight;
	}
	
	T getWeight() const {
		return weight_;
	}

	void setNumSubSamples(long long int num) {
		nSubSamples_ = num;
	}

	long long int getNumSubSamples() const {
		return nSubSamples_;
	}

	const std::vector<T>& getNormConstWeights() const {
		return normConstRatioWeights_;
	}

	const std::vector<T>& getInvNormConstWeights() const {
		return invNormConstRatioWeights_;
	}

	const RestModel* getRestrictedModel() const {
		return &(modelSample_->getRestrictedModel());
	}
	
	T getUnnormProb() const {
		return modelSample_->getUnnormProb();
	}
	//virtual ~Sample() {}
};



/*
template <typename SampleSpace, typename T>
class Sample {
private:
	typename SampleSpace::Instance instance_;
	T unProb_;
	RestrictedModel<T, SampleSpace>* restrictedModel_;

public:
	Sample(
		const typename SampleSpace::Instance& instance,
		T unProb)
	:
		instance_(instance),
		unProb_(unProb),
		restrictedModel_(nullptr)
	{}

	void set(const typename SampleSpace::Instance& instance, T unProb) {
		instance_ = instance;
		unProb_ = unProb;
	}

	void setRestrictedModel(RestrictedModel<T, SampleSpace>* condModel) {
		restrictedModel_ = condModel;
	}

	const typename SampleSpace::Instance& getInstance() const {
		return instance_;
	}
	
	T getUnnormProb() const {
		return unProb_;
	}
	
	//virtual T getWeight() const = 0;
	
	void writeTo(std::ostream& os) const {
		os << instance_;
		os << " " << log(unProb_);
	}
	
	void readFrom(std::istream& is) {
		is >> instance_;
		double tmp;
		is >> tmp;
		unProb_.setLog(tmp);
	}
	
	virtual ~Sample() {}
};

template <class SampleSpace, class T>
class ISSample : public Sample<SampleSpace, T> {
private:
	T weight_;
	typedef Sample<SampleSpace, T> Parent;
public:
	ISSample(
		const typename SampleSpace::Instance& instance,
		T unProb,
		T weight)
	:
		Parent(instance, unProb),
		weight_(weight)
	{}
	
	void set(const typename SampleSpace::Instance& instance, T unProb, T weight) {
		Parent::set(instance, unProb);
		weight_ = weight;
	}
	
	T getWeight() const {
		return weight_;
	}
	
	void writeTo(std::ostream& os) const {
		Parent::writeTo(os);
		os << " " << log(weight_);
	}
	
	void readFrom(std::istream& is) {
		Parent::readFrom(is);
		double tmp;
		is >> tmp;
		weight_.setLog(tmp);
	}
};
*/

/*template <class SampleSpace, class T>
class NonWeightedSample : public Sample<SampleSpace, T> {
private:
	typedef Sample<SampleSpace, T> Parent;
public:
	ISSample(const typename SampleSpace::Instance& instance, T unProb)
		: Parent(instance, unProb)
	{
	}
	
	void set(const typename SampleSpace::Instance& instance, T unProb) {
		Parent::set(instance, unProb);
	}
	
	T getWeight() const {
		return T(1.0);
	}
	
	void writeTo(std::ostream& os) const {
		Parent::writeTo(os);
		os << " " << log(weight_);
	}
	
	void readFrom(std::istream& is) {
		Parent::readFrom(is);
		double tmp;
		is >> tmp;
		weight_.setLog(tmp);
	}
};*/


template <typename SampleSpace, typename Model, typename T>
class SampleSource {
protected:
	const SampleSpace& sampleSpace_;
	const Model* model_;
	const int nLevels_;
	const int nHarmonicLevels_;
public:
	SampleSource(const SampleSpace& sampleSpace, const Model* model, int nLevels, int nHarmonicLevels)
		: sampleSpace_(sampleSpace), model_(model), nLevels_(nLevels), nHarmonicLevels_(nHarmonicLevels) {}
	virtual ~SampleSource() {}
	virtual std::unique_ptr<Sample<SampleSpace, Model, T>> getSample() = 0;
	virtual bool eof() const = 0;
	virtual const Model* getModel() {
		return model_;
	};
	const SampleSpace& getSampleSpace() {
		return sampleSpace_;
	};
	int nLevels() const {
		return nLevels_;
	}
	int nHarmonicLevels() const {
		return nHarmonicLevels_;
	}

	virtual void setInterruptCallback(std::function<bool()> interruptCallback) {
	}
};

/*mplate <typename T, typename SampleSpace, typename SubSampleSpace>
class NestedSampleSource : SampleSource<T, SampleSpace> {
protected:
	const SampleSpace& sampleSpace_;
public:
	SampleSource(const SampleSpace& sampleSpace)
		: sampleSpace_(sampleSpace) {}
	virtual SampleSource<T, SubSampleSpace>* getSubSampleSource() = 0;
};*/

template <typename SampleSpace, typename Model, typename T>
class SampleSink {
public:
	virtual ~SampleSink() {}
	virtual void init() {};
	virtual void addSample(const Sample<SampleSpace, Model, T>*) = 0;
};







/*template <typename SampleSpace, typename Model, typename T>
std::ostream& operator<<(std::ostream& os, const Sample<SampleSpace, Model, T>& sample) {
	sample.writeTo(os);
	return os;
}

template <typename SampleSpace, typename Model, typename T>
std::istream& operator>>(std::istream& is, Sample<SampleSpace, Model, T>& sample) {
	sample.readFrom(is);
	return is;
}*/


/*template <class SampleSpace, class T>
class Sample {
public:
	virtual const typename SampleSpace::Instance& getInstance() = 0;
	virtual T getUnnormProb() = 0;
	virtual void writeTo(std::ostream& os) const = 0;
	virtual void readFrom(std::istream& is) = 0;
};*/



/*template <class SampleSpace, class T>
class ISSample : public Sample<SampleSpace, T> {
private:
	T drawProb_;
	T weight_;
	bool weightComputed_;
	T unProb_;
	bool unProbComputed_;
	typename SampleSpace::Instance instance_;
	const ParentsetMap<T>& scores_;
public:
	ISSample(const ParentsetMap<T>& scores, SampleSpace sampleSpace)
		: scores_(scores), instance_(sampleSpace)
	{
		unProbComputed_ = false;
		weightComputed_ = false;
	}
	
	void set(const typename SampleSpace::Instance& instance, T drawProb) {
		instance_ = instance;
		drawProb_ = drawProb;
		unProbComputed_ = false;
		weightComputed_ = false;
	}
	
	const typename SampleSpace::Instance& getInstance() {
		return instance_;
	}
	
	T getUnnormProb() {
		if (!unProbComputed_) {
			unProb_ = calcUnnormProb(instance_.family, scores_, instance_, NullArc);
			unProbComputed_ = true;
		}
		return unProb_;
	}
	
	T getWeight() {
		if (!weightComputed_) {
			weight_ = drawProb_ / getUnnormProb();
			weightComputed_ = true;
		}
		return weight_;
	}
	
	void writeTo(std::ostream& os) const {
		os << log(drawProb_) << " ";
		if (unProbComputed_)
			os << true << " " << log(unProb_);
		else
			os << false << " " << 0.0 / 0.0;
		os << " " << instance_;
	}
	
	void readFrom(std::istream& is) {
		double tmp;
		is >> tmp;
		drawProb_.setLog(tmp);
		is >> unProbComputed_ >> tmp;
		unProb_.setLog(tmp);
		is >> instance_;
		weightComputed_ = false;
	}
};/**/

#endif

