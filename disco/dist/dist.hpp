/*
 *  BEANDisco: classes representing distributions (sample generators)
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


#include "../sample.hpp"
#include "../estimator.hpp"

#ifndef DIST_HPP
#define DIST_HPP


/**
 * Abstract distribution
 */
/*template <class SampleSpace, class T>
class Dist {
protected:
	const SampleSpace& sampleSpace_;
	int nSamples_;
public:
	Dist(const SampleSpace& sampleSpace) 
		: sampleSpace_(sampleSpace), nSamples_(0) {
	}
	const SampleSpace& getSampleSpace() {
		return sampleSpace_;
	}
	virtual void rand() = 0;
	//virtual void getSample(Sample<SampleSpace, T>& sample) = 0;
	virtual Sample<SampleSpace, Model, T>* getSample() = 0;
	int getNSamples() {
		return nSamples_;
	}
	virtual ~Dist() {}
};*/

/**
 * Abstract distribution with probability model parameter
 */
/*template <class SampleSpace, class Model, class T>
class ModelDist : public Dist<SampleSpace, T> {
protected:
	const Model* model_;
	using Dist<SampleSpace, T>::sampleSpace_;
	using Dist<SampleSpace, T>::nSamples_;
public:
	ModelDist(const SampleSpace& sampleSpace, const Model* model) :
			Dist<SampleSpace, T>(sampleSpace),
			model_(model)
	{}
	//virtual void setModel(Model* model) {
	//	model_ = model;
	//}
};*/


/**
 * Abstract distribution with probability model parameter
 */
template <class SampleSpace, class Model, class T>
class ModelDist {
public:
	//using SampleModel = typename Model::template RestrictedModel<SampleSpace>;
	using SampleModel = RestrictedModel<Model, SampleSpace>;
protected:
	const SampleSpace& sampleSpace_;
	const Model* model_;
	long long int nSamples_;

	std::function<bool()> interruptCallback_;
public:
	ModelDist(const SampleSpace& sampleSpace, const Model* model) :
			sampleSpace_(sampleSpace),
			model_(model),
			nSamples_(0),
			interruptCallback_(nullptr)
	{}

	virtual ~ModelDist() {}

	const SampleSpace& getSampleSpace() {
		return sampleSpace_;
	}

	virtual std::unique_ptr<Sample<SampleSpace, Model, T>> rand() = 0;

	//virtual void setModel(Model* model) {
	//	model_ = model;
	//}
	const Model* getModel() {
		return model_;
	}

	long long int getNSamples() {
		return nSamples_;
	}

	virtual int nLevels() const {
		return 0;
	}
	virtual int nHarmonicLevels() const {
		return 0;
	}

	virtual std::unique_ptr<ISEstimator<SampleSpace, Model, T, T>>
	newDirectUnweightedNormConstEstimator(double meanFraction) const {
		return nullptr;
	}

	virtual std::unique_ptr<ISEstimator<SampleSpace, Model, T, T>>
	newHarmonicUnweightedNormConstEstimator(double meanFraction) const {
		return nullptr;
	}

	void setInterruptCallback(std::function<bool()> interruptCallback) {
		interruptCallback_ = interruptCallback;
	}
};


/**
 * Abstract distribution with multiple temperature levels
 */
template <typename SampleSpace, typename Model, typename T>
class MultiLevelModelDist : public ModelDist<SampleSpace, Model, T> {
protected:
	//using Dist<SampleSpace, T>::sampleSpace_;
	using ModelDist<SampleSpace, Model, T>::nSamples_;
	using ModelDist<SampleSpace, Model, T>::sampleSpace_;
	using ModelDist<SampleSpace, Model, T>::model_;
	using ModelDist<SampleSpace, Model, T>::interruptCallback_;

	int nLevels_;
	int nHarmonicLevels_;
public:
	
	MultiLevelModelDist(const SampleSpace& sampleSpace, const Model* model, int nLevels, int nHarmonicLevels) :
		ModelDist<SampleSpace, Model, T>(sampleSpace, model),
		nLevels_(nLevels),
		nHarmonicLevels_(nHarmonicLevels)
	{
	}
	
	int nLevels() const {
		return nLevels_;
	}
	int nHarmonicLevels() const {
		return nHarmonicLevels_;
	}
	//virtual void getLevelRatios(std::vector<T>& levelRatios) = 0;
	//virtual double acceptRatio(int level) = 0;
	//virtual double swapAcceptRatio(int level) = 0;
};



#endif