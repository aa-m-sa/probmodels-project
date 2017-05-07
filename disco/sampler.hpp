/*
 *  BEANDisco: samplers
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


#include <unordered_set>

#include "sample.hpp"

#ifndef SAMPLER_HPP
#define SAMPLER_HPP

template <typename OuterSampleSpace, typename OuterModel>
class AbstractSubSampler {
public:
	virtual long long int run(const Sample<OuterSampleSpace, OuterModel, Real>* sample) = 0;
	virtual ~AbstractSubSampler() {}
};

template <typename SampleSpace, typename Model>
class Sampler { //: public AbstractSampler {
private:
	//std::shared_ptr<SampleSpace> sampleSpace_;
	//const SampleSpace& sampleSpace_
	//std::shared_ptr<SampleSource<SampleSpace, Model, Real>> source_;
	//SampleSource<SampleSpace, Model, Real>* source_;
	//const Model* model_;
	const int nLevels_;
	const int nHarmonicLevels_;
	std::shared_ptr<AbstractSubSampler<SampleSpace, Model>> subSampler_;
	//ISNormConstRatioEstimator<Real, GeneralModel, SampleSpace>* normConstRatioEstimator_;
	//std::vector<ISEstimator<SampleSpace, Model, Real>*> estimators_;
	std::vector<std::shared_ptr<SampleSink<SampleSpace, Model, Real>>> sinks_;
	//std::vector<SampleSink<SampleSpace, Model, Real>*> sinks_;
	std::function<bool(int sampleNum)> sampleCallback_;
public:
	//Sampler(const Model* model, SampleSource<Real, SampleSpace>* source)
	Sampler(
			//std::shared_ptr<SampleSpace> sampleSpace,
			//const SampleSpace& sampleSpace,
			//std::shared_ptr<SampleSource<SampleSpace, Model, Real>> source = nullptr)
			//SampleSource<SampleSpace, Model, Real>* source = nullptr)
			int nLevels,
			int nHarmonicLevels
			)
	:
		//sampleSpace_(sampleSpace),
		//model_ = model;
		nLevels_(nLevels),
		nHarmonicLevels_(nHarmonicLevels),
		//source_(source),
		subSampler_(nullptr),
		sampleCallback_()
	{
	}
	
	//void setNormConstRatioEstimator(NormConstRatioISEstimator<SampleSpace, Model>* estimator) {
	//	normConstRatioEstimator_ = estimator;
	//}

	/*void setSource(std::shared_ptr<SampleSource<SampleSpace, Model, Real>> source) {
		source_ = source;
	}*/

	int nLevels() const {
		return nLevels_;
	}
	int nHarmonicLevels() const {
		return nHarmonicLevels_;
	}

	void setSubSampler(std::shared_ptr<AbstractSubSampler<SampleSpace, Model>> subSampler) {
		subSampler_ = subSampler;
	}
	
	/*void addEstimator(ISEstimator<SampleSpace, Model, Real>* estimator) {
		estimators_.push_back(estimator);
	}*/
	
	void addSink(std::shared_ptr<SampleSink<SampleSpace, Model, Real>> sink) {
		sinks_.push_back(sink);
	}

	void setCallback(std::function<bool(int sampleNum)> sampleCallback) {
		sampleCallback_ = sampleCallback;
	}
	
	long long int run(SampleSource<SampleSpace, Model, Real>* source, long long int nSamples = 0,
			double samplingTime = 0) {
		//using SampleType = Sample<SampleSpace, Model, Real>;
		using SampleType = Sample<SampleSpace, Model, Real>;

		if (nSamples == 0)
			nSamples = std::numeric_limits<decltype(nSamples)>::max();
		if (samplingTime == 0)
			samplingTime = std::numeric_limits<decltype(samplingTime)>::infinity();

		for (auto sink : sinks_)
			sink->init();

		Timer timer;
		timer.start();

		for (long long int i = 0; i < nSamples; ++i) {

			if (source->eof())
				return i;

			std::unique_ptr<SampleType> sample(source->getSample());

			if (subSampler_) {
				auto nSubSampled = subSampler_->run(sample.get());
				if (nSubSampled < 0)
					return i;
				sample->setNumSubSamples(nSubSampled);
			}

			//for (auto estimator : estimators_)
			//	estimator->addSample(sample.get());
			
			for (auto sink : sinks_)
				sink->addSample(sample.get());

			if (sampleCallback_) {
				bool cont = sampleCallback_(i);
				if (!cont)
					return i + 1;
			}

			if (timer.elapsed() >= samplingTime)
				return i + 1;
		}
		return nSamples;
	}
};

template <typename OuterSampleSpace, typename OuterModel, typename SubSampleSpace, typename SubModel>
class SubSampler : public AbstractSubSampler<OuterSampleSpace, OuterModel> {
private:
	std::shared_ptr<SubSampleSource<SubSampleSpace, SubModel, Real>> subSourceSource_;
	std::shared_ptr<Sampler<SubSampleSpace, SubModel>> sampler_;
	long long int nSamples_;
	double samplingTime_;
	double dryRunFraction_;
public:
	SubSampler(
			std::shared_ptr<SubSampleSource<SubSampleSpace, SubModel, Real>> subSourceSource,
			std::shared_ptr<Sampler<SubSampleSpace, SubModel>> sampler,
			long long int nSamples,
			double samplingTime,
			double dryRunFraction = 0)
	:
		subSourceSource_(subSourceSource),
		sampler_(sampler),
		nSamples_(nSamples),
		samplingTime_(samplingTime),
		dryRunFraction_(dryRunFraction)
	{}

	long long int run(const Sample<OuterSampleSpace, OuterModel, Real>* sample) {
		auto restModel = sample->getRestrictedModel();
		std::shared_ptr<SampleSource<SubSampleSpace, SubModel, Real>> subSource =
			subSourceSource_->newSubSource(restModel);
		//sampler_->setSource(subSource);
		//sampler_->run(nSamples_);
		if (samplingTime_) {
			if (dryRunFraction_) {
				Timer timer;
				timer.start();
				auto nSamples = sampler_->run(subSource.get(), 0, samplingTime_ * dryRunFraction_);
				double dryRunTime = timer.elapsed();
				if (nSamples < 0)
					return nSamples;
				//nSamples *= (1.0 - dryRunFraction_) / dryRunFraction_;
				nSamples *= (samplingTime_ / dryRunTime - 1.0);
				if (nSamples < 1)
					nSamples = 1;
				return sampler_->run(subSource.get(), nSamples);
			}
			else {
				return sampler_->run(subSource.get(), 0, samplingTime_);
			}
		}
		else {
			return sampler_->run(subSource.get(), nSamples_);
		}
	}
};


using EWDag = ParentSetMap<Real>::ParentSets;

template <typename SampleSpace, typename Model, typename T>
class EWUniqueSampleSource : public SampleSource<SampleSpace, Model, T> {
private:
	using SampleSource<SampleSpace, Model, T>::sampleSpace_;
	using SampleSource<SampleSpace, Model, T>::model_;
	std::shared_ptr<SampleSource<SampleSpace, Model, T>> baseSource_;

	//std::unordered_set<typename SampleSpace::Instance>& allUniqueSamples_;
	std::unordered_set<EWDag, EWDag::hash>& allUniqueSamples_;
	size_t maxTotalSubSamples_;
	mutable bool outOfSpace_;

	//mutable std::unordered_set<typename SampleSpace::Instance> uniqueSamples_;
	mutable std::unordered_set<EWDag, EWDag::hash> uniqueSamples_;
	const Real requiredTotalProb_;

	// TODO: sämplättyjen dagien kokonaismäärän laskenta?

	mutable Real totalProb_;
	mutable std::unique_ptr<Sample<SampleSpace, Model, T>> nextSample_;
public:
	EWUniqueSampleSource(
			std::shared_ptr<SampleSource<SampleSpace, Model, T>> baseSource,
			//std::unordered_set<typename SampleSpace::Instance>& allUniqueSamples,
			std::unordered_set<EWDag, EWDag::hash>& allUniqueSamples,
			size_t maxTotalSubSamples,
			Real requiredTotalProb) :
		SampleSource<SampleSpace, Model, T>(baseSource->getSampleSpace(), baseSource->getModel(), 0, 0),
		baseSource_(baseSource),
		allUniqueSamples_(allUniqueSamples),
		maxTotalSubSamples_(maxTotalSubSamples),
		requiredTotalProb_(requiredTotalProb)
	{
		outOfSpace_ = false;
		totalProb_ = 0;
		nextSample_ = nullptr;

		// hack that avoids divide by zero in estimator if no samples
		// TODO: should be repaired in estimator
		nextSample_ = std::unique_ptr<Sample<SampleSpace, Model, T>>(
				new Sample<SampleSpace, Model, T>(std::make_shared<const ModelSample<SampleSpace, Model, T>>(
					baseSource->getModel(), typename SampleSpace::Instance(baseSource->getSampleSpace())), 0.0));
	}
	
	std::unique_ptr<Sample<SampleSpace, Model, T>> getSample() {
		return std::move(nextSample_);
	}

	bool eof() const {
		if (nextSample_)
			return false;
		while (true) {
			if (totalProb_ >= requiredTotalProb_)
				return true;

			if (allUniqueSamples_.size() + uniqueSamples_.size() >= maxTotalSubSamples_) {
				outOfSpace_ = true;
				return true;
			}

			EWDag sampleParentSets(
					*(model_->getParentModel()->getParentSetFactors()));
			do {
				if (baseSource_->eof()) {
					// TODO: is this okay?
					assert(0);
				}
				nextSample_ = baseSource_->getSample();
				sampleParentSets = nextSample_->getInstance();
			} while (uniqueSamples_.count(sampleParentSets));
			//} while (uniqueSamples_.count(nextSample_->getInstance()));

			uniqueSamples_.insert(sampleParentSets);
			//uniqueSamples_.insert(nextSample_->getInstance());
			Real prob = nextSample_->getRestrictedModel()->getModularUnnormProb();
			nextSample_->setWeight(prob);
			totalProb_ += prob;

			//if (!allUniqueSamples_.count(nextSample_->getInstance())) {
			//	allUniqueSamples_.insert(nextSample_->getInstance());
			if (!allUniqueSamples_.count(sampleParentSets)) {
				allUniqueSamples_.insert(sampleParentSets);
				return false;
			}
			nextSample_ = nullptr;
		}
	}

	bool outOfSpace() {
		return outOfSpace_;
	}
};


template <typename OuterSampleSpace, typename OuterModel, typename SubSampleSpace, typename SubModel>
class EWSubSampler : public AbstractSubSampler<OuterSampleSpace, OuterModel> {
private:
	std::shared_ptr<SubSampleSource<SubSampleSpace, SubModel, Real>> subSourceSource_;
	std::shared_ptr<Sampler<SubSampleSpace, SubModel>> sampler_;
	double totalProbFractionEps_;
	size_t maxTotalSubSamples_;
	//std::unordered_set<typename SubSampleSpace::Instance> uniqueSubSamples_;
	std::unordered_set<EWDag, EWDag::hash> uniqueSubSamples_;
public:
	EWSubSampler(
			std::shared_ptr<SubSampleSource<SubSampleSpace, SubModel, Real>> subSourceSource,
			std::shared_ptr<Sampler<SubSampleSpace, SubModel>> sampler,
			double totalProbFractionEps,
			size_t maxTotalSubSamples)
	:
		subSourceSource_(subSourceSource),
		sampler_(sampler),
		totalProbFractionEps_(totalProbFractionEps),
		maxTotalSubSamples_(maxTotalSubSamples)
	{}

	long long int run(const Sample<OuterSampleSpace, OuterModel, Real>* sample) {
		if (sample->getWeight() != Real(1)) {
			// TODO: korjaa että AISissa painot takaisin 1:ksi
			printf("Error: EW-sampling with weighted (partial?)order samples is not implemented.\n");
			exit(1);
		}
		auto restModel = sample->getRestrictedModel();
		std::shared_ptr<SampleSource<SubSampleSpace, SubModel, Real>> subSource =
			subSourceSource_->newSubSource(restModel);

		Real requiredTotalProb = restModel->getUnnormProb() * (1 - totalProbFractionEps_);
		EWUniqueSampleSource<SubSampleSpace, SubModel, Real>
				uniqueSubSource(subSource, uniqueSubSamples_, maxTotalSubSamples_, requiredTotalProb);

		auto nSampled = sampler_->run(&uniqueSubSource);

		if (uniqueSubSource.outOfSpace())
			return -1;

		return nSampled;
	}
};


/*template <typename SampleSpace, typename SourceModel, typename TargetModel, typename T>
class SampleBiaser : public  SampleSource<SampleSpace, TargetModel, T> {
private:
	using SampleSource<SampleSpace, TargetModel, T>::sampleSpace_;
	std::shared_ptr<SampleSource<SourceModel, SourceModel, T>> source_;
	//SampleEstimator<SampleSpace, SourceModel, T>* weightEstimator_;
	//Estimator<RestrictedModel<SourceModel, SampleSpace>, T>* weightEstimator_;
	Estimator<typename SampleSpace::Instance, T>* weightEstimator_;
	const TargetModel* targetModel_;
public:
	SampleBiaser(
			std::shared_ptr<SampleSource<SampleSpace, SourceModel, T>> source,
			//SampleEstimator<SampleSpace, Model, T>* weightEstimator,
			//Estimator<RestrictedModel<SourceModel, SampleSpace>, T>* weightEstimator,
			Estimator<typename SampleSpace::Instance, T>* weightEstimator,
			const TargetModel* targetModel) :
		source_(source),
		weightEstimator_(weightEstimator),
		targetModel_(targetModel)
	{
	}

	std::unique_ptr<Sample<SampleSpace, TargetModel, T>> getSample() {
		auto sourceSample = source_->getSample();
		//weightEstimator_->addSample(sourceSample);
		//T weight = weightEstimator_->getEstimate();
		//T weight = weightEstimator_->getEstimate(sourceSample.getRestrictedModel());
		T weight = weightEstimator_->getEstimate(sourceSample.getInstance());
		weight *= sourceSample.getWeight();
		auto modelSample = std::make_shared<ModelSample<SampleSpace, TargetModel, T>>(
				targetModel_, sourceSample.getInstance());
		auto sample = new Sample<SampleSpace, TargetModel, T>(modelSample, weight);
		return sample;
	}
};*/

#endif

