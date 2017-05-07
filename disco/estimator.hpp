/*
 *  BEANDisco: estimators
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
#include <vector>

#include "sample.hpp"
#include "model.hpp"
#include "feature.hpp"
#include "bucketorder.hpp"
#include "parentsetmap.hpp"
//#include "integrator.hpp"
//#include "bodist.hpp"
#include "io.hpp"
#include "linext.hpp"
#include "partialsums.hpp"

#ifndef ESTIMATOR_HPP
#define ESTIMATOR_HPP



template <typename T>
class Mean {
public:
	virtual ~Mean() {}
	virtual void reset() = 0;
	virtual void add(const T& x) = 0;
	virtual T get() const = 0;
};

template <typename T>
class FullMean : public Mean<T> {
private:
	T zero_;
	T sum_;
	size_t n_;
public:
	FullMean(const T& zero) :
		zero_(zero),
		sum_(zero_),
		n_(0)
	{
	}
	void reset() {
		sum_ = zero_;
		n_ = 0;
	}
	void add(const T& x) {
		sum_ += x;
		++n_;
	}
	T get() const {
		return sum_ / n_;
	}
};


/*template <typename T>
class TailMean : public Mean<T>{
private:
	double tailFraction_;
	T zero_;
	std::vector<T> cumSum_;
public:
	TailMean(double tailFraction, const T& zero) :
		tailFraction_(tailFraction),
		zero_(zero)
	{
		reset();
	}
	void reset() {
		cumSum_.clear();
		cumSum_.push_back(zero_);
	}
	void add(const T& x) {
		cumSum_.push_back(cumSum_.back() + x);
	}
	T get() const {
		size_t n = size_t(ceil((cumSum_.size() - 1) * tailFraction_));
		return (cumSum_.back() - cumSum_[cumSum_.size() - 1 - n]) / n;
	}
};*/

template <typename T, typename PartialSums>
class TailMean : public Mean<T>{
private:
	double tailFraction_;
	PartialSums partialSums_;
public:
	TailMean(double tailFraction, const T& zero) :
		tailFraction_(tailFraction),
		partialSums_(zero)
	{
		reset();
	}
	void reset() {
		partialSums_.reset();
	}
	void add(const T& x) {
		partialSums_.add(x);
	}
	T get() const {
		size_t end = partialSums_.size();
		size_t begin = end - size_t(ceil(end * tailFraction_));
		return partialSums_.get(begin, end) / (end - begin);
	}
};


/*template <typename T>
class TailMean2 : public Mean<T>{
private:
	double tailFraction_;
	T zero_;
	PartialSums<T> partialSums_;
	PartialSums<T>::SlidingSum slidingSum_;
	size_t n;
public:
	TailMean(double tailFraction, const T& zero) :
		tailFraction_(tailFraction),
		partialSums_(zero),
		slidingSum_(partialSums_)
	{
		//reset();
	}
	void reset() {
		partialSums_.reset();
		slidingSum_.reset();
	}
	void add(const T& x) {
		partialSums_.add(x);
		int n = partialSums_.size();
		slidingSum_.moveEnd();
		if (ceil((n - 1) * tailFraction_) > ceil((n - 2) * tailFraction_))
			slidingSum_.moveBegin();
	}
	T get() const {
		return slidingSum_.get();
	}
};/**/


template <typename T>
std::unique_ptr<Mean<T>> newMean(double tailFraction, const T& zero) {
	assert(0 < tailFraction && tailFraction <= 1.0);
	if (tailFraction == 1.0) {
		return std::unique_ptr<Mean<T>>(new FullMean<T>(zero));
	}
	else {
		//return std::unique_ptr<Mean<T>>(new TailMean<T>(tailFraction, zero));
		//return std::unique_ptr<Mean<T>>(new TailMean<T, PartialSumsNaive<T>>(tailFraction, zero));
		//return std::unique_ptr<Mean<T>>(new TailMean<T, PartialSumsSubtract<T>>(tailFraction, zero));
		return std::unique_ptr<Mean<T>>(new TailMean<T, PartialSumsTree<T>>(tailFraction, zero));
	}
}


/*template <typename T>
class MeanInstance {
public:
	virtual void reset() = 0;
	virtual void add(const T& x) = 0;
	virtual T get() const = 0;
};

template <typename T>
class MeanType {
public:
	virtual MeanInstance<T>* newInstance() const = 0;
};


template <class T>
class FullMean : public MeanType<T> {
public:
	class Instance : public MeanInstance<T> {
	private:
		T sum_;
		size_t n_;
	public:
		Instance() {
			reset();
		}
		void reset() {
			sum_ = T(0);
			n_ = 0;
		}
		void add(const T& x) {
			sum_ += x;
			++n_;
		}
		T get() const {
			return sum_ / n_;
		}
	};
	MeanInstance<T>* newInstance() const {
		return new Instance();
	}
};


template <class T>
class TailMean : public MeanType<T>{
private:
	double tailFraction_;
public:
	TailMean(double tailFraction) {
		tailFraction_ = tailFraction;
	}
	class Instance : public MeanInstance<T> {
	private:
		double tailFraction_;
		std::vector<T> cumSum_;
	public:
		Instance(double tailFraction) {
			tailFraction_ = tailFraction;
			reset();
		}
		void reset() {
			cumSum_.clear();
			cumSum_.push_back(T(0));
		}
		void add(const T& x) {
			cumSum_.push_back(cumSum_.back() + x);
		}
		T get() const {
			int n = int(ceil((cumSum_.size() - 1) * tailFraction_));
			return (cumSum_.back() - cumSum_[cumSum_.size() - 1 - n]) / n;
		}
	};
	MeanInstance<T>* newInstance() const {
		return new Instance();
	}
};*/







template <typename FT>
class Estimator {
public:
	virtual FT getEstimate() const = 0;
};

template <typename SampleSpace, typename Model, typename FT, typename T>
class ISEstimator :
		public Estimator<FT>,
		public SampleSink<SampleSpace, Model, T> {
public:
	//using SampleType = Sample<SampleSpace, Model, T>;
	//typedef Sample<SampleSpace, Model, T> SampleType;
};

template <typename SampleSpace, typename Model, typename FT, typename T>
class SampleEstimator :
		public Estimator<FT>,
		public SampleSink<SampleSpace, Model, T> {
public:
};

template <typename SampleSpace, typename Model, typename FT, typename T>
class WeightedEstimator : public SampleEstimator<SampleSpace, Model, FT, T> {
private:
	using SampleType = Sample<SampleSpace, Model, T>;
	std::shared_ptr<const Estimator<FT>> baseEstimator_;
	std::shared_ptr<const Estimator<T>> weightEstimator_;
	FT estimate_;
public:
	WeightedEstimator(
			std::shared_ptr<const Estimator<FT>> baseEstimator,
			std::shared_ptr<const Estimator<T>> weightEstimator,
			const FT& zeroValue) :
		baseEstimator_(baseEstimator),
		weightEstimator_(weightEstimator),
		estimate_(zeroValue)
	{}
	void addSample(const SampleType* sample) {
		estimate_ = baseEstimator_->getEstimate();
		estimate_ *= weightEstimator_->getEstimate();
	}
	FT getEstimate() const {
		return estimate_;
	}
};


template <typename SampleSpace, typename Model, typename T>
class ExactLinExtCountInvEstimator : public SampleEstimator<SampleSpace, Model, T, T> {
private:
	T value_;
public:
	ExactLinExtCountInvEstimator()
	{}
	void addSample(const Sample<SampleSpace, Model, T>* sample) {
		const DagFamily::Instance& dag = sample->getInstance();
		int n = dag.family.n;
		SquareMatrix<bool> adjMat(n);
		dag.getAdjMat(adjMat);
		SquareMatrix<bool> predMat(n);
		adjMatToPredMat(adjMat, predMat);
		if (n <= 20) {
			value_ = T(1.0) / countLinExtsExact<long long unsigned int>(predMat);
		}
		else if (n <= 170) {
			value_ = T(1.0) / countLinExtsExact<double>(predMat);
		}
		else {
			value_ = T(1.0) / countLinExtsExact<T>(predMat);
		}
	}
	T getEstimate() const {
		return value_;
	};
};



/*template <typename T>
class NormConstRatioEstimator : public Estimator<T> {
};*/


template <typename T>
class NormConstEstimator : public Estimator<T> {
private:
	std::shared_ptr<const Estimator<T>> unweightedNormConstEstimator_;
	std::shared_ptr<const Estimator<T>> normConstRatioEstimator_;
public:
	NormConstEstimator(
			std::shared_ptr<const Estimator<T>> unweightedNormConstEstimator,
			std::shared_ptr<const Estimator<T>> normConstRatioEstimator) :
		unweightedNormConstEstimator_(unweightedNormConstEstimator),
		normConstRatioEstimator_(normConstRatioEstimator)
	{}

	T getEstimate() const {
		return unweightedNormConstEstimator_->getEstimate() * normConstRatioEstimator_->getEstimate();
	}
};

template <typename FT, typename T>
class FeatureProbEstimator : public Estimator<FT> {
private:
	std::shared_ptr<const Estimator<FT>> unnormFeatureProbEstimator_;
	std::shared_ptr<const Estimator<T>> normConstRatioEstimator_;
public:
	FeatureProbEstimator(
			std::shared_ptr<const Estimator<FT>> unnormFeatureProbEstimator,
			std::shared_ptr<const Estimator<T>> normConstRatioEstimator) :
		unnormFeatureProbEstimator_(unnormFeatureProbEstimator),
		normConstRatioEstimator_(normConstRatioEstimator)
	{}

	FT getEstimate() const {
		return unnormFeatureProbEstimator_->getEstimate() / normConstRatioEstimator_->getEstimate();
	}
};


/*template <typename T>
class ArcProbsEstimator : public Estimator<ArcMap<T>> {
private:
	std::shared_ptr<const Estimator<FT>> unnormFeatureProbEstimator_;
	std::shared_ptr<const Estimator<T>> normConstRatioEstimator_;
public:
	FeatureProbEstimator(
			std::shared_ptr<const Estimator<FT>> unnormFeatureProbEstimator,
			std::shared_ptr<const Estimator<T>> normConstRatioEstimator) :
		unnormFeatureProbEstimator_(unnormFeatureProbEstimator),
		normConstRatioEstimator_(normConstRatioEstimator)
	{}

	ArcMap<T> getEstimate() const {
		return unnormFeatureProbEstimator_->getEstimate() / normConstRatioEstimator_->getEstimate();
	}
};*/


template <typename SampleSpace, typename Model, typename T>
class MultiLevelNormConstEstimator : public SampleEstimator<SampleSpace, Model, T, T> {
public:
	MultiLevelNormConstEstimator()
	{}
	void addSample(const Sample<SampleSpace, Model, T>* sample) {
	}
	T getEstimate() const {
		return T(1.0);
	};
};



template <typename SampleSpace, typename Model, typename T>
class ExactNormConstRatioEstimator : public SampleEstimator<SampleSpace, Model, T, T> {
public:
	ExactNormConstRatioEstimator()
	{}
	void addSample(const Sample<SampleSpace, Model, T>* sample) {
	}
	T getEstimate() const {
		return T(1.0);
	};
};

template <typename SampleSpace, typename Model, typename Feature, typename T>
class ExactFeatureProbEstimator : public SampleEstimator<SampleSpace, Model, typename Feature::ValueType, T> {
private:
	Feature feature_;
	typename Feature::ValueType prob_;
public:
	ExactFeatureProbEstimator(const Feature& feature) :
		feature_(feature),
		prob_(feature.zeroValue())
	{}
	void addSample(const Sample<SampleSpace, Model, T>* sample) {
		prob_ = sample->getRestrictedModel()->getFeatureProb(feature_);
	}
	typename Feature::ValueType getEstimate() const {
		return prob_;
	}
};


template <typename SampleSpace, typename Model, typename T>
class ExactArcProbsEstimator : public SampleEstimator<SampleSpace, Model, ArcMap<T>, T> {
private:
	ArcMap<T> probs_;
public:
	ExactArcProbsEstimator(const SampleSpace& sampleSpace) :
		probs_(sampleSpace.n)
	{}
	void addSample(const Sample<SampleSpace, Model, T>* sample) {
		sample->getRestrictedModel()->getArcProbs(probs_);
	}
	ArcMap<T> getEstimate() const {
		return probs_;
	}
};


template <typename SampleSpace, typename Model, typename T>
class ISNormConstRatioEstimator : public ISEstimator<SampleSpace, Model, T, T> {
private:
	//std::shared_ptr<const MeanType<T>> estimateType_;
	//std::unique_ptr<MeanInstance<T>> estimate_;
	std::unique_ptr<Mean<T>> estimate_;
	std::shared_ptr<const Estimator<T>> sampleEstimator_;
public:
	ISNormConstRatioEstimator(
			//std::shared_ptr<const MeanType<T>> meanType,
			double meanFraction,
			std::shared_ptr<const Estimator<T>> sampleEstimator) :
		//estimateType_(meanType),
		//estimate_(estimateType_->newInstance()),
		estimate_(newMean<T>(meanFraction, T(0))),
		sampleEstimator_(sampleEstimator)
	{
	}
	
	void init() {
		estimate_->reset();
	}
	
	void addSample(const Sample<SampleSpace, Model, T>* sample) {
		estimate_->add(sample->getWeight() * sampleEstimator_->getEstimate());
	}
	
	T getEstimate() const {
		return estimate_->get();
	}
};


template <typename SampleSpace, typename Model, typename FT, typename T>
class ISUnnormFeatureProbEstimator : public ISEstimator<SampleSpace, Model, FT, T> {
private:
	//std::shared_ptr<const MeanType<FT>> estimateType_;
	//std::unique_ptr<MeanInstance<FT>> estimate_;
	std::unique_ptr<Mean<FT>> estimate_;
	std::shared_ptr<const Estimator<FT>> sampleEstimator_;
public:
	ISUnnormFeatureProbEstimator(
			//std::shared_ptr<const MeanType<FT>> meanType,
			double meanFraction,
			std::shared_ptr<const Estimator<FT>> sampleEstimator,
			const FT& zeroValue) :
		//estimateType_(meanType),
		//estimate_(estimateType_->newInstance()),
		estimate_(newMean<FT>(meanFraction, zeroValue)),
		sampleEstimator_(sampleEstimator)
	{
	}
	
	void init() {
		estimate_->reset();
	}
	
	void addSample(const Sample<SampleSpace, Model, T>* sample) {
		estimate_->add(sampleEstimator_->getEstimate() * sample->getWeight());
	}
	
	FT getEstimate() const {
		return estimate_->get();
	}
};


template <typename SampleSpace, typename Model, typename T>
class ISMultiLevelNormConstEstimator : public ISEstimator<SampleSpace, Model, T, T> {
private:
	std::vector<std::unique_ptr<Mean<T>>> estimates_;
public:
	ISMultiLevelNormConstEstimator(
			int nLevels,
			double meanFraction) :
		estimates_(nLevels)
	{
		for (int i = 0; i < nLevels; ++i)
			estimates_[i] = newMean<T>(meanFraction, T(0));
	}
	
	void init() {
		for (auto& estimate : estimates_)
			estimate->reset();
	}
	
	void addSample(const Sample<SampleSpace, Model, T>* sample) {
		auto& weights = sample->getNormConstWeights();
		assert(estimates_.size() == weights.size());
		for (int i = 0; i < estimates_.size(); ++i)
			estimates_[i]->add(weights[i]);
	}
	
	T getEstimate() const {
		T prodEst = 1.0;
		for (auto& estimate : estimates_)
			prodEst *= estimate->get();
		return prodEst;
	}
};

template <typename SampleSpace, typename Model, typename T>
class ISMultiLevelHarmonicNormConstEstimator : public ISEstimator<SampleSpace, Model, T, T> {
private:
	std::vector<std::unique_ptr<Mean<T>>> estimates_;
public:
	ISMultiLevelHarmonicNormConstEstimator(
			int nLevels,
			double meanFraction) :
		estimates_(nLevels)
	{
		for (int i = 0; i < nLevels; ++i)
			estimates_[i] = newMean<T>(meanFraction, T(0));
	}
	
	void init() {
		for (auto& estimate : estimates_)
			estimate->reset();
	}
	
	void addSample(const Sample<SampleSpace, Model, T>* sample) {
		auto& weights = sample->getInvNormConstWeights();
		assert(estimates_.size() == weights.size());
		for (int i = 0; i < estimates_.size(); ++i)
			estimates_[i]->add(weights[i]);
	}
	
	T getEstimate() const {
		T prodEst = 1.0;
		for (auto& estimate : estimates_)
			prodEst /= estimate->get();
		return prodEst;
	}
};


/*
template <typename T, typename Model, typename Feature, class SampleSpace>
class ISUnnormalizedFeatureProbEstimator : public FeatureProbEstimator<T, Model, Feature>,
		public ISEstimator<T, SampleSpace> {
public:
	using RestModel = RestrictedModel<T, SampleSpace>;
	using SampleType = Sample<SampleSpace, T>;
	using SampleEstimator = FeatureProbEstimator<T, RestModel, Feature>;
private:
	MeanType<T> estimateType_;
	std::unique_ptr<MeanInstance<T>> estimate_;
	const SampleEstimator& sampleEstimator_;
public:
	ISFeatureProbEstimator(const MeanType<T>& meanType, const SampleEstimator& sampleEstimator) :
		estimateType_(meanType),
		estimate_(estimateType_->newInstance()),
		sampleEstimator_(sampleEstimator)
	{
	}
	
	void reset() {
		estimate_->reset();
	}
	
	void addSample(const SampleType* sample) {
		estimate_->add(sample->getWeight() * sampleEstimator_.getEstimate(sample->getRestrictedModel()));
	}
	
	T getUnnormEstimate(const Model& model) {
		return estimate_->get();
	}
};
*/























/*

template <typename T, typename Model>
class NormConstRatioEstimator {
public:
	virtual T getEstimate(const Model& model) = 0;
};

template <typename T, typename Model>
class NormConstEstimator {
protected:
	NormConstRatioEstimator<T, Model> ncre;
public:
	//virtual T getUnnormEstimate(const Model& model) = 0;
	T getEstimate(const Model& model) {
		return ncre.getEstimate(model);
	}
};

template <typename T, typename Model, typename Feature>
class FeatureProbEstimator {
protected:
	NormConstRatioEstimator<T, Model> ncre;
public:
	virtual T getUnnormEstimate(const Model& model) = 0;
	T getEstimate(const Model& model) {
		return getUnnormEstimate(model) / ncre.getEstimate(model);
	}
};

template <typename T, typename Model>
class ExactNormConstRatioEstimator : public NormConstRatioEstimator<T, Model> {
public:
	T getEstimate(const Model* model) {
		return T(1.0);
	}
};

template <typename T, typename Model, typename Feature>
class ExactFeatureProbEstimator : public FeatureProbEstimator<T, Model, Feature> {
private:
	const Feature& feature_;
public:
	//T getUnnormEstimate(const Model& model) {
	//	return model.template getUnnormFeatureProb<Feature>(feature_);
	//}
	T getEstimate(const Model* model) {
		return model->template getFeatureProb<Feature>(feature_);
	}
};


template <typename T, typename Model, class SampleSpace>
class ISNormConstRatioEstimator : public NormConstRatioEstimator<T, Model>,
		public ISEstimator<T, SampleSpace> {
public:
	using RestModel = RestrictedModel<T, SampleSpace>;
	using SampleType = Sample<SampleSpace, T>;
	using SampleEstimator = NormConstRatioEstimator<T, RestModel>;
private:
	MeanType<T> estimateType_;
	std::unique_ptr<MeanInstance<T>> estimate_;
	const SampleEstimator& sampleEstimator_;
public:
	ISNormConstRatioEstimator(const MeanType<T>& meanType, const SampleEstimator& sampleEstimator) :
		estimateType_(meanType),
		estimate_(estimateType_->newInstance()),
		sampleEstimator_(sampleEstimator)
	{
	}
	
	void reset() {
		estimate_->reset();
	}
	
	void addSample(const SampleType* sample) {
		estimate_->add(sample->getWeight() * sampleEstimator_.getEstimate(sample->getRestrictedModel()));
	}
	
	T getEstimate(const Model& model) {
		return estimate_->get();
	}
};


template <typename T, typename Model, typename Feature, class SampleSpace>
class ISFeatureProbEstimator : public FeatureProbEstimator<T, Model, Feature>,
		public ISEstimator<T, SampleSpace> {
public:
	using RestModel = RestrictedModel<T, SampleSpace>;
	using SampleType = Sample<SampleSpace, T>;
	using SampleEstimator = FeatureProbEstimator<T, RestModel, Feature>;
private:
	MeanType<T> estimateType_;
	std::unique_ptr<MeanInstance<T>> estimate_;
	const SampleEstimator& sampleEstimator_;
public:
	ISFeatureProbEstimator(const MeanType<T>& meanType, const SampleEstimator& sampleEstimator) :
		estimateType_(meanType),
		estimate_(estimateType_->newInstance()),
		sampleEstimator_(sampleEstimator)
	{
	}
	
	void reset() {
		estimate_->reset();
	}
	
	void addSample(const SampleType* sample) {
		estimate_->add(sample->getWeight() * sampleEstimator_.getEstimate(sample->getRestrictedModel()));
	}
	
	T getUnnormEstimate(const Model& model) {
		return estimate_->get();
	}
};
*/

/*
template <class SampleSpace, class T>
class ISEstimate {
protected:
	OutStream* resStream_;
	
public:
	ISEstimate(OutStream* resStream) : resStream_(resStream) {}
	virtual ~ISEstimate() {}
	virtual void reset() = 0;
	virtual void addSample(const Sample<SampleSpace, T>* sample) = 0;
	virtual void updateCurrentEstimate() = 0;
};



template <class SampleSpace, class T, class Mean>
class LogNormConstEstimate : public ISEstimate<SampleSpace, T> {
private:
	Mean sumType_;
	typename Mean::Instance totalWeight_;
	size_t nSamples_;
public:
	using ISEstimate<SampleSpace, T>::resStream_;
	
	LogNormConstEstimate(OutStream* resStream, Mean sumType) :
		ISEstimate<SampleSpace, T>(resStream), sumType_(sumType)
	{
		reset();
	}
	
	void reset() {
		totalWeight_.reset();
	}
	
	void addSample(const Sample<SampleSpace, T>* sample) {
		totalWeight_.add(sumType_, getWeight(sample));
	}
	
	void updateCurrentEstimate() {
		resStream_->newField();
		resStream_->stream() << log(totalWeight_.get(sumType_)) << std::endl;
		resStream_->endField();
	}
};

template <class SampleSpace, class T, class Mean>
class MultiLevelLogNormConstEstimate {
private:
	OutStream* resStream_;
	Mean sumType_;
	std::vector<typename Mean::Instance> totalWeights_;
	size_t nSamples_;
public:
	MultiLevelLogNormConstEstimate(OutStream* resStream, size_t nLevels, Mean sumType) :
		resStream_(resStream), totalWeights_(nLevels), sumType_(sumType)
	{
		reset();
	}
	
	void reset() {
		for (int i = 0; i < totalWeights_.size(); ++i)
			totalWeights_.reset();
	}
	
	void addRatioSamples(const std::vector<T> ratioSamples) {
		assert(totalWeights_.size() == ratioSamples.size());
		for (int i = 0; i < totalWeights_.size(); ++i)
			totalWeights_[i].add(sumType_, ratioSamples[i]);
		++nSamples_;
	}
	
	void updateCurrentEstimate() {
		T estimate = T(1.0);
		for (int i = 0; i < totalWeights_.size(); ++i)
			estimate *= totalWeights_[i].get(sumType_);
		resStream_->newField();
		resStream_->stream() << log(estimate) << std::endl;
		resStream_->endField();
	}
};


template <class SampleSpace, class T, class Feature, class Mean>
class FeatureProbEstimate : public ISEstimate<SampleSpace, T> {
private:
	Feature feature_;
	Mean sumType_;
	typename Mean::Instance featureWeight_;
	typename Mean::Instance totalWeight_;
	const SampleSpace& sampleSpace_;
	const Model<T>& model_;
public:
	using ISEstimate<SampleSpace, T>::resStream_;
	
	FeatureProbEstimate(OutStream* resStream,
	                 const SampleSpace& sampleSpace,
	                 const Model<T>& model,
	                 const Feature& feature,
	                 Mean sumType) :
		ISEstimate<SampleSpace, T>(resStream),
		sampleSpace_(sampleSpace),
		model_(model),
		feature_(feature),
		sumType_(sumType)
	{
		reset();
	}
	
	void reset() {
		featureWeight_.reset();
		totalWeight_.reset();
	}
	
	void addSample(const Sample<SampleSpace, T>* sample) {
		T weight = getWeight(sample);
		T mp = sample->getUnnormProb();
		T fp = model_.template getUnnormFeatureProb<SampleSpace, Feature>(
				sample->getInstance(), feature_);
		T featProb = fp / mp;
		featureWeight_.add(sumType_, weight * featProb);
		totalWeight_.add(sumType_, weight);
	}
	
	void updateCurrentEstimate() {
		T estimate = featureWeight_.get(sumType_) / totalWeight_.get(sumType_);
		resStream_->newField();
		resStream_->stream() << to<double>(estimate) << std::endl;
		resStream_->endField();
	}
};


template <class SampleSpace, class T, class Mean>
class ArcProbsEstimate : public ISEstimate<SampleSpace, T> {
private:
	Mean sumType_;
	ArcMap<typename Mean::Instance> featureWeights_;
	typename Mean::Instance totalWeight_;
	const SampleSpace& sampleSpace_;
	const Model<T>& model_;
public:
	using ISEstimate<SampleSpace, T>::resStream_;
	
	ArcProbsEstimate(OutStream* resStream,
	                  const SampleSpace& sampleSpace,
	                  const Model<T>& model,
	                  Mean sumType) :
		ISEstimate<SampleSpace, T>(resStream),
		sampleSpace_(sampleSpace),
		model_(model),
		featureWeights_(sampleSpace.n),
		sumType_(sumType)
	{
		reset();
	}
	
	void reset() {
		for (auto arc : Arcs(sampleSpace_.n))
			featureWeights_[arc].reset();
		totalWeight_.reset();
	}
	
	void addSample(const Sample<SampleSpace, T>* sample) {
		T weight = getWeight(sample);
		T mp = sample->getUnnormProb();
		ArcMap<T> fps(sampleSpace_.n);
		model_.template getUnnormArcProbs<SampleSpace>(sample->getInstance(), fps);
		
		for (auto arc : Arcs(sampleSpace_.n)) {
			T featProb = fps[arc] / mp;
			featureWeights_[arc].add(sumType_, weight * featProb);
		}
		totalWeight_.add(sumType_, weight);
	}
	
	void updateCurrentEstimate() {
		resStream_->newField();
		for (auto arc : Arcs(sampleSpace_.n)) {
			T estimate = featureWeights_[arc].get(sumType_) / totalWeight_.get(sumType_);
			resStream_->stream() << to<double>(estimate) << "\n";
		}
		resStream_->endField();
	}
};
*/



#endif

