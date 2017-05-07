/*
 *  BEANDisco: AIS distribution
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


#include "dist.hpp"
#include "../common.hpp"

#ifndef AISDIST_HPP
#define AISDIST_HPP

/**
 * Annealed importance sampling (AIS) distribution
 */
template <typename SampleSpace, typename Model, typename T>
class AnnealedDist : public ModelDist<SampleSpace, Model, T> {
private:
	using SampleModel = typename ModelDist<SampleSpace, Model, T>::SampleModel;
	//const Model& model_;
	using ModelDist<SampleSpace, Model, T>::nSamples_;
	using ModelDist<SampleSpace, Model, T>::sampleSpace_;
	using ModelDist<SampleSpace, Model, T>::model_;
	using ModelDist<SampleSpace, Model, T>::interruptCallback_;

	std::shared_ptr<MCProposalDist<SampleSpace>> proposalDist_;
	const std::vector<double> powers_;
	const size_t samplesPerAnnealing_;
	const size_t stepsBetweenSamples_;

	//constexpr size_t interruptTestMask_ = 0xFF;
	#define interruptTestMask_ 0xFF

	std::vector<long long int> nAccepts_;
	long long int nCoolAccepts_;

	size_t samplesDrawnCurrent_;
	typename SampleSpace::Instance instance_;
	//SampleModel* sampleModel_;
	T prob_;
	T weight_;
public:
	
	AnnealedDist(
			const SampleSpace& sampleSpace,
			const Model* model,
			std::vector<double> powers,
			std::shared_ptr<MCProposalDist<SampleSpace>> proposalDist,
			size_t samplesPerAnnealing = 1,
			size_t stepsBetweenSamples = 1)
	:
		ModelDist<SampleSpace, Model, T>(sampleSpace, model),
		//model_(model),
		proposalDist_(proposalDist),
		powers_(powers),
		samplesPerAnnealing_(samplesPerAnnealing),
		stepsBetweenSamples_(stepsBetweenSamples),
		nAccepts_(powers.size()),
		samplesDrawnCurrent_(samplesPerAnnealing),
		instance_(sampleSpace)
	{
		for (int i = 0; i < nAccepts_.size(); ++i)
			nAccepts_[i] = 0;
		nCoolAccepts_ = 0;
	}

	std::unique_ptr<Sample<SampleSpace, Model, T>>
	rand() {
		typename SampleSpace::Instance newInstance(sampleSpace_);

		// new annealing
		if (samplesDrawnCurrent_ >= samplesPerAnnealing_) {
			T invWeight;
			setLog(invWeight, -sampleSpace_.logSize());
			instance_.rand();
			//std::unique_ptr<ModelSample<SampleSpace, Model, T>> modelSample_(
			//		new ModelSample<SampleSpace, Model, T>(model_, instance));
			prob_ = model_->template getUnnormProb<SampleSpace>(instance_);
			//SampleModel* oldSampleModel = model_->template newRestrictedModel<SampleSpace>(oldInstance);
			//T oldProb = oldSampleModel->getUnnormProb();
			for (size_t i = 0; i < powers_.size(); ++i) {
				if (interruptCallback_ && (i & interruptTestMask_) == 0 && interruptCallback_())
					return nullptr;
				newInstance = instance_;
				proposalDist_->randStep(newInstance);
				T newProb = model_->template getUnnormProb<SampleSpace>(newInstance);
				//SampleModel* newSampleModel = model_->template newRestrictedModel<SampleSpace>(oldInstance);
				//T newProb = newSampleModel->getUnnormProb();
				T relProb = pow(newProb / prob_, powers_[i]);
				if (randu() < to<double>(relProb)) {
					instance_ = newInstance;
					prob_ = newProb;
					invWeight *= relProb;
					++nAccepts_[i];
				}
			}
			weight_ = prob_ / invWeight;
			samplesDrawnCurrent_ = 1;
		}
		else {
			for (size_t i = 0; i < stepsBetweenSamples_; ++i) {
				if (interruptCallback_ && (i & interruptTestMask_) == 0 && interruptCallback_())
					return nullptr;
				newInstance = instance_;
				proposalDist_->randStep(newInstance);
				T newProb = model_->template getUnnormProb<SampleSpace>(newInstance);
				T relProb = newProb / prob_;
				if (randu() < to<double>(relProb)) {
					instance_ = newInstance;
					prob_ = newProb;
					++nCoolAccepts_;
				}
			}
			++samplesDrawnCurrent_;
		}
		auto modelSample = std::make_shared<ModelSample<SampleSpace, Model, T>>(model_, instance_);
		++nSamples_;

		return std::unique_ptr<Sample<SampleSpace, Model, T>>(
				new Sample<SampleSpace, Model, T>(modelSample, weight_));
	}
	
	double acceptRatio(int level) {
		return nAccepts_[level] / (double) nSamples_;
	}

	virtual std::unique_ptr<ISEstimator<SampleSpace, Model, T, T>>
	newDirectUnweightedNormConstEstimator(double meanFraction) const {
		return std::unique_ptr<ISEstimator<SampleSpace, Model, T, T>>(
			new ISMultiLevelNormConstEstimator<SampleSpace, Model, T>(1, meanFraction));
	}
};



#endif