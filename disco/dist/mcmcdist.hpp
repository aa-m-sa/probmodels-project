/*
 *  BEANDisco: MCMC distribution
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

#ifndef MCMCDIST_HPP
#define MCMCDIST_HPP


/**
 * Markov chain Monte Carlo (MCMC) distribution
 */
template <typename SampleSpace, typename Model, typename T>
class MCMCDist : public ModelDist<SampleSpace, Model, T> {
private:
	using SampleModel = typename ModelDist<SampleSpace, Model, T>::SampleModel;
	//const Model& model_;
	using ModelDist<SampleSpace, Model, T>::nSamples_;
	using ModelDist<SampleSpace, Model, T>::sampleSpace_;
	using ModelDist<SampleSpace, Model, T>::model_;
	using ModelDist<SampleSpace, Model, T>::interruptCallback_;

	std::shared_ptr<MCProposalDist<SampleSpace>> proposalDist_;

	//constexpr size_t interruptTestMask_ = 0xFF;
	#define interruptTestMask_ 0xFF

	typename SampleSpace::Instance instance_;
	//std::unique_ptr<SampleModel> sampleModel_;
	T prob_;

	T uniformInvNormConst;

	int nAccepts_;
	int nSteps_;
	int nStepsPerSample_;
public:
	//using Dist<SampleSpace, T>::nSamples;
	//using Dist<SampleSpace, T>::sampleSpace_;
	
	MCMCDist(
			const SampleSpace& sampleSpace,
			const Model* model,
			std::shared_ptr<MCProposalDist<SampleSpace>> proposalDist,
			//int burnin,
			int stepsPerSample)
	:
		ModelDist<SampleSpace, Model, T>(sampleSpace, model),
		proposalDist_(proposalDist),
		instance_(sampleSpace)
		//model_(model)
	{
		setLog(uniformInvNormConst, -sampleSpace_.logSize());
		nAccepts_ = 0;
		nSteps_ = 0;
		nStepsPerSample_ = stepsPerSample;
		instance_.rand();
		prob_ = model_->template getUnnormProb<SampleSpace>(instance_);
		//sampleModel_ = std::unique_ptr<SampleModel>(
		//		model_->template newRestrictedModel<SampleSpace>(instance_));
		//prob_ = sampleModel_->getUnnormProb();
		//for (int i = 0; i < burnin; ++i)
		//	step();
		nAccepts_ = 0;
		nSteps_ = 0;
	}

	void step() {
		typename SampleSpace::Instance newInstance(sampleSpace_);
		newInstance = instance_;
		//newInstance.randSwap();
		proposalDist_->randStep(newInstance);
		T newProb = model_->template getUnnormProb<SampleSpace>(newInstance);
		//std::unique_ptr<SampleModel> newSampleModel = std::unique_ptr<SampleModel>(
		//		model_->template newRestrictedModel<SampleSpace>(newInstance));
		//T newProb = newSampleModel->getUnnormProb();
		if (randu() < to<double>(newProb / prob_)) {
			instance_ = newInstance;
			//sampleModel_ = std::move(newSampleModel);
			prob_ = newProb;
			++nAccepts_;
		}
		++nSteps_;
	}

	void burn(int steps) {
		for (int i = 0; i < steps; ++i)
			step();
	}
	
	std::unique_ptr<Sample<SampleSpace, Model, T>>
	rand() {
		for (int j = 0; j < nStepsPerSample_; ++j) {
			if (interruptCallback_ && (j & interruptTestMask_) == 0 && interruptCallback_())
				return nullptr;
			step();
		}
		++nSamples_;
		auto modelSample = std::make_shared<ModelSample<SampleSpace, Model, T>>(model_, instance_);

		std::vector<T> invProb(1);
		invProb[0] = uniformInvNormConst / prob_;

		return std::unique_ptr<Sample<SampleSpace, Model, T>>(
				new Sample<SampleSpace, Model, T>(modelSample, 1.0, std::vector<T>(), invProb));
	}
	
	/*Sample<SampleSpace, Model, T>* getSample() {
		Sample<SampleSpace, Model, T>* sample = new Sample<SampleSpace, Model, T>(instance_, prob_);
		//sample->setRestrictedModel(std::move(sampleModel_).get());
		sample->setRestrictedModel(model_->template newRestrictedModel<SampleSpace>(instance_));
		return sample;
	}*/
	
	double acceptRatio() {
		return nAccepts_ / (double) nSteps_;
	}

	int nHarmonicLevels() const {
		return 1;
	}

	virtual std::unique_ptr<ISEstimator<SampleSpace, Model, T, T>>
	newHarmonicUnweightedNormConstEstimator(double meanFraction) const {
		return std::unique_ptr<ISEstimator<SampleSpace, Model, T, T>>(
			new ISMultiLevelHarmonicNormConstEstimator<SampleSpace, Model, T>(1, meanFraction));
	}
};



#endif