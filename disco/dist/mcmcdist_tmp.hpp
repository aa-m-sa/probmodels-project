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

	std::shared_ptr<ModelSample<SampleSpace, Model, T>> modelSample_;
	//typename SampleSpace::Instance instance_;
	//std::unique_ptr<SampleModel> sampleModel_;
	//T prob_;
	State state_;

	std::vector<T> invProb_;
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
		modelSample_(nullptr),
		instance_(sampleSpace),
		invProb_(1)
		//model_(model)
	{
		setLog(uniformInvNormConst, -sampleSpace_.logSize());
		nAccepts_ = 0;
		nSteps_ = 0;
		nStepsPerSample_ = stepsPerSample;
		instance_.rand();
		//prob_ = model_->template getUnnormProb<SampleSpace>(instance_);
		//invProb_[0] = uniformInvNormConst / prob_;
		nAccepts_ = 0;
		nSteps_ = 0;
		state_ = proposalDist_->initState(instance_);
	}

	void step() {
		//typename SampleSpace::Instance newInstance(sampleSpace_);
		//newInstance = instance_;
		//proposalDist_->randStep(newInstance);
		//T newProb = model_->template getUnnormProb<SampleSpace>(newInstance);
		//if (randu() < to<double>(newProb / prob_)) {
		//	instance_ = newInstance;
		//	prob_ = newProb;
		//	++nAccepts_;
		//}
		auto proposal = proposalDist_->randProposal(state_);
		if (randu() < to<double>(proposal.getMHRatio())) {
			proposal.applyTo(state_);
			++nAccepts_;
		}
		++nSteps_;
	}

	void burn(int steps) {
		for (int i = 0; i < steps; ++i)
			step();
	}
	
	bool rand() {
		for (int j = 0; j < nStepsPerSample_; ++j) {
			if (interruptCallback_ && (j & interruptTestMask_) == 0 && interruptCallback_())
				return false;
			step();
		}
		++nSamples_;
		modelSample_ = std::make_shared<ModelSample<SampleSpace, Model, T>>(model_, state_->getInstance());
		invProb_[0] = uniformInvNormConst / modelSample_->getUnnormProb();

		return true;
	}
	
	/*Sample<SampleSpace, Model, T>* getSample() {
		Sample<SampleSpace, Model, T>* sample = new Sample<SampleSpace, Model, T>(instance_, prob_);
		//sample->setRestrictedModel(std::move(sampleModel_).get());
		sample->setRestrictedModel(model_->template newRestrictedModel<SampleSpace>(instance_));
		return sample;
	}*/
	std::unique_ptr<Sample<SampleSpace, Model, T>>
	getSample() {
		return std::unique_ptr<Sample<SampleSpace, Model, T>>(
				new Sample<SampleSpace, Model, T>(modelSample_, 1.0, std::vector<T>(), invProb_));
	}
	
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