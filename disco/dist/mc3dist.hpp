/*
 *  BEANDisco: MC3 distribution
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

#ifndef MC3DIST_HPP
#define MC3DIST_HPP


/**
 * Metropolis coupled Markov chain Monte Carlo (MC³) sampler
 */
template <typename SampleSpace, typename Model, typename T>
class MC3Sampler {
private:
	using SampleModel = RestrictedModel<Model, SampleSpace>;

	const Model* model_;
	const SampleSpace& sampleSpace_;
	
	std::shared_ptr<MCProposalDist<SampleSpace>> proposalDist_;

	int nLevels_;
	int nSwapsPerStep_;
	std::vector<double> levelExponents_;
	std::vector<typename SampleSpace::Instance> states_;
	std::vector<T> probs_;
	std::vector<int> nMHAccepts_;
	std::vector<int> nMHSteps_;
	std::vector<int> nSAccepts_;
	std::vector<int> nSSteps_;
	
public:
	MC3Sampler(
			const Model* model,
			const SampleSpace& sampleSpace,
			const std::vector<double>& levelExponents,
			int nSwapsPerStep,
			std::shared_ptr<MCProposalDist<SampleSpace>> proposalDist)
	:
		model_(model),
		sampleSpace_(sampleSpace),
		proposalDist_(proposalDist),
		nLevels_(levelExponents.size()), 
		nSwapsPerStep_(nSwapsPerStep),
		levelExponents_(levelExponents),
		states_(nLevels_ + 1, sampleSpace),
		probs_(nLevels_ + 1),
		nMHAccepts_(nLevels_ + 1),
		nMHSteps_(nLevels_ + 1),
		nSAccepts_(nLevels_ + 1),
		nSSteps_(nLevels_ + 1)
	{
		assert(levelExponents.size() >= 1);
		levelExponents_.push_back(1.0);
		for (int l = 0; l <= nLevels_; ++l) {
			states_[l].rand();
			//SampleModel* sampleModel = model_->template newRestrictedModel<SampleSpace>(states_[l]);
			//probs_[l] = sampleModel->getUnnormProb();
			//delete sampleModel;
			probs_[l] = model_->template getUnnormProb<SampleSpace>(states_[l]);
		}
		
		resetStats();
	}
	
	~MC3Sampler() {
	}

	void resetStats() {
		for (int l = 0; l <= nLevels_; ++l) {
			nMHAccepts_[l] = 0;
			nMHSteps_[l] = 0;
			nSAccepts_[l] = 0;
			nSSteps_[l] = 0;
		}
	}
	
	int nLevels() const {
		return nLevels_;
	}
	
	double getMHAcceptRatio(int level) const {
		return nMHAccepts_[level] / (double) nMHSteps_[level];
	}

	double getSAcceptRatio(int level) const {
		return nSAccepts_[level] / (double) nSSteps_[level];
	}
	
	T getLevelUnnormProbRatio(int sampleLevel, int otherLevel) const {
		return pow(probs_[sampleLevel],
				levelExponents_[otherLevel] - levelExponents_[sampleLevel]);
	}
	
	T getLevelUnnormProb(int level) {
		return pow(probs_[level], levelExponents_[level]);
	}
	
	const typename SampleSpace::Instance& getCurrentState(int level = 0) const {
		return states_[level];
	}
	
	void mcmcStep(int level, int nSwaps = 1) {
		typename SampleSpace::Instance newState(sampleSpace_);
		newState = states_[level];
		for (int i = 0; i < nSwaps; ++i) {
			proposalDist_->randStep(newState);
			//newState.randSwap();
		}

		//SampleModel* newSampleModel = model_->template newRestrictedModel<SampleSpace>(newState);
		//newProb = newSampleModel->getUnnormProb();
		//delete newSampleModel;
		T newProb = model_->template getUnnormProb<SampleSpace>(newState);
		T probRatio = pow(newProb / probs_[level], levelExponents_[level]);
		if (randu() < to<double>(probRatio)) {
			states_[level] = newState;
			probs_[level] = newProb;
			++nMHAccepts_[level];
		}
		++nMHSteps_[level];
	}
	
	void swapStep(int level) {
		T probRatio = probs_[level + 1] / probs_[level];
		probRatio = pow(probRatio, levelExponents_[level] - levelExponents_[level + 1]);
		if (randu() < to<double>(probRatio)) {
			std::swap(states_[level], states_[level + 1]);
			std::swap(probs_[level], probs_[level + 1]);
			++nSAccepts_[level];
		}
		++nSSteps_[level];
	}
	
	void step() {
		for (int l = 0; l <= nLevels_; ++l) {
			mcmcStep(l);
		}
		for (int i = 0; i < nSwapsPerStep_; ++i) {
			int l = int(randu() * (nLevels_));
			l = (l == nLevels_ ? nLevels_ - 1 : l);
			swapStep(l);
		}
	}
};



/**
 * Metropolis coupled Markov chain Monte Carlo (MC³) distribution
 */
template <typename SampleSpace, typename Model, typename T>
class MC3Dist : public MultiLevelModelDist<SampleSpace, Model, T> {
private:
	//const Model& model_;
	using MultiLevelModelDist<SampleSpace, Model, T>::nSamples_;
	using MultiLevelModelDist<SampleSpace, Model, T>::sampleSpace_;
	using MultiLevelModelDist<SampleSpace, Model, T>::model_;
	using MultiLevelModelDist<SampleSpace, Model, T>::nLevels_;
	using MultiLevelModelDist<SampleSpace, Model, T>::interruptCallback_;

	MC3Sampler<SampleSpace, Model, T>* sampler_;
	int stepsPerSample_;
	
	std::vector<T> levelProbRatios;
	std::vector<T> invLevelProbRatios;

	T uniformNormConst;
public:
	MC3Dist(
			const SampleSpace& sampleSpace,
			const Model* model,
			std::vector<double> levelExponents,
			std::shared_ptr<MCProposalDist<SampleSpace>> proposalDist,
			//int burnInSteps,
			int stepsPerSample,
			int levelSwapsPerStep)
	:
		MultiLevelModelDist<SampleSpace, Model, T>(sampleSpace, model, levelExponents.size() + 1, levelExponents.size() + 1),
		levelProbRatios(nLevels_),
		invLevelProbRatios(nLevels_)
		//model_(model)
	{
		setLog(uniformNormConst, sampleSpace_.logSize());
		levelExponents.insert(levelExponents.begin(), 0.0);
		sampler_ = new MC3Sampler<SampleSpace, Model, T>(model, sampleSpace, levelExponents, levelSwapsPerStep, proposalDist);
		//for (int i = 0; i < burnInSteps; ++i)
		//	sampler_->step();
		stepsPerSample_ = stepsPerSample;
		sampler_->resetStats();
	}
	
	~MC3Dist() {
		delete sampler_;
	}

	void burn(int steps) {
		for (int i = 0; i < steps; ++i)
			sampler_->step();
	}
	
	std::unique_ptr<Sample<SampleSpace, Model, T>>
	rand() {
		for (int i = 0; i < stepsPerSample_; ++i) {
			if (interruptCallback_ && interruptCallback_())
				return nullptr;
			sampler_->step();
		}
		++nSamples_;
		auto instance = sampler_->getCurrentState(nLevels_);
		auto modelSample = std::make_shared<ModelSample<SampleSpace, Model, T>>(model_, instance);

		for (int l = 0; l < nLevels_; ++l) {
			levelProbRatios[l] = sampler_->getLevelUnnormProbRatio(l, l + 1);
			invLevelProbRatios[l] = sampler_->getLevelUnnormProbRatio(l + 1, l);
		}
		levelProbRatios[0] *= uniformNormConst;
		invLevelProbRatios[0] /= uniformNormConst;

		return std::unique_ptr<Sample<SampleSpace, Model, T>>(
				new Sample<SampleSpace, Model, T>(modelSample, 1.0, levelProbRatios, invLevelProbRatios));
	}

	/*void getLevelRatios(std::vector<T>& levelRatios) {
		assert(levelRatios.size() == nLevels + 1);
		levelRatios[0] = sampleSpace_.logSize();
		for (int l = 0; l < nLevels; ++l)
			levelRatios[l + 1] = sampler_->getLevelUnnormProbRatio(l, l + 1);
	}*/
	
	double acceptRatio(int level) {
		return sampler_->getMHAcceptRatio(level);
	}
	double swapAcceptRatio(int level) {
		return sampler_->getSAcceptRatio(level);
	}

	virtual std::unique_ptr<ISEstimator<SampleSpace, Model, T, T>>
	newDirectUnweightedNormConstEstimator(double meanFraction) const {
		return std::unique_ptr<ISEstimator<SampleSpace, Model, T, T>>(
			new ISMultiLevelNormConstEstimator<SampleSpace, Model, T>(nLevels_, meanFraction));
	}

	virtual std::unique_ptr<ISEstimator<SampleSpace, Model, T, T>>
	newHarmonicUnweightedNormConstEstimator(double meanFraction) const {
		return std::unique_ptr<ISEstimator<SampleSpace, Model, T, T>>(
			new ISMultiLevelHarmonicNormConstEstimator<SampleSpace, Model, T>(nLevels_, meanFraction));
	}
};




#endif