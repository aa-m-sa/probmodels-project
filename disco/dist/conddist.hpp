/*
 *  BEANDisco: conditional distributions for generating DAGs from linear orders and bucket orders
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

#include <vector>

#include "dist.hpp"
#include "../common.hpp"

#include "../model.hpp"

#include "../dag.hpp"
#include "../order.hpp"
#include "../bucketorder.hpp"

#ifndef CONDDIST_HPP
#define CONDDIST_HPP


template <typename SampleSpace, typename Model, typename ConditionFamily, typename T>
class ExactCondDist;


template <typename T>
class ExactCondDist<DagFamily, OrderModularModel<T>, OrderFamily, T> : public ModelDist<DagFamily, RestrictedModel<OrderModularModel<T>, OrderFamily>, T> {
private:
	using Model = RestrictedModel<OrderModularModel<T>, OrderFamily>;
	//using Parent = ModelDist<DagFamily, RestModel, T>;
	using SampleModel = typename ModelDist<DagFamily, Model, T>::SampleModel;
	using ModelDist<DagFamily, Model, T>::nSamples_;
	using ModelDist<DagFamily, Model, T>::sampleSpace_;
	using ModelDist<DagFamily, Model, T>::model_;

	std::vector<DiscreteDist<double>*> parentSetIndexDists_;
	std::vector<std::vector<typename SubsetMap<T>::SubsetIndex>> parentSetIndices_;

public:
	ExactCondDist(const DagFamily& sampleSpace, const Model* model) :
			ModelDist<DagFamily, Model, T>(sampleSpace, model),
			parentSetIndexDists_(sampleSpace.n),
			parentSetIndices_(sampleSpace.n)
	{
		auto parentSetScores = model_->getParentModel()->getParentSetFactors();
		for (int node = 0; node < model_->getNumVariables(); ++node) {
			SortedArraySubset preds(model_->getNumVariables());
			model_->getInstance().getPreds(node, preds);
			parentSetIndices_[node].clear();
			// compute the number of potential parent sets and their total probabilities (scores)
			int nParentSets = 0;
			T totalScore = 0;
			for (const auto& subsetScore : parentSetScores->forNode(node)->forSubsetsOf(preds, DummySet())) {
				const auto& score = subsetScore.second;
				nParentSets += 1;
				totalScore += score;
			}
			// compute and store the (normalized) probabilities of parent sets and store the
			// corresponding parent set indices
			std::vector<double> probs;
			probs.reserve(nParentSets);
			for (const auto& subsetScore : parentSetScores->forNode(node)->forSubsetsOf(preds,
					typename SubsetMap<T>::Subset(*(parentSetScores->forNode(node))))) {
				const auto& parentSet = subsetScore.first;
				const auto& score = subsetScore.second;
				parentSetIndices_[node].push_back(parentSet.getIndex());
				double prob = double(score / totalScore);
				probs.push_back(prob);
			}
			// create a parent set distribution
			parentSetIndexDists_[node] = new DiscreteDist<double>(probs);
		}
	}
	
	~ExactCondDist() {
		for (int node = 0; node < model_->getNumVariables(); ++node) {
			delete parentSetIndexDists_[node];
		}
	}

	std::unique_ptr<Sample<DagFamily, Model, T>> rand() {
		DagFamily::Instance instance(sampleSpace_);
		//T prob = 1.0;
		auto parentSetScores = model_->getParentModel()->getParentSetFactors();
		//auto predSetScores = model_->getParentModel()->getPredSetFactors();
		SortedArraySubset parents(model_->getNumVariables());
		for (int node = 0; node < model_->getNumVariables(); ++node) {
			int i = parentSetIndexDists_[node]->rand();
			auto parentSetIndex = parentSetIndices_[node][i];
			parentSetScores->forNode(node)->getElements(parentSetIndex, parents);
			//parentSet.getElements(parents);
			instance.setParents(node, parents);
			SortedArraySubset preds(model_->getNumVariables());
			model_->getInstance().getPreds(node, preds);
			//prob_ *= (*predSetScores->forNode(node))[preds];
			//prob_ *= (*parentSetScores->forNode(node))[parentSetIndex];
		}
		++nSamples_;
		auto modelSample = std::make_shared<ModelSample<DagFamily, Model, T>>(model_, instance);

		return std::unique_ptr<Sample<DagFamily, Model, T>>(
				new Sample<DagFamily, Model, T>(modelSample));
	}
};



template <typename POF, typename T>
class ExactCondDist<DagFamily, OrderModularModel<T>, POF, T> : public ModelDist<DagFamily, RestrictedModel<OrderModularModel<T>, POF>, T> {
private:
	using Model = RestrictedModel<OrderModularModel<T>, POF>;
	using SampleModel = typename ModelDist<DagFamily, Model, T>::SampleModel;
	using ModelDist<DagFamily, Model, T>::nSamples_;
	using ModelDist<DagFamily, Model, T>::sampleSpace_;
	using ModelDist<DagFamily, Model, T>::model_;

	typename POF::template PredSetMap<T>::ZetaTransformData const* predSetPartialSums_;
	typename POF::template PredSetMap<T> const* predSetSums_;
	typename POF::template DownSetMap<T> const* forwardProbs_;

	struct PredSetDatum {
		int nParentSets;
		T totalProb;
		std::vector<typename SubsetMap<T>::SubsetIndex> parentSetIndices;
		std::vector<double> probs;
		DiscreteDist<double>* dist;
	};
	typename POF::template PredSetMap<PredSetDatum> predSetData_;

public:
	ExactCondDist(const DagFamily& sampleSpace, const Model* model) :
			ModelDist<DagFamily, Model, T>(sampleSpace, model),
			predSetData_(model->getFamily())
	{
		predSetPartialSums_ = model_->getPredSetPartialSums();
		predSetSums_ = model_->getPredSetSums();
		forwardProbs_ = model_->getForwardProbs();

		const auto& poInstance = model_->getInstance();
		auto parentSetScores = model_->getParentModel()->getParentSetFactors();
		for (int node = 0; node < model_->getNumVariables(); ++node) {
			SortedArraySubset potPreds(model_->getNumVariables());
			poInstance.getPotPreds(node, potPreds);
			// init counters
			for (auto predSet : poInstance.predSets(node)) {
				predSetData_[predSet].nParentSets = 0;
				predSetData_[predSet].totalProb = 0;
			}
			// compute the number of parent sets and total probabilities (scores) for each potential
			// predecessor set
			for (const auto& subsetScore : parentSetScores->forNode(node)->forSubsetsOf(potPreds,
					typename POF::PredSet(poInstance, node))) {
				const auto& predSet = subsetScore.first;
				const auto& score = subsetScore.second;
				predSetData_[predSet].nParentSets += 1;
				predSetData_[predSet].totalProb += score;
			}
			// reserve required space
			for (auto predSet : poInstance.predSets(node)) {
				predSetData_[predSet].parentSetIndices.reserve(predSetData_[predSet].nParentSets);
				predSetData_[predSet].probs.reserve(predSetData_[predSet].nParentSets);
			}
			// compute and store (the normalized) probabilities for parent sets and store the
			// corresponding parent set indices
			for (const auto& subsetScore : parentSetScores->forNode(node)->forSubsetsOf(potPreds,
					SetPair<typename SubsetMap<T>::Subset, typename POF::PredSet>(
						typename SubsetMap<T>::Subset(*(parentSetScores->forNode(node))),
						typename POF::PredSet(poInstance, node)))) {
				const auto& parentSet = subsetScore.first.first;
				const auto& predSet = subsetScore.first.second;
				const auto& score = subsetScore.second;
				predSetData_[predSet].parentSetIndices.push_back(parentSet.getIndex());
				double prob = (double)(score / predSetData_[predSet].totalProb);
				predSetData_[predSet].probs.push_back(prob);
			}
			// create a parent set distribution for each potential predecessor set
			for (auto predSet : poInstance.predSets(node)) {
				predSetData_[predSet].dist = new DiscreteDist<double>(predSetData_[predSet].probs);
			}
		}
	}
	
	~ExactCondDist() {
		const auto& poInstance = model_->getInstance();
		for (int node = 0; node < model_->getNumVariables(); ++node) {
			for (auto predSet : poInstance.predSets(node)) {
				delete predSetData_[predSet].dist;
			}
		}
	}

	std::unique_ptr<Sample<DagFamily, Model, T>> rand() {
		// draw the order
		const auto& poInstance = model_->getInstance();
		std::vector<int> order(sampleSpace_.n);
		forwardProbs_->template backTrackForwardSumOrder<RandSelector<T>>(*predSetSums_, poInstance, order);
		//std::cout << order << std::endl;

		// draw the parents
		//auto predSetScores = model_->getParentModel()->getPredSetFactors();
		auto parentSetScores = model_->getParentModel()->getParentSetFactors();
		SortedArraySubset parents(model_->getNumVariables());
		//typename POF::DownSet downSet = poInstance.getFullDownSet();
		//for (int i = model_->getNumVariables(); i >= 0; --i) {
		//	int node = order[i];
			//downSet.remove(node);
		DagFamily::Instance instance(sampleSpace_);
		typename POF::DownSet downSet(poInstance);
		//prob_ = 1.0;
		for (int node : order) {
			typename POF::PredSet predSet(node, downSet);
			typename POF::PredSet tailSet(poInstance, node);
			predSetSums_->template backTrackZetaTransformSet<RandSelector<T>>(
					poInstance, node, predSetPartialSums_, predSet, tailSet);
			int i = predSetData_[tailSet].dist->rand();
			assert(0 <= i && i < predSetData_[tailSet].parentSetIndices.size());
			auto parentSetIndex = predSetData_[tailSet].parentSetIndices[i];
			parentSetScores->forNode(node)->getElements(parentSetIndex, parents);
			//parentSet.getElements(parents);
			//std::cout << "add: " << node << " <- " << parents << " ⊆ " << downSet.getElements() << std::endl;
			assert(isSubsetOf(parents, downSet.getElements()));
			instance.setParents(node, parents);
			//prob_ *= (*predSetScores->forNode(node))[predSet];
			//prob_ *= (*parentSetScores->forNode(node))[parentSetIndex];
			downSet.insert(node);
		}
		//RestrictedModel<T, OrderFamily, SumProdOp<T>> orderModel();
		//ExactOCondDagDist<T> oCondDagDist(dagFamily, orderModel);
		//model_->get
		//T sampleSpaceSize;
		//setLog(sampleSpaceSize, model_.condition_.family.logNumLinExts());
		//prob_ = T(1.0) / sampleSpaceSize;
		auto modelSample = std::make_shared<ModelSample<DagFamily, Model, T>>(model_, instance);
		++nSamples_;

		return std::unique_ptr<Sample<DagFamily, Model, T>>(
				new Sample<DagFamily, Model, T>(modelSample));
	}
};



/*template <typename SampleSpace, typename Model, typename T>
class ISDist : public ModelDist<SampleSpace, Model, T>{
	Sample<SampleSpace, T>* getSample() {
		new Sample<DagFamily, T>(instance_, prob_, weight_);
	}
};*/



#endif
