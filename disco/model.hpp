/*
 *  BEANDisco: the model
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

#include "data.hpp"
#include "scores.hpp"
#include "parentsetmap.hpp"
//#include "stacksubset.hpp"
#include "lognum.hpp"

//#include "sample.hpp"

#include "bucketorder.hpp"
#include "order.hpp"
#include "dag.hpp"

#include "feature.hpp"

#ifndef MODEL_HPP
#define MODEL_HPP


/**
 * Computes the scores for all node-parentset pairs
 */
/*template <class T>
void computeScores(const DataView* dataView, const ScoreFun* scoreFun, ParentSetMap<T>& scores) {
	for (int node = 0; node < scores.nNodes; ++node) {
		for (auto& parentsetScore : scores[node]->forSubsets(SortedArraySubset(scores.maxParents))) {
			const auto& parents = parentsetScore.first;
			auto& score = parentsetScore.second;
			double logScore = computeScore(dataView, parents, node, scoreFun);
			setLog(score, logScore);
		}
	}
}*/
/*template <class T>
void computeScores(const DataView* dataView, const ScoreFun* scoreFun, ParentSetMap<T>& scores) {
	StackSubset parents(scores.maxParents);
	for (int node = 0; node < scores.nNodes; ++node) {
		parents.clear();
		do {
			if (parents.contains(node))
				continue;
			double logScore = computeScore(dataView, parents, node, scoreFun);
			setLog(scores(node, parents), logScore);
		} while (parents.next(0, scores.nNodes, scores.maxParents));
	}
}*/


template <typename T>
std::unique_ptr<ParentSetMap<T>> buildScoresFromData(const DataView* dataView, const ScoreFun* scoreFun,
		int maxParents) {
	ParentSetMap<T>* scores = new ParentSetMap<T>(dataView->getNumVariables(), maxParents);
	for (int node = 0; node < scores->nNodes; ++node) {
		for (auto parentsetScore : scores->forNode(node)->forSubsets(SortedArraySubset(scores->maxParents))) {
			const auto& parents = parentsetScore.first;
			auto& score = parentsetScore.second;
			double logScore = computeScore(dataView, parents, node, scoreFun);
			setLog(score, logScore);
		}
	}
	return std::unique_ptr<ParentSetMap<T>>(scores);
}


/**
 * Divides scores by the number of parent sets with same size.
 */
template <class T>
void applyIndegreeUniformPrior(ParentSetMap<T>* scores) {
	auto n = scores->nNodes;

	// precompute inverse binomial coefficients
	std::vector<T> invBinoms(n - 1);
	for (int k = 0; k < n - 1; ++k)
		invBinoms[k] = T(1.0) / binom<T>(n - 1, k);

	// apply the prior
	for (int node = 0; node < n; ++node) {
		for (auto subsetScore : scores->forNode(node)->forSubsets(
				SortedArraySubset(scores->maxParents))) {
			const auto& parentSet = subsetScore.first;
			auto& score = subsetScore.second;
			score *= invBinoms[parentSet.size()];
		}
	}
}

//template <typename T>
//class StructurePrior {
//protected:
//	const DagFamily& dagFamily_;
//public:
//	//virtual T getUnnormProb(DagFamily::Instance dag) = 0;
//};
//
//template <typename T>
//class ModularStructurePrior : public StructurePrior<T> {
//public:
//	/*virtual T getUnnormProb(const DagFamily::Instance& dag) {
//		T unnormProb = 1.0;
//		StackSubset parentset;
//		for (int v = 0; v < dagFamily_.n; ++v) {
//			dag.getParents(v, parentset);
//			unnormProb *= getParentsetFactor(v, parentset);
//		}
//	}*/
//	virtual T getParentsetFactor(int node, const StackSubset& parentset) const = 0;
//};
//
//template <typename T>
//class UniformModularStructurePrior : public ModularStructurePrior<T> {
//	T getParentsetFactor(int node, const StackSubset& parentset) const {
//		return T(1.0);
//	}
//};
//
//template <typename T>
//class OrderPrior {
//protected:
//	const OrderFamily& orderFamily_;
//public:
//	//virtual T getUnnormProb(OrderFamily::Instance order) = 0;
//};
//
//template <typename T>
//class ModularOrderPrior : public OrderPrior<T> {
//public:
//	/*virtual T getUnnormProb(const OrderFamily::Instance& order) {
//		T unnormProb = 1.0;
//		StackSubset predSet;
//		for (int v = 0; v < orderFamily_.n; ++v) {
//			order.getPredSet(node, predSet);
//			unnormProb *= getPredSetFactor(v, predSet);
//		}
//	}*/
//	virtual T getPredSetFactor(int node, const StackSubset& predSet) const = 0;
//};
//
//template <typename T>
//class UniformModularOrderPrior : public ModularOrderPrior<T> {
//	T getPredSetFactor(int node, const StackSubset& predSet) const {
//		return T(1.0);
//	}
//};
//
//template <typename T>
//class OrderModularStructurePrior : public StructurePrior<T> {
//public:
//	OrderModularStructurePrior(ModularStructurePrior<T> _structurePrior, ModularOrderPrior<T> _orderPrior)
//		: structurePrior(_structurePrior), orderPrior(_orderPrior)
//	{}
//	ModularStructurePrior<T> structurePrior;
//	ModularOrderPrior<T> orderPrior;
//};
//
//template <typename T>
//class UniformOrderModularStructurePrior : public OrderModularStructurePrior<T> {
//	UniformOrderModularStructurePrior()
//		: OrderModularStructurePrior<T>(UniformModularStructurePrior<T>(), UniformModularStructurePrior<T>())
//	{}
//};




template <typename T>
class PredSetMap {
};

template <typename T>
class StructureModel {
protected:
	const DagFamily& dagFamily_;
public:
	StructureModel(const DagFamily& dagFamily) :
		dagFamily_(dagFamily)
	{}

	virtual ~StructureModel()
	{}

	const DagFamily& getDagFamily() const {
		return dagFamily_;
	}

	int getNumVariables() const {
		return dagFamily_.n;
	}

};


template <typename OuterModel, typename SampleSpace>
class RestrictedModel;

template <typename T>
class ModularModel : public StructureModel<T> {
private:
	using This = ModularModel<T>;
	const ParentSetMap<T>* parentSetFactors_;
public:
	ModularModel(
			const DagFamily& dagFamily,
			const ParentSetMap<T>* parentSetFactors) :
		StructureModel<T>(dagFamily),
		parentSetFactors_(parentSetFactors)
	{
	}

	const ParentSetMap<T>* getParentSetFactors() const {
		return parentSetFactors_;
	}

	template <typename SampleSpace>
	std::unique_ptr<RestrictedModel<This, SampleSpace>>
	newRestrictedModel(const typename SampleSpace::Instance& instance) const {
		return std::unique_ptr<RestrictedModel<This, SampleSpace>>(
				new RestrictedModel<This, SampleSpace>(this, instance));
	}

	template <typename SampleSpace>
	T getUnnormProb(const typename SampleSpace::Instance instance) const {
		auto restrictedModel = newRestrictedModel<SampleSpace>(instance);
		return restrictedModel->getUnnormProb();
	}
};



template <typename T>
class OrderModularModel : public StructureModel<T> {
protected:
	using StructureModel<T>::dagFamily_;
	const OrderFamily& orderFamily_;
private:
	using This = OrderModularModel<T>;
	const ParentSetMap<T>* parentSetFactors_;
	const PredSetMap<T>* predSetFactors_;
public:
	OrderModularModel(
			const DagFamily& dagFamily,
			const OrderFamily& orderFamily,
			const ParentSetMap<T>* parentSetFactors,
			const PredSetMap<T>* predSetFactors) :
		StructureModel<T>(dagFamily),
		orderFamily_(orderFamily),
		parentSetFactors_(parentSetFactors),
		predSetFactors_(predSetFactors)
	{
		assert(dagFamily_.n == orderFamily_.n);
	}

	const OrderFamily& getOrderFamily() const {
		return orderFamily_;
	}

	const ParentSetMap<T>* getParentSetFactors() const {
		return parentSetFactors_;
	}

	const PredSetMap<T>* getPredSetFactors() const {
		return predSetFactors_;
	}

	template <typename SampleSpace>
	std::unique_ptr<RestrictedModel<This, SampleSpace>>
	newRestrictedModel(const typename SampleSpace::Instance& instance) const {
		return std::unique_ptr<RestrictedModel<This, SampleSpace>>(
				new RestrictedModel<This, SampleSpace>(this, instance));
	}

	template <typename SampleSpace>
	T getUnnormProb(const typename SampleSpace::Instance instance) const {
		auto restrictedModel = newRestrictedModel<SampleSpace>(instance);
		return restrictedModel->getUnnormProb();
	}
};



/*
template <typename T>
class RestrictedModel<ModularModel<T>, DagFamily> {
private:
	const ModularModel<T>& parentModel_;
	const DagFamily::Instance& instance_;
	
	bool probComputed_;
	T prob_;
public:
	RestrictedModel(
		const ModularModel<T>& parentModel,
		const DagFamily::Instance& instance)
	:
		parentModel_(parentModel),
		instance_(instance),
		probComputed_(false)
	{
	}

	T getUnnormProb() {
		if (!probComputed_) {
			prob_ = T(1.0);
			auto psProbs = parentModel_.getParentsetFactors();
			for (int node = 0; node < parentModel_.n; ++node)
				prob_ *= psProbs(node, instance_.getParents(node));
			probComputed_ = true;
		}
		return prob_;
	}

	template <typename Feature>
	typename Feature::ValueType getFeatureProb(const Feature& feature) const {
		return feature.value(instance_);
	}

	//template <typename Feature>
	//T getUnnormFeatureProb(const Feature& feature) const {
	//	return getUnnormProb() * feature.value(instance_);
	//}

	void getArcProbs(ArcMap<T>& probs) {
		probs.setAll(0);
		for (int node = 0; node < family_.n; ++node)
			for (auto parent : instance_.getParents(node))
				probs[Arc(parent, node)] = 1;
		//for (auto arc : allArcs(parentModel_.n))
		//	probs[arc] = instance_.contains(arc) ? T(1.0) : T(0.0);
	}

	//void getUnnormArcProbs(ArcMap<T>& probs)
	//{
	//	getArcProbs(probs);
	//	for (auto arc : allArcs(parentModel_.n))
	//		probs[arc] *= getUnnormProb();
	//}
};*/




#endif

