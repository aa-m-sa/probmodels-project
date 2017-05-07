/*
 *  BEANDisco: the bucketorder-restricted model
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

#include "common.hpp"
#include "model.hpp"
#include "bucketorder.hpp"

#ifndef POMODEL_HPP
#define POMODEL_HPP


/*template <typename POF, typename T>
class SampleModel<OrderModularModel<T>, POF> {
	typename POF::Instance instance_;
	RestrictedModel<OrderModularModel<T>, POF>
public:
};*/

template <typename POF, typename T>
class RestrictedModel<OrderModularModel<T>, POF> {
public:
	using PredSet = typename POF::PredSet;
	using DownSet = typename POF::DownSet;
	using PSMap = typename POF::template PredSetMap<T>;
	using DSMap = typename POF::template DownSetMap<T>;
private:
	using SPOp = SumProdOp<T>;
	//using RestrictedModel<T>::parentModel_;
	//using RestrictedModel<T>::instance_;
	const OrderModularModel<T>* parentModel_;
	//const ParentSetMap<T>& scores_;
	const POF& family_;
	const typename POF::Instance& instance_;
	//typename POF::Instance instance_;

	mutable PSMap* predSetProbs_;
	mutable PSMap* predSetSums_;
	mutable typename PSMap::ZetaTransformData* predSetPartialSums_;
	mutable DSMap* forwardProbs_;
	mutable DSMap* backwardProbs_;

	template <typename Op>
	void calcPredSetProbs(
			const ParentSetMap<T>* scores,
			PSMap& predSetSums) const
	{
		predSetSums.setAll(0.0);
		for (int node = 0; node < family_.n; ++node) {
			//printf("%d: \n", node);
			SortedArraySubset potPreds(family_.n);
			instance_.getPotPreds(node, potPreds);
			// assume potPreds sorted?
			for (auto subsetScore : scores->forNode(node)->forSubsetsOf(potPreds,
					PredSet(instance_, node))) {
				const auto& predSet = subsetScore.first;
				const auto& score = subsetScore.second;
				Op::add(predSetSums(predSet), score);
				//printf("*");
			}
			//printf("\n");
		}
	}



	template <typename Op, typename Feature>
	void calcPredSetProbs(
			const ParentSetMap<T>* scores,
			const Feature& feature,
			PSMap& predSetSums) const
	{
		predSetSums.setAll(0.0);
		for (int node = 0; node < family_.n; ++node) {
			SortedArraySubset potPreds(family_.n);
			instance_.getPotPreds(node, potPreds);
			// assume potPreds sorted?
			for (auto subsetScore : scores->forNode(node)->forSubsetsOf(potPreds,
					SetPair<PredSet, typename Feature::Subset>(
					PredSet(instance_, node),
					typename Feature::Subset(feature)))) {
				const auto& predSet = subsetScore.first.first;
				const auto& parentSet = subsetScore.first.second;
				const auto& score = subsetScore.second;
				if (feature.holds(node, parentSet))
					Op::add(predSetSums(predSet), score);
			}
		}
	}
	
	/**
	 * Computes gammas from forward and backward sums.
	 */
	template <typename Op>
	void calcForwardBackwardSums(
			const DSMap& forwardProbs,
			const DSMap& backwardProbs,
			PSMap& forwardBackwardSums) const
	{
		// for each variable
		for (int node = 0; node < family_.n; ++node) {
			// iterate over all predecessor sets
			for (auto predSet : instance_.predSets(node)) {
				DownSet downSet(predSet);
				T fp = forwardProbs[downSet];
				downSet.insertWithinBucket(node);
				T bp = backwardProbs[downSet];
				forwardBackwardSums(predSet) = Op::prod(fp, bp);
			}
			forwardBackwardSums.template upZetaTransform<Op>(instance_, node);
		}
	}

	/**
	 * Computes the final (unnormalized) probabilities for each arc from gammas and local scores.
	 */
	template <typename Op>
	void addForwardBackwardSums(
			const PSMap& forwardBackwardSums,
			const ParentSetMap<T>* scores,
			ArcMap<T>& sums) const
	{
		Arc arc;
		for (arc.head = 0; arc.head < family_.n; ++arc.head) {
			SortedArraySubset potPreds(family_.n);
			instance_.getPotPreds(arc.head, potPreds);
			for (auto subsetScore : scores->forNode(arc.head)->forSubsetsOf(potPreds,
					SetPair<SortedArraySubset, PredSet>(SortedArraySubset(family_.n),
					typename POF::PredSet(instance_, arc.head)))) {
				const auto& parentSet = subsetScore.first.first;
				const auto& predSet = subsetScore.first.second;
				const auto& score = subsetScore.second;
				T prod = Op::prod(score, forwardBackwardSums(predSet));
				for (int i = 0; i < parentSet.size(); ++i) {
					arc.tail = parentSet[i];
					Op::add(sums[arc], prod);
				}
			}
		}
	}

public:
	RestrictedModel(
		const OrderModularModel<T>* parentModel,
		const typename POF::Instance& instance)
	:
		parentModel_(parentModel),
		family_(instance.family),
		instance_(instance)
	{
		predSetProbs_ = nullptr;
		predSetSums_ = nullptr;
		predSetPartialSums_ = nullptr;
		forwardProbs_ = nullptr;
		backwardProbs_ = nullptr;
	}

	~RestrictedModel() {
		delete backwardProbs_;
		delete forwardProbs_;
		delete predSetPartialSums_;
		delete predSetSums_;
		delete predSetProbs_;
	}

	const DagFamily& getDagFamily() const {
		return parentModel_->getDagFamily();
	}
	const OrderFamily& getOrderFamily() const {
		return parentModel_->getOrderFamily();
	}
	int getNumVariables() const {
		return parentModel_->getNumVariables();
	}
	
	const OrderModularModel<T>* getParentModel() const {
		return parentModel_;
	}
	
	const POF& getFamily() const {
		return family_;
	}
	const typename POF::Instance& getInstance() const {
		return instance_;
	}


	const PSMap* getPredSetProbs() const {
		if (predSetProbs_ == nullptr) {
			predSetProbs_ = new PSMap(family_);
			calcPredSetProbs<SPOp>(parentModel_->getParentSetFactors(), *predSetProbs_);
		}
		return predSetProbs_;
	}

	const PSMap* getPredSetSums() const {
		if (predSetSums_ == nullptr) {
			predSetSums_ = new PSMap(*getPredSetProbs());
			for (int node = 0; node < family_.n; ++node)
				predSetSums_->template zetaTransform<SPOp>(instance_, node);
		}
		return predSetSums_;
	}

	const typename PSMap::ZetaTransformData* getPredSetPartialSums() const {
		if (predSetPartialSums_ == nullptr) {
			delete predSetSums_;
			predSetSums_ = new PSMap(*getPredSetProbs());
			predSetPartialSums_ = new typename PSMap::ZetaTransformData(*predSetSums_);
			for (int node = 0; node < family_.n; ++node)
				predSetSums_->template zetaTransform<SPOp>(instance_, node, predSetPartialSums_);
		}
		return predSetPartialSums_;
	}

	const DSMap* getForwardProbs() const {
		if (forwardProbs_ == nullptr) {
			forwardProbs_ = new DSMap(family_);
			forwardProbs_->template forwardSum<SPOp>(*getPredSetSums());
		}
		return forwardProbs_;
	}

	const DSMap* getBackwardProbs() const {
		if (backwardProbs_ == nullptr) {
			backwardProbs_ = new DSMap(family_);
			backwardProbs_->template backwardSum<SPOp>(*getPredSetSums());
		}
		return backwardProbs_;
	}
	
	/**
	 * Computes the (unnormalized) probability of given partial order.
	 */
	T getUnnormProb() const {
		return getForwardProbs()->getFull();
	}

	/**
	 * Computes the (unnormalized) probability of given feature in given partial order.
	 */
	template <typename Feature>
	typename Feature::ValueType getUnnormFeatureProb(const Feature& feature) const {
		PSMap predSetProbs(family_);
		calcPredSetProbs<SPOp, Feature>(parentModel_->getParentSetFactors(), feature, predSetProbs);

		// compute alphas
		PSMap predSetSums(predSetProbs);
		for (int node = 0; node < family_.n; ++node)
			predSetSums.template zetaTransform<SPOp>(instance_, node);
		
		// compute the probability
		DSMap forwardProbs(family_);
		forwardProbs.template forwardSum<SPOp>(predSetSums);
		
		return forwardProbs.getFull();
	}

	template <typename Feature>
	typename Feature::ValueType getFeatureProb(const Feature& feature) const {
		return getUnnormFeatureProb(feature) / getUnnormProb();
	}


	/**
	 * Computes the (unnormalized) probabilities of all arc simultaneously in given partial order.
	 */
	void getUnnormArcProbs(ArcMap<T>& probs) const {
		// compute alphas
		//const PSMap* predSetSums = getPredSetSums();
		
		// compute forward and backward functions
		const DSMap* forwardProbs = getForwardProbs();
		const DSMap* backwardProbs = getBackwardProbs();
		//DSMap backwardProbs(family_);
		//backwardProbs.template backwardSum<SPOp>(*predSetSums);
		
		// compute gammas
		PSMap forwardBackwardSums(family_);
		calcForwardBackwardSums<SPOp>(*forwardProbs, *backwardProbs, forwardBackwardSums);
		
		// compute all arc probs at once
		probs.setAll(0.0);
		addForwardBackwardSums<SPOp>(forwardBackwardSums, parentModel_->getParentSetFactors(), probs);
	}

	void getArcProbs(ArcMap<T>& probs) const {
		getUnnormArcProbs(probs);
		probs /= getUnnormProb();
	}
};





template <typename POF, typename T>
class RestrictedModel<RestrictedModel<OrderModularModel<T>, POF>, DagFamily> {
private:
	using ParentModel = RestrictedModel<OrderModularModel<T>, POF>;
	const ParentModel* parentModel_;
	const DagFamily& family_;
	const DagFamily::Instance& instance_;
	
	mutable bool probComputed_;
	mutable T prob_;
public:
	RestrictedModel(
		const ParentModel* parentModel,
		const DagFamily::Instance& instance)
	:
		parentModel_(parentModel),
		family_(instance.family),
		instance_(instance),
		probComputed_(false)
	{
	}

	/*T getUnnormProb() {
		if (!probComputed_) {
			prob_ = T(1.0);
			auto parentSetProbs = parentModel_.getParentsetFactors();
			for (int node = 0; node < parentModel_.n; ++node)
				prob_ *= parentSetProbs(node, instance_.getParents(node));
			probComputed_ = true;
		}
		return prob_;
	}*/

	T getModularUnnormProb() const {
		if (!probComputed_) {
			prob_ = T(1.0);
			auto parentSetProbs = parentModel_->getParentModel()->getParentSetFactors();
			for (int node = 0; node < family_.n; ++node)
				prob_ *= (*parentSetProbs)[node][instance_.getParents(node)];
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

	void getArcProbs(ArcMap<T>& probs) const {
		probs.setAll(0);
		for (int node = 0; node < family_.n; ++node)
			for (auto parent : instance_.getParents(node))
				probs[Arc(parent, node)] = 1;
		//for (auto arc : allArcs(family_.n))
		//	probs[arc] = instance_.contains(arc) ? T(1.0) : T(0.0);
	}

	//void getUnnormArcProbs(ArcMap<T>& probs)
	//{
	//	getArcProbs(probs);
	//	for (auto arc : allArcs(parentModel_.n))
	//		probs[arc] *= getUnnormProb();
	//}
};



//RestrictedModel


/*
template <typename T>
class OrderModularRestrictedModel<T, BucketOrderFamily> {
private:
	OrderModularPartialOrderModel<T, BucketOrderFamily> impl_;
public:
	OrderModularPartialOrderModel(
			const OrderModularModel<T>& parentModel,
			const typename POF::Instance& instance) :
		impl_(parentModel, instance)
	{}
	T getUnnormProb() {
		return impl_.getUnnormProb();
	}

	template <typename Feature>
	T getUnnormFeatureProb(const Feature& feature) const {
		return impl_.getUnnormFeatureProb<Feature>(feature);
	}

	void getUnnormArcProbs(ArcMap<T>& probs)
	{
		impl_.getUnnormArcProbs(probs);
	}
};*/

#endif
