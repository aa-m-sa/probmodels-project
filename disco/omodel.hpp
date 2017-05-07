/*
 *  BEANDisco: the order-restricted model
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
#include "order.hpp"

#ifndef OMODEL_HPP
#define OMODEL_HPP

template <typename T>
class RestrictedModel<OrderModularModel<T>, OrderFamily> {
private:
	using SPOp = SumProdOp<T>;
	const OrderModularModel<T>* parentModel_;
	const OrderFamily& family_;
	const OrderFamily::Instance& instance_;

	mutable T* predSetProbs_;
	mutable T* forwardProbs_;
	mutable T* backwardProbs_;

	template <typename Op>
	void calcPredSetProbs(
			const ParentSetMap<T>* scores,
			T* predSetSums) const
	{
		for (int node = 0; node < family_.n; ++node) {
			predSetSums[node] = 0.0;
			SortedArraySubset potPreds(family_.n);
			instance_.getPreds(node, potPreds);
			for (const auto& subsetScore : scores->forNode(node)->forSubsetsOf(potPreds, DummySet())) {
				const auto& score = subsetScore.second;
				Op::add(predSetSums[node], score);
			}
		}
	}



	template <typename Op, typename Feature>
	void calcPredSetProbs(
			const ParentSetMap<T>* scores,
			const Feature& feature,
			T* predSetSums) const
	{
		for (int node = 0; node < family_.n; ++node) {
			predSetSums[node] = 0.0;
			SortedArraySubset potPreds(family_.n);
			instance_.getPreds(node, potPreds);
			for (const auto& subsetScore : scores->forNode(node)->forSubsetsOf(potPreds,
					typename Feature::Subset(feature))) {
				const auto& parentSet = subsetScore.first;
				const auto& score = subsetScore.second;
				if (feature.holds(node, parentSet))
					Op::add(predSetSums[node], score);
			}
		}
	}
	
	/**
	 * Computes gammas from forward and backward sums.
	 */
	template <typename Op>
	void calcForwardBackwardSums(
			const T* forwardProbs,
			const T* backwardProbs,
			T* forwardBackwardSums) const
	{
		// for each variable
		for (int node = 0; node < family_.n; ++node) {
			forwardBackwardSums[node] = Op::prod(forwardProbs[node], backwardProbs[node + 1]);
		}
	}

	/**
	 * Computes the final (unnormalized) probabilities for each arc from gammas and local scores.
	 */
	template <typename Op>
	void addForwardBackwardSums(
			const T* forwardBackwardSums,
			const ParentSetMap<T>* scores,
			ArcMap<T>& sums) const
	{
		Arc arc;
		for (arc.head = 0; arc.head < family_.n; ++arc.head) {
			SortedArraySubset preds(family_.n);
			instance_.getPreds(arc.head, preds);
			for (const auto& subsetScore : scores->forNode(arc.head)->forSubsetsOf(preds, SortedArraySubset(family_.n))) {
				const auto& parentSet = subsetScore.first;
				const auto& score = subsetScore.second;
				T prod = Op::prod(score, forwardBackwardSums[arc.head]);
				for (int i = 0; i < parentSet.size(); ++i) {
					arc.tail = parentSet[i];
					Op::add(sums[arc], prod);
				}
			}
		}
	}

	T* getPredSetProbs() const {
		if (predSetProbs_ == nullptr) {
			predSetProbs_ = new T[family_.n];
			calcPredSetProbs<SPOp>(parentModel_->getParentSetFactors(), predSetProbs_);
		}
		return predSetProbs_;
	}



	T* getForwardProbs() const {
		if (forwardProbs_ == nullptr) {
			T* predSetProbs = getPredSetProbs();
			forwardProbs_ = new T[family_.n + 1];
			forwardProbs_[0] = 1.0;
			for (int node = 0; node < family_.n; ++node)
				forwardProbs_[node + 1] = SPOp::prod(forwardProbs_[node], predSetProbs[node]);
		}
		return forwardProbs_;
	}

	T* getBackwardProbs() const {
		if (backwardProbs_ == nullptr) {
			T* predSetProbs = getPredSetProbs();
			backwardProbs_ = new T[family_.n + 1];
			backwardProbs_[family_.n] = 1.0;
			for (int node = family_.n - 1; node > 0; --node)
				backwardProbs_[node] = SPOp::prod(backwardProbs_[node + 1], predSetProbs[node]);
		}
		return backwardProbs_;
	}

	RestrictedModel(const RestrictedModel& other) = delete;
	RestrictedModel& operator=(const RestrictedModel& other) = delete;
public:
	RestrictedModel(
		const OrderModularModel<T>* parentModel,
		const OrderFamily::Instance& instance)
	:
		parentModel_(parentModel),
		family_(instance.family),
		instance_(instance)
	{
		predSetProbs_ = nullptr;
		forwardProbs_ = nullptr;
		backwardProbs_ = nullptr;
	}

	~RestrictedModel() {
		delete[] predSetProbs_;
		delete[] forwardProbs_;
		delete[] backwardProbs_;
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

	const OrderFamily& getFamily() const {
		return family_;
	}
	const OrderFamily::Instance& getInstance() const {
		return instance_;
	}
	
	/**
	 * Computes the (unnormalized) probability of given partial order.
	 */
	T getUnnormProb() const {
		return getForwardProbs()[family_.n];
	}

	/**
	 * Computes the (unnormalized) probability of given feature in given partial order.
	 */
	template <typename Feature>
	typename Feature::ValueType getUnnormFeatureProb(const Feature& feature) const {
		T* predSetProbs = new T[family_.n];
		calcPredSetProbs<SPOp, Feature>(parentModel_->getParentSetFactors(), feature, predSetProbs);

		// compute the probability
		T prob = 1.0;
		for (int node = 0; node < family_.n; ++node)
			prob = SPOp::prod(prob, predSetProbs[node]);

		delete[] predSetProbs;
		return prob;
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
		//T* predSetProbs = getPredSetProbs();
		
		// compute forward and backward functions
		T* forwardProbs = getForwardProbs();
		T* backwardProbs = getBackwardProbs();

		
		// compute gammas
		T* forwardBackwardSums = new T[family_.n];
		calcForwardBackwardSums<SPOp>(forwardProbs, backwardProbs, forwardBackwardSums);
		
		// compute all arc probs at once
		probs.setAll(0.0);
		addForwardBackwardSums<SPOp>(forwardBackwardSums, parentModel_->getParentSetFactors(), probs);

		delete[] forwardBackwardSums;
	}

	void getArcProbs(ArcMap<T>& probs) const {
		getUnnormArcProbs(probs);
		probs /= getUnnormProb();
	}
};


#endif
