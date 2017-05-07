/*
 *  BEANDisco: exact full model computations
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


#include "model.hpp"

#ifndef EXACT_HPP
#define EXACT_HPP


template <typename T>
T computeNormConstDP(const OrderModularModel<T>* model) {
	using SPOp = SumProdOp<T>;
	int n = model->getNumVariables();
	auto scores = model->getParentSetFactors();

	logger.println(2, "    Computing predecessor set probabilities (alphas)...");
	std::vector<FullSubsetMap<T>*> predSetProbs(n);
	SortedArraySubset allNodes = SortedArraySubset::fullSet(n);
	SortedArraySubset potParents(n, allNodes);
	for (int node = 0; node < n; ++node) {
		potParents.remove(node);
		predSetProbs[node] = new FullSubsetMap<T>(potParents);
		predSetProbs[node]->setAll(0);
		for (auto subsetScore : scores->forNode(node)->forSubsets(
				typename FullSubsetMap<T>::Subset(*predSetProbs[node]))) {
			const auto& parentSet = subsetScore.first;
			const auto& score = subsetScore.second;
			(*predSetProbs[node])[parentSet] = score;
		}
		predSetProbs[node]->template zetaTransform<SPOp>();
		/*for (predSetProb : predSetProbs[node]->forSubsets()) {
			const auto& predSet = predSetProb.first;
			auto& score = score.second;
			score *= (*model->getPredSetFactors()->forNode(node))[predSet];
		}*/
		potParents.insert(node);
	}

	logger.println(2, "    Computing forward sums...");
	FullSubsetMap<T> forwardProbs(allNodes);
	forwardProbs.template forwardSum<SPOp>(predSetProbs);
	
	return forwardProbs.getFull();
}


/*emplate <typename T>
Real computeNormConstTianHe(const ModularModel<T>* model) {
	int n = model->getNumVariables();
	auto scores = model->getParentSetFactors();
	std::vector<FullSubsetMap<Real>*> predSetProbs(n);
	SortedArraySubset allNodes = SortedArraySubset::fullSet(n);
	SortedArraySubset potParents(allNodes);
	for (int node = 0; node < n; ++node) {
		potParents.remove(node);
		predSetProbs[node] = new FullSubsetMap<Real>(potParents);
		predSetProbs[node]->setAll(0);
		for (auto subsetScore : scores->forNode(node)->forSubsets(
				FullSubsetMap<Real>::Subset(*predSetProbs[node]))) {
			const auto& parentSet = subsetScore.first;
			const auto& score = subsetScore.second;
			(*predSetProbs[node])[parentSet] = score;
		}
		predSetProbs[node]->zetaTransform<SumProdOp<Real>>();
		potParents.insert(node);
	}


	restrictedProbs[0] = 0;

	I nonrootMask = 0;
	int nonrootSize = 0;
	while (nonrootMask < (I(1) << n) - I(1)) {
		{
			int elem = 0;
			while (true) {
				I eMask = I(1) << elem;
				nonrootMask ^= eMask
				sign = -sign;
				if (subrootMask & eMask)
					break;
				++elem;
				--nonrootSize;
			}
			nonrootElements[nonrootSize] = elem;
			++nonrootSize;
		}

		I subrootMask = 0;
		T sign = -1;
		T val = 0;
		for (I subrootMaskComp = 1; subrootMaskComp < (I(1) << nonrootSize); ++subrootMaskComp) {
			{
				int i = 0;
				while (true) {
					I eMask = I(1) << nonrootElements[i];
					subrootMask ^= eMask
					sign = -sign;
					if (subrootMask & eMask) {
						subrootProbs[subrootMask] = subrootProbs[subrootMask ^ eMask] *
								(*predSetProbs[nonrootElements[i]])[]
						break;
					}
					++i;
				}
			}
			assert((subrootMask & ~nonrootMask) == 0);
			val += sign * restrictedProbs[nonrootMask ^ subrootMask];
		}
		restrictedProbs[nonrootMask] = val;
	}

	return restrictedProbs[(I(1) << n) - I(1)];
}/**/



template <typename T>
ArcMap<T> computeArcProbsDP(const OrderModularModel<T>* model) {
	int n = model->getNumVariables();
	auto scores = model->getParentSetFactors();

	logger.println(2, "    Computing predecessor set probabilities (alphas)...");
	std::vector<FullSubsetMap<T>*> predSetProbs(n);
	SortedArraySubset allNodes = SortedArraySubset::fullSet(n);
	SortedArraySubset potParents(n, allNodes);
	for (int node = 0; node < n; ++node) {
		potParents.remove(node);
		predSetProbs[node] = new FullSubsetMap<T>(potParents);
		predSetProbs[node]->setAll(0);
		for (auto subsetScore : scores->forNode(node)->forSubsets(
				typename FullSubsetMap<T>::Subset(*predSetProbs[node]))) {
			const auto& parentSet = subsetScore.first;
			const auto& score = subsetScore.second;
			(*predSetProbs[node])[parentSet] = score;
		}
		predSetProbs[node]->template zetaTransform<SumProdOp<T>>();
		potParents.insert(node);
	}

	logger.println(2, "    Computing forward sums...");
	FullSubsetMap<T> forwardProbs(allNodes);
	forwardProbs.template forwardSum<SumProdOp<T>>(predSetProbs);

	logger.println(2, "    Computing backward sums...");
	FullSubsetMap<T> backwardProbs(allNodes);
	backwardProbs.template backwardSum<SumProdOp<T>>(predSetProbs);

	logger.println(2, "    Computing forward-backward sums...");
	std::vector<FullSubsetMap<T>*> forwardBackwardSums(n);
	for (int node = 0; node < n; ++node) {
		potParents.remove(node);
		forwardBackwardSums[node] = new FullSubsetMap<T>(potParents);
		for (auto subsetSum : forwardBackwardSums[node]->forSubsets(
				typename FullSubsetMap<T>::Subset(forwardProbs))) {
			const auto& downSet = subsetSum.first;
			auto& fpsum = subsetSum.second;
			T fp = forwardProbs[downSet];
			auto downSetPlusNode = downSet;
			downSetPlusNode.insert(node);
			T bp = backwardProbs[downSetPlusNode];
			fpsum = SumProdOp<T>::prod(fp, bp);
		}
		forwardBackwardSums[node]->template upZetaTransform<SumProdOp<T>>();
		potParents.insert(node);
	}
	
	// compute all arc probs at once
	logger.println(2, "    Computing arc probs...");
	ArcMap<T> probs(n);
	probs.setAll(0.0);
	for (int node = 0; node < n; ++node) {
		for (auto subsetScore : scores->forNode(node)->forSubsets(
				SetPair<SortedArraySubset, typename FullSubsetMap<T>::Subset>(SortedArraySubset(n),
				typename FullSubsetMap<T>::Subset(*forwardBackwardSums[node])))) {
			const auto& parentSetElems = subsetScore.first.first;
			const auto& parentSet = subsetScore.first.second;
			const auto& score = subsetScore.second;
			T prod = SumProdOp<T>::prod(score, (*forwardBackwardSums[node])[parentSet]);
			for (auto parent : parentSetElems) {
				SumProdOp<T>::add(probs[Arc(parent, node)], prod);
			}
			
		}
	}
	//addForwardBackwardSums<SPOp>(forwardBackwardSums, parentModel_->getParentSetFactors(), probs);

	return probs / forwardProbs.getFull();
}/**/


#endif

