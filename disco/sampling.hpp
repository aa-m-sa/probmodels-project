/*
 *  BEANDisco: sampling
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

#include <cassert>

#ifndef SAMPLING_HPP
#define SAMPLING_HPP

std::vector<double> createAnnealingPowers(int nLevels, double levelParam, const std::string& levelScheme) {
	std::vector<double> annealingPowers;
	if (levelParam > nLevels + 1)
		levelParam = nLevels + 1;
	else if (levelParam < -(nLevels + 1))
		levelParam = -(nLevels + 1);
	
	for (int i = 0; i < nLevels; ++i) {
		if (levelScheme == "linear")
			annealingPowers.push_back((double)(i+1) / (double)(nLevels + 1));
		else if (levelScheme == "autogeometric")
			annealingPowers.push_back(exp((i - nLevels) / sqrt(nLevels)));
		else if (levelScheme == "geometric") {
			if (levelParam < 0) {
				double power = (1.0 - pow(1.0 + levelParam/(nLevels+1), i+1)) / 
						(1.0 - pow(1.0 + levelParam/(nLevels+1), nLevels+1));
				annealingPowers.push_back(power);
			} else if (levelParam == 0) {
				annealingPowers.push_back((double)(i+1) / (double)(nLevels + 1));
			} else {
				double power = 1.0 - (1.0 - pow(1.0 - levelParam/(nLevels+1), nLevels-i)) / 
						(1.0 - pow(1.0 - levelParam/(nLevels+1), nLevels+1));
				annealingPowers.push_back(power);
			}
		}
		else {
			throw Exception("Invalid level scheme '%s'.", levelScheme);
		}
	}
	//for (int i = 0; i < aisTargetDistSteps; ++i)
	//	annealingPowers.push_back(1.0);
	
	return annealingPowers;
}


template <typename SampleSpace>
std::unique_ptr<MCProposalDist<SampleSpace>>
createMCProposalDist(const SampleSpace& sampleSpace, const std::string& name);


template <>
std::unique_ptr<MCProposalDist<BucketOrderFamily>>
createMCProposalDist<BucketOrderFamily>(const BucketOrderFamily& sampleSpace, const std::string& name) {
	MCProposalDist<BucketOrderFamily>* mcProposalDist;
	if (name == "swap") {
		mcProposalDist = new BucketOrderFamily::UniformRandSwap(sampleSpace);
	} else {
		throw Exception("Invalid proposal distribution '%s' for bucket orders.", name);
	}
	return std::unique_ptr<MCProposalDist<BucketOrderFamily>>(mcProposalDist);
}

template <>
std::unique_ptr<MCProposalDist<OrderFamily>>
createMCProposalDist<OrderFamily>(const OrderFamily& sampleSpace, const std::string& name) {
	MCProposalDist<OrderFamily>* mcProposalDist;
	if (name == "swap") {
		mcProposalDist = new OrderFamily::UniformRandSwap(sampleSpace);
	} else {
		throw Exception("Invalid proposal distribution '%s' for orders.", name);
	}
	return std::unique_ptr<MCProposalDist<OrderFamily>>(mcProposalDist);
}



template <typename SampleSpace, typename Model>
std::unique_ptr<ModelDist<SampleSpace, Model, Real>> createDist(const SampleSpace& sampleSpace, const Model* model) {
//ModelDist<SampleSpace, Model, Real>* createDist(const SampleSpace& sampleSpace, const Model* model) {
	ModelDist<SampleSpace, Model, Real>* sampleDist = nullptr;
	if (options.computationMethod == "ais") {
		std::vector<double> annealingPowers; 
		/*int levels = options.aisNLevels;
		if (options.aisNLevels.timesDataSize) {
			levels *= model->
		}*/
		annealingPowers = createAnnealingPowers(options.aisNLevels,
				options.aisLevelParam, options.aisLevelScheme);
		std::unique_ptr<MCProposalDist<SampleSpace>> proposalDist(
				createMCProposalDist<SampleSpace>(sampleSpace, options.aisProposalDist));
		sampleDist = new AnnealedDist<SampleSpace, Model, Real>(sampleSpace, model,
				annealingPowers, std::move(proposalDist), options.aisTargetDistSamples,
				options.aisSampleSpacing);
	}
	else if (options.computationMethod == "mc3") {
		std::vector<double> annealingPowers;
		annealingPowers = createAnnealingPowers(options.mc3NLevels, 
				options.mc3LevelParam, options.mc3LevelScheme);
		//annealingPowers.insert(annealingPowers.begin(), 0.0);
		std::unique_ptr<MCProposalDist<SampleSpace>> proposalDist(
				createMCProposalDist<SampleSpace>(sampleSpace, options.mc3ProposalDist));
		auto mc3Dist = new MC3Dist<SampleSpace, Model, Real>(sampleSpace, model,
				annealingPowers, std::move(proposalDist), options.mc3SampleSpacing,
				options.mc3LevelSwapsPerStep);
		if (options.mc3BurninSteps > 0) {
			logger.println(2, "  Running burn-in...");
			Timer burninTimer;
			burninTimer.start();
			mc3Dist->burn(options.mc3BurninSteps);
			double samplerBurninTime = burninTimer.elapsed();	
			logger.printfln(1, "    Elapsed %s.", prettyDuration(samplerBurninTime));
			logStream->stream() << "sampler_burnin_time = " << samplerBurninTime << std::endl;
		}
		sampleDist = mc3Dist;
	}
	else if (options.computationMethod == "mcmc") {
		std::unique_ptr<MCProposalDist<SampleSpace>> proposalDist(
				createMCProposalDist<SampleSpace>(sampleSpace, options.mcmcProposalDist));
		auto mcmcDist = new MCMCDist<SampleSpace, Model, Real>(sampleSpace, model,
				std::move(proposalDist), options.mcmcSampleSpacing);
		if (options.mcmcBurninSteps > 0) {
			logger.println(2, "  Running burn-in...");
			Timer burninTimer;
			burninTimer.start();
			mcmcDist->burn(options.mcmcBurninSteps);
			double samplerBurninTime = burninTimer.elapsed();	
			logger.printfln(1, "    Elapsed %s.", prettyDuration(samplerBurninTime));
			logStream->stream() << "sampler_burnin_time = " << samplerBurninTime << std::endl;
		}
		sampleDist = mcmcDist;
	}
	else {
		assert(0);
		throw Exception("Invalid computation method '%s'.", options.computationMethod);
	}
	return std::unique_ptr<ModelDist<SampleSpace, Model, Real>>(sampleDist);
	//return sampleDist;
}






template <typename SampleSpace, typename Model>
std::unique_ptr<SampleSource<SampleSpace, Model, Real>> getSampleSource(
		const SampleSpace& sampleSpace, const Model* model) {
	SampleSource<SampleSpace, Model, Real>* sampleSource = nullptr;
	if (options.computationMethod == "mc") {
		logger.println(2, "  Attaching the sample input file...");
		if (!inOutStreams.sampleInStream) {
			throw Exception("Sample input file not given.");
		}
		auto sampleFileSource = new SampleFileSource<SampleSpace, Model, Real>(sampleSpace, model, inOutStreams.sampleInStream.get());
		sampleFileSource->skip(options.mcBurninSteps);
		sampleSource = std::move(sampleFileSource);
	}
	else {
		logger.println(2, "  Creating the sample generator...");
		auto boDist = createDist<SampleSpace, Model>(sampleSpace, model);
		sampleSource = new SampleDistSource<SampleSpace, Model, Real>(std::move(boDist));
	}
	return std::unique_ptr<SampleSource<SampleSpace, Model, Real>>(sampleSource);
}


template <typename Model>
struct Estimators {
	std::shared_ptr<Estimator<Real>> unweightedNormConstEstimator;
	std::shared_ptr<Estimator<Real>> unweightedHarmonicNormConstEstimator;
	std::shared_ptr<Estimator<Real>> normConstRatioEstimator;
	//std::shared_ptr<Estimator<Real>> normConstEstimator;
	std::shared_ptr<Estimator<Real>> unArcProbEstimator;
	std::shared_ptr<Estimator<ArcMap<Real>>> unArcProbsEstimator;
};

template <typename SampleSpace, typename Model>
std::unique_ptr<Estimators<RestrictedModel<Model, SampleSpace>>>
addSampleEstimators(
		const SampleSpace& sampleSpace,
		Sampler<SampleSpace, Model>* sampler) {
	std::unique_ptr<Estimators<RestrictedModel<Model, SampleSpace>>> estimators(
		new Estimators<RestrictedModel<Model, SampleSpace>>());

	auto normConstRatioEstimator =
			std::make_shared<ExactNormConstRatioEstimator<SampleSpace, Model, Real>>();
	sampler->addSink(normConstRatioEstimator);
	estimators->normConstRatioEstimator = normConstRatioEstimator;

	if (contains(options.taskType, "arcprobs")) {
		if (options.arcs.size() == 1) {
			ArcProbFeature<Real> feature(options.arcs[0]);
			auto unArcProbEstimator = std::make_shared<
					ExactFeatureProbEstimator<SampleSpace, Model, ArcProbFeature<Real>, Real>>(feature);
			sampler->addSink(unArcProbEstimator);
			estimators->unArcProbEstimator = unArcProbEstimator;
		}
		else {
			auto unArcProbsEstimator = std::make_shared<
					ExactArcProbsEstimator<SampleSpace, Model, Real>>(sampleSpace);
			sampler->addSink(unArcProbsEstimator);
			estimators->unArcProbsEstimator = unArcProbsEstimator;
		}
	}

	return estimators;
}

template <typename SampleSpace, typename Model>
std::unique_ptr<Estimators<RestrictedModel<Model, SampleSpace>>>
addWeightedEstimators(
		const SampleSpace& sampleSpace,
		Sampler<SampleSpace, Model>* sampler,
		const Estimators<RestrictedModel<Model, SampleSpace>>* sourceEstimators,
		std::shared_ptr<Estimator<Real>> weightEstimator) {
	std::unique_ptr<Estimators<RestrictedModel<Model, SampleSpace>>> estimators(
		new Estimators<RestrictedModel<Model, SampleSpace>>());

	auto normConstRatioEstimator =
			std::make_shared<WeightedEstimator<SampleSpace, Model, Real, Real>>(
					sourceEstimators->normConstRatioEstimator, weightEstimator, Real(0));
	sampler->addSink(normConstRatioEstimator);
	estimators->normConstRatioEstimator = normConstRatioEstimator;

	if (sourceEstimators->unArcProbEstimator) {
		auto unArcProbEstimator =
			std::make_shared<WeightedEstimator<SampleSpace, Model, Real, Real>>(
					sourceEstimators->unArcProbEstimator, weightEstimator, Real(0));
		sampler->addSink(unArcProbEstimator);
		estimators->unArcProbEstimator = unArcProbEstimator;
	}

	if (sourceEstimators->unArcProbsEstimator) {
		auto unArcProbsEstimator =
			std::make_shared<WeightedEstimator<SampleSpace, Model, ArcMap<Real>, Real>>(
					sourceEstimators->unArcProbsEstimator, weightEstimator, ArcMap<Real>(sampleSpace.n, Real(0)));
		sampler->addSink(unArcProbsEstimator);
		estimators->unArcProbsEstimator = unArcProbsEstimator;
	}

	return estimators;
}

template <typename SampleSpace, typename Model>
std::unique_ptr<Estimators<Model>>
addISEstimators(
		const SampleSpace& sampleSpace,
		Sampler<SampleSpace, Model>* sampler,
		//std::shared_ptr<const MeanType<Real>> meanType,
		double meanFraction,
		const Estimators<RestrictedModel<Model, SampleSpace>>* sampleEstimators) {
	std::unique_ptr<Estimators<Model>> estimators(new Estimators<Model>());

	auto normConstRatioEstimator = 
			std::make_shared<ISNormConstRatioEstimator<SampleSpace, Model, Real>>(meanFraction, 
				sampleEstimators->normConstRatioEstimator);
	sampler->addSink(normConstRatioEstimator);
	estimators->normConstRatioEstimator = normConstRatioEstimator;

	if (contains(options.taskType, "normconst")) {
		auto unweightedNormConstEstimator = std::make_shared<
				ISMultiLevelNormConstEstimator<SampleSpace, Model, Real>>(sampler->nLevels(), meanFraction);
		sampler->addSink(unweightedNormConstEstimator);
		estimators->unweightedNormConstEstimator = unweightedNormConstEstimator;
	}

	if (contains(options.taskType, "normconst-harmonic")) {
		auto unweightedHarmonicNormConstEstimator = std::make_shared<
				ISMultiLevelHarmonicNormConstEstimator<SampleSpace, Model, Real>>(sampler->nHarmonicLevels(), meanFraction);
		sampler->addSink(unweightedHarmonicNormConstEstimator);
		estimators->unweightedHarmonicNormConstEstimator = unweightedHarmonicNormConstEstimator;
	}

	if (sampleEstimators->unArcProbEstimator) {
		auto unArcProbEstimator = std::make_shared<
				ISUnnormFeatureProbEstimator<SampleSpace, Model, Real, Real>>(meanFraction, 
					sampleEstimators->unArcProbEstimator, Real(0));
		sampler->addSink(unArcProbEstimator);
		estimators->unArcProbEstimator = unArcProbEstimator;
	}

	if (sampleEstimators->unArcProbsEstimator) {
		auto unArcProbsEstimator = std::make_shared<
				ISUnnormFeatureProbEstimator<SampleSpace, Model, ArcMap<Real>, Real>>(meanFraction, 
					sampleEstimators->unArcProbsEstimator, ArcMap<Real>(sampleSpace.n, Real(0)));
		sampler->addSink(unArcProbsEstimator);
		estimators->unArcProbsEstimator = unArcProbsEstimator;
	}

	return estimators;
}


template <typename SampleSpace, typename Model>
void runSampler(const SampleSpace& sampleSpace, const Model* model, bool weightDags = false) {
		//std::shared_ptr<SampleEstimator<DagFamily, RestrictedModel<Model, SampleSpace>, Real>> weightEstimator = nullptr) {
	using SampleModel = RestrictedModel<Model, SampleSpace>;

	Timer timer;
	timer.start();

	// get/build sample source
	std::shared_ptr<SampleSource<SampleSpace, Model, Real>> sampleSource =
			getSampleSource(sampleSpace, model);

	// build the sampler
	logger.println(2, "  Building the sampler...");
	auto sampler = std::make_shared<Sampler<SampleSpace, Model>>(sampleSource->nLevels(), sampleSource->nHarmonicLevels());

	// add possible sample file sink
	if (inOutStreams.sampleOutStream) {
		logger.println(2, "  Attaching the sample output file...");
		auto sampleSink = std::make_shared<SampleFileSink<SampleSpace, Model, Real>>(
				inOutStreams.sampleOutStream.get());
		sampler->addSink(sampleSink);
	}

	// add possible subsample count sink
	if (inOutStreams.subSampleCountOutStream) {
		sampler->addSink(std::make_shared<SubSampleCountFileSink<SampleSpace, Model, Real>>(
				inOutStreams.subSampleCountOutStream.get()));
	}

	// add possible sample weight sink
	if (inOutStreams.sampleWeightOutStream) {
		sampler->addSink(makeFileSink<SampleSpace, Model, Real>(
			inOutStreams.sampleWeightOutStream.get(),
			[] (OutStream* out, const Sample<SampleSpace, Model, Real>* sample) {
				out->stream() << std::setiosflags(std::ios::scientific);
				out->stream() << std::setprecision(10);
				out->stream() << log(sample->getWeight());
				out->stream() << std::endl;
			}));
	}

	// add possible sample probability sink
	if (inOutStreams.sampleProbOutStream) {
		sampler->addSink(makeFileSink<SampleSpace, Model, Real>(
			inOutStreams.sampleProbOutStream.get(),
			[] (OutStream* out, const Sample<SampleSpace, Model, Real>* sample) {
				out->stream() << std::setiosflags(std::ios::scientific);
				out->stream() << std::setprecision(10);
				out->stream() << log(sample->getUnnormProb());
				out->stream() << std::endl;
			}));
	}

	// add possible multilevel weight sink
	if (inOutStreams.levelWeightOutStream) {
		sampler->addSink(makeFileSink<SampleSpace, Model, Real>(
			inOutStreams.levelWeightOutStream.get(),
			[] (OutStream* out, const Sample<SampleSpace, Model, Real>* sample) {
				out->stream() << std::setiosflags(std::ios::scientific);
				out->stream() << std::setprecision(10);
				for (auto w : sample->getNormConstWeights())
					out->stream() << log(w) << " ";
				out->stream() << std::endl;
			},
			options.levelWeightOutSpacing));
	}

	// add possible multilevel inverse weight sink
	if (inOutStreams.levelInvWeightOutStream) {
		sampler->addSink(makeFileSink<SampleSpace, Model, Real>(
			inOutStreams.levelInvWeightOutStream.get(),
			[] (OutStream* out, const Sample<SampleSpace, Model, Real>* sample) {
				out->stream() << std::setiosflags(std::ios::scientific);
				out->stream() << std::setprecision(10);
				for (auto w : sample->getInvNormConstWeights())
					out->stream() << log(w) << " ";
				out->stream() << std::endl;
			},
			options.levelInvWeightOutSpacing));
	}

	std::unique_ptr<Estimators<SampleModel>> sampleEstimators;

	if (options.subSampleDags) {
		using SubSampleSpace = DagFamily;
		using SubDist = ExactCondDist<SubSampleSpace, Model, SampleSpace, Real>;
		const SubSampleSpace& subSampleSpace = model->getDagFamily();
		auto subSourceSource = std::make_shared<SubSampleDistSource<SubSampleSpace, SampleModel, SubDist, Real>>(
				subSampleSpace);
		auto nestedSampler = std::make_shared<Sampler<SubSampleSpace, SampleModel>>(0, 0);

		// add possible subsample file sink
		if (inOutStreams.subSampleOutStream) {
			logger.println(2, "  Attaching the subsample output file...");
			auto subSampleSink = std::make_shared<SampleFileSink<SubSampleSpace, SampleModel, Real>>(
					inOutStreams.subSampleOutStream.get());
			nestedSampler->addSink(subSampleSink);
		}

		// create and add estimators
		auto subSampleEstimators = addSampleEstimators(subSampleSpace, nestedSampler.get());
		if (weightDags) {
			auto weightEstimator = std::make_shared<ExactLinExtCountInvEstimator<DagFamily, SampleModel, Real>>();
			nestedSampler->addSink(weightEstimator);
			subSampleEstimators = addWeightedEstimators(subSampleSpace, nestedSampler.get(), std::move(subSampleEstimators).get(), weightEstimator);
		}

		//sampleEstimators = addISEstimators(nestedSampler.get(), std::make_shared<FullMean<Real>>(), subSampleEstimators.get());
		sampleEstimators = addISEstimators(subSampleSpace, nestedSampler.get(), 1.0, subSampleEstimators.get());

		// create the subsampler instance
		if (options.subSamplingType == "weighted") {
			auto subSampler = std::make_shared<SubSampler<SampleSpace, Model, SubSampleSpace, SampleModel>>(
					subSourceSource, nestedSampler, options.nSubSamples, options.subSamplingTimeTarget,
					options.subSamplingTimingFraction);
			sampler->setSubSampler(subSampler);
		}
		else if (options.subSamplingType == "ellis-wong") {
			size_t maxTotalSubSamples = (size_t)-1;
			if (options.maxTotalEWSubSamples > 0) {
				maxTotalSubSamples = options.maxTotalEWSubSamples;
			}
			else if (options.maxTotalEWSubSampleMem > 0) {
				maxTotalSubSamples = options.maxTotalEWSubSampleMem /
						((2 + sampleSpace.n) * sizeof(void*));			// size of EWDag
			}

			auto subSampler = std::make_shared<EWSubSampler<SampleSpace, Model, SubSampleSpace, SampleModel>>(
					subSourceSource, nestedSampler, options.subSamplingTotalProbFractionEps, maxTotalSubSamples);
			sampler->setSubSampler(subSampler);
		}
		else {
			throw Exception("Invalid subsampling type '%s'", options.subSamplingType);
		}
	}
	else {
		sampleEstimators = addSampleEstimators(sampleSpace, sampler.get());
	}

	// add possible sample probability sink
	if (inOutStreams.sampleNormConstRatioOutStream) {
		sampler->addSink(makeFileSink<SampleSpace, Model, Real>(
			inOutStreams.sampleNormConstRatioOutStream.get(),
			[&] (OutStream* out, const Sample<SampleSpace, Model, Real>* sample) {
				out->stream() << std::setiosflags(std::ios::scientific);
				out->stream() << std::setprecision(10);
				out->stream() << log(sampleEstimators->normConstRatioEstimator->getEstimate());
				out->stream() << std::endl;
			}));
	}

	//logger.println(2, "  Creating estimators...");
	//std::shared_ptr<MeanType<Real>> meanType = std::make_shared<FullMean<Real>>();
	//auto estimators = addISEstimators(sampler.get(), meanType, sampleEstimators.get());
	double meanFraction = 1.0 - options.burninFraction;
	auto estimators = addISEstimators(sampleSpace, sampler.get(), meanFraction, sampleEstimators.get());

	std::vector<std::function<void()>> perSampleEstimateCallbacks;
	std::vector<std::function<void()>> finalEstimateCallbacks;

	auto addOutputCallback = [&] (std::function<void()> func, const std::string& type, const std::string& optionName) {
		if (type == "final") {
			finalEstimateCallbacks.push_back(func);
		}
		else if (type == "update" || type == "cumulate") {
			perSampleEstimateCallbacks.push_back(func);
		}
		else {
			throw Exception("Invalid %s output type '%s'", optionName, type);
		}
	};

	std::shared_ptr<Estimator<Real>> normConstEstimator(nullptr);
	if (contains(options.taskType, "normconst")) {
		normConstEstimator = std::make_shared<NormConstEstimator<Real>>(estimators->unweightedNormConstEstimator, estimators->normConstRatioEstimator);
		auto outputFunc = [&] {
				inOutStreams.writeNormConst(normConstEstimator->getEstimate());
			};
		addOutputCallback(outputFunc, options.normConstOutType, "normconst");
	}

	std::shared_ptr<Estimator<Real>> harmonicNormConstEstimator(nullptr);
	if (contains(options.taskType, "normconst-harmonic")) {
		harmonicNormConstEstimator = std::make_shared<NormConstEstimator<Real>>(estimators->unweightedHarmonicNormConstEstimator, estimators->normConstRatioEstimator);
		auto outputFunc = [&] {
				inOutStreams.writeNormConstHarmonic(harmonicNormConstEstimator->getEstimate());
			};
		addOutputCallback(outputFunc, options.normConstHarmonicOutType, "normconst-harmonic");
	}

	std::shared_ptr<Estimator<Real>> arcProbEstimator(nullptr);
	std::shared_ptr<Estimator<ArcMap<Real>>> arcProbsEstimator(nullptr);
	if (contains(options.taskType, "arcprobs")) {
		if (options.arcs.size() == 1) {
			arcProbEstimator = std::make_shared<FeatureProbEstimator<Real, Real>>(estimators->unArcProbEstimator,
					estimators->normConstRatioEstimator);
			auto outputFunc = [&] {
					inOutStreams.arcProbsOutStream->newField();
					inOutStreams.arcProbsOutStream->stream()
							<< options.arcs[0] << " "
							<< to<double>(arcProbEstimator->getEstimate()) << "\n";
					inOutStreams.arcProbsOutStream->endField();
				};
			addOutputCallback(outputFunc, options.arcProbsOutType, "arcprobs");
		}
		else {
			arcProbsEstimator = std::make_shared<FeatureProbEstimator<ArcMap<Real>, Real>>(estimators->unArcProbsEstimator,
					estimators->normConstRatioEstimator);
			auto outputFunc = [&] {
					auto estimate = arcProbsEstimator->getEstimate();
					inOutStreams.writeArcProbs(estimate);
				};
			addOutputCallback(outputFunc, options.arcProbsOutType, "arcprobs");
		}
	}

	Timer samplingTimer;

	// add sampling time output
	if (contains(options.taskType, "sampling-time")) {
		auto outputFunc = [&] {
				inOutStreams.writeSamplingTime(samplingTimer.elapsed());
			};
		addOutputCallback(outputFunc, options.samplingTimeOutType, "sampling time");
	}

	if (contains(options.taskType, "accept-ratio")) {
		// TODO!!
		throw Exception("Accept ratio output not implemented.");
		/*auto outputFunc = [&] {
				inOutStreams.writeAcceptRatio();
				samplingTimeOutStream->newField();
				for (int i = 0; i < acceptRatios.size(); ++i) {
					samplingTimeOutStream->stream()
							<< std::setiosflags(std::ios::fixed)
							<< std::setprecision(3)
							<< acceptRatios << "\n";
				}
				samplingTimeOutStream->endField();
			};
		addOutputCallback(outputFunc, options.normConstOutType, "accept ratio");*/
	}

	// set interrupt handler
	bool interrupt = false;
	auto interruptHandler = [&] (int sig) {
		interrupt = true;
		if (options.interruptMidSample)
			logger.printfln(-1, "\r    INTERRUPT: Stopping as soon as possible and then quit."
					" Press CTRL+C again to force quit.");
		else
			logger.printfln(-1, "\r    INTERRUPT: Waiting for current sample to finish and then quit."
					" Press CTRL+C again to force quit.");
		resetInterruptHandler();
	};
	resetInterruptHandler(interruptHandler);

	// define per sample callback
	double lastStatusUpdate = 0;
	sampler->setCallback(
		[&] (int sampleNum) -> bool {
			if ((sampleNum + 1) % options.estimateOutputSpacing == 0) {
				for (auto& callback : perSampleEstimateCallbacks) {
					callback();
				}
			}
			double elapsedTime = samplingTimer.elapsed();
			// output elapsed time
			//inOutStreams.samplingTimeOutStream->writelnField(elapsedTime);
			// print elapsed and etimated remaining time
			double remainingTime = 0;
			if (options.samplingTimeTarget)
				remainingTime = options.samplingTimeTarget - elapsedTime;
			if (options.nSamples)
				remainingTime = (options.nSamples - sampleNum - 1) * elapsedTime / (sampleNum + 1);
			if (elapsedTime >= lastStatusUpdate + 0.1) {
				logger.printf(2, "    Samples drawn: %d%s",
						sampleNum + 1, options.nSamples ? format("/%d", options.nSamples).str() : "");
				logger.printf(2, ", Elapsed time: %s", prettyDuration(elapsedTime));
				if (remainingTime)
					logger.printf(2, ", Estimated time remaining: %s", prettyDuration(remainingTime));
				logger.printf(2, "\r");
				lastStatusUpdate = elapsedTime;
			}
			return !interrupt;
		});
	
	if (options.interruptMidSample) {
		sampleSource->setInterruptCallback(
			[&] {
				return interrupt || (options.samplingTimeTarget &&
						samplingTimer.elapsed() > options.samplingTimeTarget);
			});
	}

	// run the sampler
	logger.println(2, "  Sampling...");
	samplingTimer.start();
	auto nSampled = sampler->run(sampleSource.get(), options.nSamples, options.samplingTimeTarget);
	double samplerSamplingTime = samplingTimer.elapsed();

	if (nSampled < options.nSamples && sampleSource->eof())
		logger.printfln(2, "    Ran out of samples.");
	logger.printfln(2, "    Total samples drawn: %d%s, Elapsed %s.", nSampled,
			options.nSamples ? format("/%d", options.nSamples).str() : "",
			prettyDuration(samplerSamplingTime));
	logStream->stream() << "sampler-samples-drawn = " << nSampled << std::endl;
	logStream->stream() << "sampler-sampling-time = " << samplerSamplingTime << std::endl;

	// reset default interrupt handler
	resetInterruptHandler();

	// output final estimates
	logger.println(2, "  Outputting the final estimates...");
	for (auto& callback : finalEstimateCallbacks) {
		callback();
	}

	// elapsed total time
	double samplerTime = timer.elapsed();	
	logger.printfln(1, "  Elapsed %s.", prettyDuration(samplerTime));
	logStream->stream() << "sampler-total-time = " << samplerTime << std::endl;
}

template <typename ModelType>
int selectBucketSize(ModelType model) {
	if (options.maxBucketSize.isSpecial()) {
		unsigned long totNumParentSets = 0; 
		int n = model->getNumVariables();
		for (int i = 0; i < n; ++i)
			totNumParentSets += model->getParentSetFactors()->forNode(i)->numSubsets();
		int b = 1;
		while (((unsigned long)1 << b) * (unsigned long)b * (unsigned long)n <= totNumParentSets)
			++b;
		if (b > 1)
			--b;
		logger.printfln(1, "  Select bucket size %d", b);
		logStream->stream() << "automatic-bucket-size = " << b << std::endl;
		return b;
	}
	else {
		return options.maxBucketSize.value();
	}
}

//template <typename Model>
//std::unique_ptr<AbstractSampler> createSampler(Model* model) {
//void selectAndRunSampler(const Model* model) {
void selectAndRunSampler(const StructureModel<Real>* model) {
	if (options.sampleType == "bucketorder") {
		//std::shared_ptr<BucketOrderFamily> bof(new BucketOrderFamily(model->getNumVariables(), options.maxBucketSize));
		//int maxBucketSize = selectBucketSize(model);
		//BucketOrderFamily bof(model->getNumVariables(), maxBucketSize);
		if (options.subSampleDags) {
			if (auto orderModularModel = dynamic_cast<const OrderModularModel<Real>*>(model)) {
				int maxBucketSize = selectBucketSize(orderModularModel);
				BucketOrderFamily bof(model->getNumVariables(), maxBucketSize);
				runSampler(bof, orderModularModel);
			}
			else if (auto modularModel = dynamic_cast<const ModularModel<Real>*>(model)) {
				PredSetMap<Real> orderPrior;
				OrderFamily orderFamily(modularModel->getNumVariables());
				auto orderModularModel = std::make_shared<OrderModularModel<Real>>(modularModel->getDagFamily(),
						orderFamily, modularModel->getParentSetFactors(), &orderPrior);
				int maxBucketSize = selectBucketSize(orderModularModel);
				BucketOrderFamily bof(model->getNumVariables(), maxBucketSize);
				//auto weightEstimator = std::make_shared<ExactLinExtCountInvEstimator<DagFamily, RestrictedModel<OrderModularModel, BucketOrderFamily>, Real>>();
				bool weightDags = (options.subSamplingType == "weighted");
				runSampler(bof, orderModularModel.get(), weightDags);
				/*std::unique_ptr<Estimator<DagFamily::Instance, Real>> weightEstimator(
						new ExactLinExtCountEstimator<Real>());*/
				//auto subSampleSource = new SampleBiaser<DagFamily, RestrictedModel<OrderModularModel<Real>, BucketOrderFamily>, ModularModel<Real>, Real>(
				//		subSampleSource, weightEstimator.get(), modularModel);
				//auto targetWeighting = make_shared<OrderModularVsModularRatio<Real>>(modularModel);
				//runSampler(bof, sampligModel, targetWeighting.get());
			}
			else {
				assert(0);
			}
		}
		else {
			if (auto orderModularModel = dynamic_cast<const OrderModularModel<Real>*>(model)) {
				int maxBucketSize = selectBucketSize(orderModularModel);
				BucketOrderFamily bof(model->getNumVariables(), maxBucketSize);
				runSampler(bof, orderModularModel);
			}
			else if (dynamic_cast<const ModularModel<Real>*>(model)) {
				throw Exception("Basic partial order sampling not possible with modular prior.");
			}
			else {
				assert(0);
			}
		}
	}
	else if (options.sampleType == "order") {
		if (options.subSampleDags) {
			if (auto orderModularModel = dynamic_cast<const OrderModularModel<Real>*>(model)) {
				runSampler(orderModularModel->getOrderFamily(), orderModularModel);
			}
			else if (auto modularModel = dynamic_cast<const ModularModel<Real>*>(model)) {
				PredSetMap<Real> orderPrior;
				OrderFamily orderFamily(modularModel->getNumVariables());
				auto orderModularModel = std::make_shared<OrderModularModel<Real>>(modularModel->getDagFamily(),
						orderFamily, modularModel->getParentSetFactors(), &orderPrior);
				bool weightDags = (options.subSamplingType == "weighted");
				runSampler(orderModularModel->getOrderFamily(), orderModularModel.get(), weightDags);
			}
			else {
				assert(0);
			}
		}
		else {
			if (auto orderModularModel = dynamic_cast<const OrderModularModel<Real>*>(model)) {
				runSampler(orderModularModel->getOrderFamily(), orderModularModel);
			}
			else if (dynamic_cast<const ModularModel<Real>*>(model)) {
				throw Exception("Basic order sampling not possible with modular prior.");
			}
			else {
				assert(0);
			}
		}
	}
	else if (options.sampleType == "dag") {
		// TODO
		throw Exception("Sample type 'dag' not implemented.");
		//DagFamily df(model.getNumVariables());
		//runSampler(df, model);
	}
	else {
		throw Exception("Invalid sample type '%s'.", options.sampleType);
		//sampler = nullptr;
	}
	//return std::unique_ptr<AbstractSampler>(sampler);
}



#endif

