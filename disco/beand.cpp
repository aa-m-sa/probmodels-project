/*
 *  BEANDisco: main program
 *  
 *  Copyright 2011-2015 Teppo Niinimäki <teppo.niinimaki(at)helsinki.fi>
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

#include <algorithm>
#include <iostream>
#include <iomanip> 
#include <fstream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <exception>

#include <cstdarg>

#include <sys/resource.h> 

#include "debug.hpp"

#include "common.hpp"
#include "format.hpp"
#include "logger.hpp"
#include "lognum.hpp"
#include "timer.hpp"
#include "data.hpp"
#include "scores.hpp"
#include "adtree.hpp"
#include "parentsetmap.hpp"
#include "bucketorder.hpp"
//#include "parbucketorder.hpp"
//#include "bodist.hpp"

#include "dist/dist.hpp"
#include "dist/mcmcdist.hpp"
#include "dist/mc3dist.hpp"
#include "dist/aisdist.hpp"
#include "dist/conddist.hpp"

#include "io.hpp"
#include "sampleio.hpp"

#include "estimator.hpp"
#include "model.hpp"
#include "pomodel.hpp"
#include "omodel.hpp"

#include "options.hpp"

#include "interrupt.hpp"

//#include "subsetmap.hpp"
//#include "sortedarraysubset.hpp"

//#define NDEBUG
#include <cassert>


#define BEAND_VERSION_STRING  "2.0"

#define BEAND_COMPILED_DATE __DATE__ " " __TIME__


// type definitions
typedef Lognum<double> Real;
//typedef double Real;
//#include <gmpxx.h>
//typedef Lognum<mpf_class> Real;
//template <typename T>
//T to(mpf_class x){
//	return T(x.get_d());
//}

#include "sampler.hpp"


// create a logger
Logger logger;

#include "exact.hpp"

Options options;

std::unique_ptr<OutStream> logStream(nullptr);


struct InOutStreams {
	// input streams
	std::shared_ptr<InStream> sampleInStream;
	std::shared_ptr<InStream> subSampleInStream;
	
	// sample output streams
	std::shared_ptr<OutStream> sampleOutStream;
	std::shared_ptr<OutStream> subSampleCountOutStream;
	std::shared_ptr<OutStream> levelWeightOutStream;
	std::shared_ptr<OutStream> levelInvWeightOutStream;
	std::shared_ptr<OutStream> subSampleOutStream;
	std::shared_ptr<OutStream> sampleProbOutStream;
	std::shared_ptr<OutStream> sampleNormConstRatioOutStream;
	std::shared_ptr<OutStream> sampleWeightOutStream;

	// estimate output streams
	std::shared_ptr<OutStream> normConstOutStream;
	std::shared_ptr<OutStream> normConstHarmonicOutStream;
	//std::shared_ptr<OutStream> featureOutStream;
	std::shared_ptr<OutStream> arcProbsOutStream;

	// accept ratios
	//std::shared_ptr<OutStream> aisAcceptRatioStream;
	//std::shared_ptr<OutStream> mc3AcceptRatioStream;
	//std::shared_ptr<OutStream> mcmcAcceptRatioStream;

	std::shared_ptr<OutStream> samplingTimeOutStream;

	// open margin stream if needed
	//unique_ptr<OutStream> marginStream(openOutStream(sampleProbOutFilename, true));
	//marginStream->stream().setf(std::ios::scientific);
	//marginStream->stream().precision(10);

	std::function<void(const ArcMap<Real>& probs)> writeArcProbs;
	std::function<void(const Real& normConst)> writeNormConst;
	std::function<void(const Real& normConst)> writeNormConstHarmonic;
	std::function<void(const double& samplingTime)> writeSamplingTime;
	std::function<void(const Real& samplingTime)> writeAcceptRatios;

	void openStreams() {
		sampleInStream = openInStream(options.sampleInFilename);

		sampleOutStream = openOutStream(options.sampleOutFilename);
		subSampleOutStream = openOutStream(options.subSampleOutFilename);
		subSampleCountOutStream = openOutStream(options.subSampleCountOutFilename);
		levelWeightOutStream = openOutStream(options.levelWeightOutFilename);
		levelInvWeightOutStream = openOutStream(options.levelInvWeightOutFilename);
		sampleProbOutStream = openOutStream(options.sampleProbOutFilename);
		sampleNormConstRatioOutStream = openOutStream(options.sampleNormConstRatioOutFilename);
		sampleWeightOutStream = openOutStream(options.sampleWeightOutFilename);


		//aisAcceptRatioStream = openOutStream(options.aisAcceptRatioOutFilename,
		//		options.aisAcceptRatioOutType == "update");
		//mc3AcceptRatioStream = openOutStream(options.mc3AcceptRatioOutFilename,
		//		options.mc3AcceptRatioOutType == "update");
		//mcmcAcceptRatioStream = openOutStream(options.mcmcAcceptRatioOutFilename,
		//		options.mcmcAcceptRatioOutType == "update");
		normConstOutStream = openOutStream(options.normConstOutFilename,
				options.normConstOutType == "update");
		normConstHarmonicOutStream = openOutStream(options.normConstHarmonicOutFilename,
				options.normConstHarmonicOutType == "update");
		arcProbsOutStream = openOutStream(options.arcProbsOutFilename,
				options.arcProbsOutType == "update");
		samplingTimeOutStream = openOutStream(options.samplingTimeOutFilename,
				options.samplingTimeOutType == "update");
		//samplingTimeOutStream->stream().setf(std::ios::fixed);
		//samplingTimeOutStream->stream().precision(3);
		//std::resetioflags(std::ios::floatfield)

		if (options.arcProbsOutFormat == "pretty") {
			writeArcProbs = [&] (const ArcMap<Real>& probs) {
				arcProbsOutStream->newField();
				for (auto arc: allArcs(probs.getSize())) {
					arcProbsOutStream->stream()
							<< arc << " "
							<< std::setiosflags(std::ios::fixed)
							<< std::setprecision(10)
							<< to<double>(probs[arc]) << "\n";
				}
				arcProbsOutStream->endField();
			};
		}
		else if (options.arcProbsOutFormat == "plain-row") {
			writeArcProbs = [&] (const ArcMap<Real>& probs) {
				arcProbsOutStream->newField();
				int first = true;
				for (auto arc: allArcs(probs.getSize())) {
					if (!first)
						arcProbsOutStream->stream() << " ";
					arcProbsOutStream->stream()
							<< std::setiosflags(std::ios::fixed)
							<< std::setprecision(10)
							<< to<double>(probs[arc]);
					first = false;
				}
				arcProbsOutStream->stream() << "\n";
				arcProbsOutStream->endField();
			};
		}
		else {
			throw Exception("Invalid arc probability output format '%s'.", options.arcProbsOutFormat);
		}

		writeNormConst = [&] (const Real& normConst) {
			normConstOutStream->newField();
			normConstOutStream->stream()
					<< std::setiosflags(std::ios::scientific)
					<< std::setprecision(10)
					<< log(normConst) << "\n";
			normConstOutStream->endField();
		};

		writeNormConstHarmonic = [&] (const Real& normConst) {
			normConstHarmonicOutStream->newField();
			normConstHarmonicOutStream->stream()
					<< std::setiosflags(std::ios::scientific)
					<< std::setprecision(10)
					<< log(normConst) << "\n";
			normConstHarmonicOutStream->endField();
		};

		writeSamplingTime = [&] (const double& samplingTime) {
			samplingTimeOutStream->newField();
			samplingTimeOutStream->stream()
					<< std::setiosflags(std::ios::fixed)
					<< std::setprecision(3)
					<< samplingTime << "\n";
			samplingTimeOutStream->endField();
		};
	}
};

InOutStreams inOutStreams;

#include "sampling.hpp"



std::string toString(const boost::any& value) {
	std::stringstream ss;
	if (auto v = boost::any_cast<int>(&value))
		ss << *v;
	else if (auto v = boost::any_cast<unsigned int>(&value))
		ss << *v;
	else if (auto v = boost::any_cast<long long int>(&value))
		ss << *v;
	else if (auto v = boost::any_cast<unsigned long long int>(&value))
		ss << *v;
	else if (auto v = boost::any_cast<size_t>(&value))
		ss << *v;
	else if (auto v = boost::any_cast<double>(&value))
		ss << *v;
	else if (auto v = boost::any_cast<bool>(&value))
		ss << *v;
	else if (auto v = boost::any_cast<std::string>(&value))
		ss << *v;
	else if (auto v = boost::any_cast<Flag>(&value))
		ss << *v;
	else if (auto v = boost::any_cast<AutoFlag>(&value))
		ss << *v;
	else if (auto v = boost::any_cast<AutoInt>(&value))
		ss << *v;
	else if (auto v = boost::any_cast<std::vector<std::string>>(&value))
		ss << *v;
	else if (auto v = boost::any_cast<std::list<std::string>>(&value))
		ss << *v;
	else if (auto v = boost::any_cast<std::vector<Arc>>(&value))
		ss << *v;
	else if (auto v = boost::any_cast<IntRangeList>(&value))
		ss << *v;
	else if (auto v = boost::any_cast<NoneOrInt>(&value))
		ss << *v;
	else
		assert(0);
	return ss.str();
}

/*
 * Main program.
 */

#include <string>
#include <iomanip>
//#include <boost/program_options.hpp>

using namespace std;
//namespace opts = boost::program_options;


int main(int argc, char** argv) {
	
	try {
		options.parse(argc, argv);
	} catch (Exception& e) {
		logger.printfln(-1, "Error: %s", e.what());
		logger.printfln(-1, "Type '%s --help' to see help message.", argv[0]);
		logger.println(-1, "Aborting.");
		return 1;
	}
	
	auto printNameAndVersion = [] {
		logger.println(-1, "BEANDisco - Bayesian Exact and Approximate Network Discovery");
		logger.printfln(-1, "Version %s (build date: %s)", BEAND_VERSION_STRING, BEAND_COMPILED_DATE);
	};

	if (options.showVersion) {
		printNameAndVersion();
		return 0;
	}

	if (options.showHelp) {
		printNameAndVersion();
		logger.println(-1);
		logger.println(-1, "Usage:");
		logger.printfln(-1, "  %s [options] [infile [outfile]]", argv[0]);
		logger.println(-1);
		logger.println(-1, options.allOpts);
		return 0;
	}

	if (options.quiet)
		logger.setVerbosity(-1);
	else
		logger.setVerbosity(options.verbosity);

	// open log stream for statistics
	logStream = openOutStream(options.logOutFilename);
	logStream->stream().setf(ios::fixed);
	logStream->stream().precision(2);
	logStream->stream() << "version = " << BEAND_VERSION_STRING << endl;
	logStream->stream() << "build-date = " << BEAND_COMPILED_DATE << endl;

	logger.println(3, "Parameters:");
	for (const auto& option : options.vm) {
		const auto& name = option.first;
		auto value = toString(option.second.value());
		logger.printfln(3, "  %s = %s", name, value);
		logStream->stream() << name << " = " << value << endl;
	}
	
	// initialize rng
	//rng.seed(rngSeed);
	rng.seed(options.rngSeed);
	
	
	// start global timer
	Timer timer;
	timer.start();
	
	// local scores
	shared_ptr<ParentSetMap<Real>> scores(nullptr);
	
	// can't read both data and scores
	if (!options.dataInFilename.empty() && !options.scoreInFilename.empty()) {
		logger.println(-1, "Error: only data or only score file should be provided.");
		logger.printfln(-1, "Type '%s --help' to see help message.", argv[0]);
		return 1;
	}
	// read the data and compute the scores
	else if (!options.dataInFilename.empty()) {
		// read data
		logger.println(1, "Reading data...");
		shared_ptr<Data> data(new Data());
		try {
			unique_ptr<InStream> dataInStream(openInStream(options.dataInFilename));
			data->read(dataInStream->stream());
		} catch (Exception& e) {
			logger.printfln(-1, "Error: While reading data file '%s': %s",
					options.dataInFilename, e.what());
			return 1;
		}
		
		// extract selected rows and columns from the data
		if (!options.dataVariables.empty() || !options.dataSamples.empty()) {
			if (options.dataVariables.empty())
				options.dataVariables.addRange(0, data->getNumVariables() - 1);
			if (options.dataSamples.empty())
				options.dataSamples.addRange(0, data->getNumSamples() - 1);
			shared_ptr<Data> filteredData(new Data());
			filteredData->selectFrom(*data,
					options.dataVariables.getIntList(0, data->getNumVariables() - 1),
					options.dataSamples.getIntList(0, data->getNumSamples() - 1));
			data = filteredData;
		}
		
		logger.printfln(2, "  Number of variables = %d", data->nVariables);
		logger.printfln(2, "  Number of samples = %d", data->nSamples);
		
		int maxIndegree = options.maxIndegree.isSpecial() ?
				data->getNumVariables() - 1 : 
				options.maxIndegree.value();

		if (!(0 <= maxIndegree && maxIndegree < data->getNumVariables())) {
			logger.println(-1, "Error: The maximum in-degree should be between 0 and n-1.");
			return 1;
		}
		
		// construct dataView (plain data or ADTree)
		shared_ptr<DataView> dataView = data;
		if (options.useADTree) {
			logger.println(1, "Constructing ADTree...");
			Timer adTreeTimer; adTreeTimer.start();
			shared_ptr<ADTree> adTree(new ADTree(*data, options.adTreeMinCount,
					options.adTreeMaxDepth, maxIndegree + 1));
			double adTreeTime = adTreeTimer.elapsed();
			logStream->stream() << "adtree_construction_time = " << adTreeTime << endl;
			logger.printfln(2, "  Number of ADNodes = %d", adTree->getNumADNodes());
			logger.printfln(1, "  Elapsed %s.", prettyDuration(adTreeTime));
			dataView = adTree;
		}
		
		// compute the maximum number of configurations for a variable set and for a single variable
		int maxNumValuesSet = 1;
		vector<int> arities(data->arities, data->arities + data->nVariables);
		sort(arities.begin(), arities.end());
		for (int i = 1; i <= maxIndegree + 1; ++i)
			maxNumValuesSet *= arities[data->nVariables - i];
		int maxNumValuesSingle = arities[data->nVariables - 1];
		
		
		// compute scores
		logger.println(1, "Computing scores...");
		Timer scoreTimer; scoreTimer.start();
		logger.println(2, "  Preparing for score computation...");
		ScoreFun* scoreFun = nullptr;
		if (options.scoreType == "k2") {
			if (options.precomputeGammas)
				scoreFun = new K2ScoreCached(data->nSamples + maxNumValuesSingle);
			else
				scoreFun = new K2Score();
		} else if (options.scoreType == "bdeu") {
			if (options.precomputeGammas)
				scoreFun = new BDeuScoreCached(options.equivalentSampleSize, data->nSamples, maxNumValuesSet);
			else
				scoreFun = new BDeuScore(options.equivalentSampleSize);
		} else if (options.scoreType == "ll") {
			scoreFun = new LLScore();
		} else if (options.scoreType == "mdl") {
			scoreFun = new MDLScore(data->nSamples);
		} else if (options.scoreType == "aic") {
			scoreFun = new AICScore();
		} else {
			logger.printfln(-1, "Error: Invalid score type '%s'", options.scoreType);
			return 1;
		}
		logger.println(2, "  Computing the actual scores...");
		//model.buildScoresFromData(dataView.get(), scoreFun, maxIndegree);
		scores = move(buildScoresFromData<Real>(dataView.get(), scoreFun, maxIndegree));
		delete scoreFun;
		double scoreTime = scoreTimer.elapsed();
		logStream->stream() << "score_computation_time = " << scoreTime << endl;
		logger.printfln(1, "  Elapsed %s.", prettyDuration(scoreTime));
		
	}
	// read the scores
	else if (!options.scoreInFilename.empty()) {
		logger.println(1, "Reading scores...");
		try {
			unique_ptr<InStream> scoreInStream(openInStream(options.scoreInFilename));
			//model.readScores(scoreInStream->stream());
			scores = move(readScores<Real>(scoreInStream->stream()));
		} catch (Exception& e) {
			logger.printfln(-1, "Error: While reading score file '%s': %s",
					options.scoreInFilename, e.what());
			return 1;
		}
	}
	// complain if neither data nor scores was given
	else {
		logger.println(-1, "Error: either data or score file should be provided.");
		logger.printfln(-1, "Type '%s --help' to see help message.", argv[0]);
		return 1;
	}
	
	// optionally write scores to file
	if (!options.scoreOutFilename.empty()) {	
		logger.println(1, "Writing scores...");
		unique_ptr<OutStream> scoreOutStream(openOutStream(options.scoreOutFilename));
		//model.writeScores(scoreOutStream->stream());
		writeScores(scoreOutStream->stream(), scores.get());
	}

	DagFamily dagFamily(scores->nNodes);
	OrderFamily orderFamily(dagFamily.n);

	// apply structure prior factor to scores
	if (options.structurePriorFunction == "uniform") {
		// nothing to be done
	}
	else if (options.structurePriorFunction == "indegree-uniform") {
		applyIndegreeUniformPrior(scores.get());
	}
	else {
		logger.println(-1, "Error: invalid structure prior factor.");
		logger.printfln(-1, "Type '%s --help' to see help message.", argv[0]);
		return 1;
	}

	// create order prior factor
	PredSetMap<Real> orderPrior;
	if (options.orderPriorFunction == "uniform") {
		// nothing to be done
	}
	else {
		logger.println(-1, "Error: invalid order prior factor.");
		logger.printfln(-1, "Type '%s --help' to see help message.", argv[0]);
		return 1;
	}


	// build the target model
	shared_ptr<StructureModel<Real>> model;
	if (options.priorType == "modular") {
		model = make_shared<ModularModel<Real>>(dagFamily, scores.get());
	}
	else if (options.priorType == "ordermodular") {
		model = make_shared<OrderModularModel<Real>>(dagFamily, orderFamily, scores.get(), &orderPrior);
	}
	else {
		logger.println(-1, "Error: invalid prior type.");
		logger.printfln(-1, "Type '%s --help' to see help message.", argv[0]);
		return 1;
	}

	
	// open result stream for writing
	//bool anytimeEstimate = (options.estimateOutType != "final");
	//bool overwriteEstimate = (options.estimateOutType == "update");
	//unique_ptr<OutStream> estimateOutStream(openOutStream(options.estimateOutFilename, overwriteEstimate));


	
	// actual computation

	if (options.computationMethod == "exact") {
		logger.println(1, "Computing exactly...");
		Timer taskTimer; taskTimer.start();
		inOutStreams.openStreams();
		if (contains(options.taskType, "normconst")) {
			logger.println(1, "  Computing the exact normalizing constant...");
			Real normConst;
			if (auto orderModularModel = dynamic_cast<const OrderModularModel<Real>*>(model.get())) {
				logger.println(2, "    Using the dynamic programming algorithm (by Koivisto & Sood).");
				normConst = computeNormConstDP(orderModularModel);
			}
			else if (auto modularModel = dynamic_cast<const ModularModel<Real>*>(model.get())) {
				// TODO
				logger.printfln(-1, "Error: exact normconst for modular model not implemented.");
				return 1;
				//logger.println(2, "    Using the inclusion-exclusion algorithm (by Ellis & Wong).");
				//normConst = computeNormConstIE(modularModel);
			}
			else {
				assert(0);
			}
			inOutStreams.writeNormConst(normConst);
		}
		if (contains(options.taskType, "arcprobs")) {
			logger.println(1, "  Computing the exact arc probabilities...");
			ArcMap<Real> arcProbs(model->getNumVariables());
			if (auto orderModularModel = dynamic_cast<const OrderModularModel<Real>*>(model.get())) {
				logger.println(2, "    Using the dynamic programming algorithm (by Koivisto & Sood).");
				arcProbs = computeArcProbsDP(orderModularModel);
			}
			else if (auto modularModel = dynamic_cast<const ModularModel<Real>*>(model.get())) {
				// TODO
				logger.printfln(-1, "Error: exact arcprobs for modular model not implemented.");
				return 1;
				//logger.println(2, "    Using the inclusion-exclusion algorithm (by Ellis & Wong).");
				//normConst = computeArcProbsIE(modularModel);
			}
			else {
				assert(0);
			}
			inOutStreams.writeArcProbs(arcProbs);
		}
		/*if (contains(options.taskType, "bestsample")) {
			logger.println(1, "  Finding the best sample...");
			// TODO
			logger.printfln(-1, "Error: finding best sample not implemented.");
			return 1;
			//Real maxUnnormProb = calcMaxPermExact(*scores);
			//resStream << log(maxUnnormProb) << std::endl;
		}*/
		if (contains(options.taskType, "bestdag")) {
			logger.println(1, "  Finding the best structure...");
			// TODO
			logger.printfln(-1, "Error: finding best dag not implemented.");
			return 1;
		}
		double elapsedTaskTime = taskTimer.elapsed();	
		logger.printfln(1, "  Elapsed %s.", prettyDuration(elapsedTaskTime));
	}
	else if (options.computationMethod == "ais" ||
			options.computationMethod == "mc3" ||
			options.computationMethod == "mcmc" ||
			options.computationMethod == "mc") {
		
		inOutStreams.openStreams();

		logger.println(1, "Running the sampler...");
		try {
			selectAndRunSampler(model.get());
		} catch (Exception& e) {
			logger.printfln(-1, "Error: %s", e.what());
			return 1;
		}
	}
	else if (options.computationMethod != "") {
		logger.printfln(-1, "Error: invalid computation method '%s'.", options.computationMethod);
		return 1;
	}

	// actual computation ends
	
	// print timing
	double elapsedTime = timer.elapsed();
	logger.printfln(0, "Elapsed %s total.", prettyDuration(elapsedTime));
	logStream->stream() << "time = " << elapsedTime << endl;
	
	struct rusage rusage;
	getrusage(RUSAGE_SELF, &rusage);

	size_t maxRSS = rusage.ru_maxrss * 1024;
	logger.printfln(1, "Used max %g of RAM.", prettySize<units::B>(maxRSS));
	logStream->stream() << "max-memory-rss = " << maxRSS << endl;

	double totalUserTime = rusage.ru_utime.tv_sec + rusage.ru_utime.tv_usec / 1e6;
	logger.printfln(2, "Used %g of user CPU time.", prettyDuration(totalUserTime));
	logStream->stream() << "total-user-time = " << totalUserTime << endl;

	double totalSystemTime = rusage.ru_stime.tv_sec + rusage.ru_stime.tv_usec / 1e6;
	logger.printfln(2, "Used %g of system CPU time.", prettyDuration(totalSystemTime));
	logStream->stream() << "total-system-time = " << totalSystemTime << endl;

	return 0;
}



