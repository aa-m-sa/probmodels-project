/*
 *  BEANDisco: command line option handling
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


#include <string>
#include <list>
#include <vector>
#include <fstream>
#include <boost/program_options.hpp>
//#include <regex>
#include <boost/regex.hpp>

#include "optionparser.hpp"

#ifndef OPTIONS_HPP
#define OPTIONS_HPP


/**
 * Flag
 */
enum Flag {NO = false, YES = true};

std::ostream& operator<<(std::ostream& os, const Flag& flag) {
	switch (flag) {
		case NO:
			os << "no";
			break;
		case YES:
			os << "yes";
			break;
		default:
			assert(0);
	}
	return os;
}


template <>
Flag validate<Flag>(const std::string& s) {
	if (s == "yes" || s == "true" || s == "1") {
		return YES;
	}
	else if (s == "no" || s == "false" || s == "0") {
		return NO;
	}
	else
		throw boost::program_options::invalid_option_value(s);
}/**/



/**
 * SpecialOr
 */
template <const char* Special, typename T>
struct SpecialOr {
	bool isSpecial_;
	T value_;
public:
	/*SpecialOr(const std::string& value) {
		if (value != Special)
			throw Exception("Invalid")
		isSpecial_ = true;
	}*/
	SpecialOr() {
		isSpecial_ = true;
	}
	SpecialOr(const T& value) {
		isSpecial_ = false;
		value_ = value;
	}
	bool isSpecial() const {
		return isSpecial_;
	}
	T value() const {
		return value_;
	}
};

template <const char* Special, typename T>
std::ostream& operator<<(std::ostream& os, const SpecialOr<Special, T>& v) {
	if (v.isSpecial())
		os << Special;
	else
		os << v.value();
	return os;
}

template <const char* S, typename T>
SpecialOr<S, T> validateSpecialOr(const std::string& s) {
	if (s == S)
		return SpecialOr<S, T>();
	else
		return SpecialOr<S, T>(validate<T>(s));
}


/**
 * AutoFlag
 */
extern const char Auto[] = "auto";
using AutoFlag = SpecialOr<Auto, Flag>;

/**
 * AutoInt
 */
using AutoInt = SpecialOr<Auto, int>;


/**
 * IntRangeList
 */
struct IntRangeList {
private:
	struct Range {
		int first;
		int last;
		Range(int f, int l) : first(f), last(l) {}
	};
	
	std::list<Range> ranges_;
	friend std::ostream& operator<<(std::ostream& os, const IntRangeList& irList);
public:
	void addRange(int first, int last) {
		ranges_.push_back(Range(first, last));
	}
	std::list<int> getIntList(int min, int max) {
		std::list<int> ints;
		for (const auto& range : ranges_)
			for (int i = std::max(min, range.first); i <= std::min(max, range.last); ++i)
				ints.push_back(i);
		return ints;
	}
	bool empty() const {
		return ranges_.empty();
	}
};

std::ostream& operator<<(std::ostream& os, const IntRangeList& irList) {
	int firstGone = false;
	for (const auto& range : irList.ranges_) {
		if (firstGone)
			os << ",";
		os << range.first << "-" << range.last;
		firstGone = true;
	}
	return os;
}


std::ostream& operator<<(std::ostream& os, const std::list<std::string>& strList) {
	int firstGone = false;
	for (const auto& str : strList) {
		if (firstGone)
			os << ",";
		os << str;
		firstGone = true;
	}
	return os;
}



/*struct AisLevels {
	bool timesDataSize;
	int value;
};

std::ostream& operator<<(std::ostream& os, const AisLevels& v) {
	os << v.value;
	if (v.timesDataSize)
		os << 'D';
	return os;
}

template <>
AisLevels validate<AisLevels>(const std::string& s) {
	static boost::regex rNumLevels(
			"^\\s*(\\d+(?:\\.\\d+)?)\\s?(D)?$");
	boost::smatch match;
	if (boost::regex_match(s, match, rNumLevels)) {
		AisLevels aisLevels;
		aisLevels.value = std::stod(match[1]);
		aisLevels.timesDataSize = false;
		if (match[2].matched) {
			aisLevels.timesDataSize = true;
		}
		return aisLevels;
	}
	else {
		throw boost::program_options::invalid_option_value(s);
	}
}*/

/**
 * Duration
 */
double validateDuration(const std::string& s) {
	static boost::regex rDuration(
			"^"
			"(\\s*(\\d+(?:\\.\\d+)?)\\s?(?:d|days))?"
			"(\\s*(\\d+(?:\\.\\d+)?)\\s?(?:h|hours))?"
			"(\\s*(\\d+(?:\\.\\d+)?)\\s?(?:m|min|minutes))?"
			"(\\s*(\\d+(?:\\.\\d+)?)\\s?(?:s|sec|seconds)?)?"
			"$");
	boost::smatch match;
	if (boost::regex_match(s, match, rDuration)) {
		double days = match[1].matched ? std::stod(match[2]) : 0;
		double hous = match[3].matched ? std::stod(match[4]) : 0;
		double mins = match[5].matched ? std::stod(match[6]) : 0;
		double secs = match[7].matched ? std::stod(match[8]) : 0;
		return ((days * 24 + hous) * 60 + mins) * 60 + secs;
	}
	else {
		throw boost::program_options::invalid_option_value(s);
	}
}

/**
 *
 */
size_t validateSize(const std::string& s) {
	static boost::regex rSize(
			"^\\s*(\\d+(?:\\.\\d+)?)\\s?(\\w+)?$");
	boost::smatch match;
	if (boost::regex_match(s, match, rSize)) {
		double size = std::stod(match[1]);
		if (match[2].matched) {
			std::string unit = match[2];
			if (boost::regex_match(unit, boost::regex("k"))) {
				size *= 1024;
			}
			else if (boost::regex_match(unit, boost::regex("M"))) {
				size *= std::pow(1024, 2);
			}
			else if (boost::regex_match(unit, boost::regex("G"))) {
				size *= std::pow(1024, 3);
			}
			else if (boost::regex_match(unit, boost::regex("T"))) {
				size *= std::pow(1024, 4);
			}
			else if (boost::regex_match(unit, boost::regex("P"))) {
				size *= std::pow(1024, 5);
			}
			else if (boost::regex_match(unit, boost::regex("E"))) {
				size *= std::pow(1024, 6);
			}
			else {
				throw boost::program_options::invalid_option_value(s);
			}
		}
		return (size_t)size;
	}
	else {
		throw boost::program_options::invalid_option_value(s);
	}
}


template <typename Set>
std::string validateOneOf(const std::string& s, const Set& set) {
	if (!contains(set, s))
		throw boost::program_options::invalid_option_value(s);
	return s;
}

template <typename Set>
std::vector<std::string> validateSubsetOf(const std::string& s, const Set& set) {
	std::vector<std::string> ss;
	for (auto elem : split(s, ",")) {
		if (!contains(set, elem) || contains(ss, elem))
			throw boost::program_options::invalid_option_value(s);
		ss.push_back(elem);
	}
	return ss;
}


namespace opts = boost::program_options;

void validate(boost::any& v, const std::vector<std::string>& values,
		Flag* target_type, int) {
	opts::validators::check_first_occurrence(v);
	const std::string& s = opts::validators::get_single_string(values);

	v = boost::any(validate<Flag>(s));
}/**/


template <const char* S, typename T>
void validate(boost::any& v, const std::vector<std::string>& values,
		SpecialOr<S, T>* target_type, int) {
	opts::validators::check_first_occurrence(v);
	const std::string& s = opts::validators::get_single_string(values);

	v = boost::any(validateSpecialOr<S, T>(s));
}/**/


void validate(boost::any& v, const std::vector<std::string>& values,
		std::list<std::string>* target_type, int) {
	opts::validators::check_first_occurrence(v);
	const std::string& s = opts::validators::get_single_string(values);

	std::list<std::string> strs;
	for (auto str : split(s, ",")) {
		strs.push_back(str);
	}
	v = boost::any(strs);
}/**/

void validate(boost::any& v, const std::vector<std::string>& values,
		std::vector<Arc>* target_type, int) {
	opts::validators::check_first_occurrence(v);
	const std::string& s = opts::validators::get_single_string(values);

	static boost::regex r("^(0|[1-9][0-9]*)->?(0|[1-9][0-9]*)$");
	//static std::regex r("^(0|[1-9][0-9]*)->?(0|[1-9][0-9]*)$");
	std::vector<Arc> arcs;
	for (auto arcStr : split(s, ",")) {
		boost::smatch match;
		if (!boost::regex_match(arcStr, match, r)) {
			throw opts::invalid_option_value(s);
			//throw opts::invalid_option_value(
			//		format("invalid arc list item '%s'", arcStr).str());
			//throw Exception("invalid arc list item '%s'", arcStr);
		}
		arcs.push_back(Arc(std::stoi(match[1]), std::stoi(match[2])));
	}
	v = boost::any(arcs);
}/**/

template <>
IntRangeList validate<IntRangeList>(const std::string& s) {
	static boost::regex rValue("^0|[1-9][0-9]*$");
	static boost::regex rRange("^(0|[1-9][0-9]*)?-(0|[1-9][0-9]*)?$");
	
	IntRangeList intRanges;
	for (auto rangeStr : split(s, ",")) {
		boost::smatch match;
		if (boost::regex_match(rangeStr, match, rValue)) {
			int value = std::stoi(match[0]);
			intRanges.addRange(value, value);
		}
		else if (boost::regex_match(rangeStr, match, rRange)) {
			int first = match[1].matched ?
					std::stoi(match[1]) :
					std::numeric_limits<int>::min();
			int last = match[2].matched ?
					std::stoi(match[2]) :
					std::numeric_limits<int>::max();
			intRanges.addRange(first, last);
		}
		else {
			throw opts::invalid_option_value(s);
			//throw opts::invalid_option_value(
			//		format("invalid number range list item '%s'", rangeStr).str());
		}
	}
	return intRanges;
}

void validate(boost::any& v, const std::vector<std::string>& values,
		IntRangeList* target_type, int) {
	opts::validators::check_first_occurrence(v);
	const std::string& s = opts::validators::get_single_string(values);

	v = boost::any(validate<IntRangeList>(s));
}/**/


extern const char None[] = "none";
using NoneOrInt = SpecialOr<None, int>;


class Options {
//private:
//	Logger& logger;
public:
	//using namespace std;
	using string = std::string;
	
	//opts::options_description desc;
	opts::options_description generalOpts;
	opts::options_description modelOpts;
	opts::options_description methodOpts;
	opts::options_description inputOpts;
	opts::options_description outputOpts;
	opts::options_description allOpts;
	opts::variables_map vm;
	

	// === general options ===
	bool showHelp;
	bool showVersion;
	int verbosity;
	bool quiet;
	string configFilename;
	unsigned int rngSeed;
	

	// === input options ===
	string dataInFilename;
	IntRangeList dataVariables;
	IntRangeList dataSamples;
	//string dataVariablesStr;
	//string dataSamplesStr;
	string scoreInFilename;
	string sampleInFilename;


	// === model options ===
	string priorType;
	string structurePriorFunction;
	string orderPriorFunction;
	//int maxIndegree;
	NoneOrInt maxIndegree;
	string scoreType;
	double equivalentSampleSize;


	// === method options ===
	//std::string taskTypeStr;
	std::vector<std::string> taskType;
	std::vector<Arc> arcs;
	string computationMethod;

	// score computation
	Flag precomputeGammas;
	Flag useADTree;
	int adTreeMinCount;
	int adTreeMaxDepth;
	
	// sampling in general
	double burninFraction;

	// mc
	int mcBurninSteps;

	// mcmc
	int mcmcSampleSpacing;
	int mcmcBurninSteps;
	string mcmcProposalDist;
	
	// mc3
	int mc3SampleSpacing;
	int mc3LevelSwapsPerStep;
	int mc3BurninSteps;
	int mc3NLevels;
	string mc3LevelScheme;
	double mc3LevelParam;
	string mc3ProposalDist;
	
	// ais
	int aisNLevels;
	string aisLevelScheme;
	double aisLevelParam;
	int aisTargetDistSamples;
	int aisSampleSpacing;
	string aisProposalDist;
	
	// per sample
	string sampleType;
	//int nBuckets;
	AutoInt maxBucketSize;
	int nSamples;
	double samplingTimeTarget; //string samplingTimeTargetStr;
	Flag interruptMidSample;
	// subsamples
	bool subSampleDags;
	AutoFlag subSampleDagsExplicit;
	int nSubSamples;
	double subSamplingTimeTarget; //string subSamplingTimeTargetStr;
	double subSamplingTimingFraction;
	string subSamplingType;
	double subSamplingTotalProbFractionEps;
	size_t maxTotalEWSubSamples; //string maxTotalEWSubSamplesStr;
	size_t maxTotalEWSubSampleMem; //string maxTotalEWSubSampleMemStr;


	// === output options ===

	// scores
	string scoreOutFilename;

	// estimates
	int estimateOutputSpacing;
	string normConstOutFilename;
	string normConstOutType;
	string normConstHarmonicOutFilename;
	string normConstHarmonicOutType;
	string arcProbsOutFilename;
	string arcProbsOutType;
	string arcProbsOutFormat;
	string featureOutFilename;
	string featureOutType;

	// per sample info
	string sampleOutFilename;
	string subSampleOutFilename;
	string levelWeightOutFilename;
	int levelWeightOutSpacing;
	string levelInvWeightOutFilename;
	int levelInvWeightOutSpacing;
	string sampleProbOutFilename;
	string sampleNormConstRatioOutFilename;
	string sampleWeightOutFilename;
	string subSampleCountOutFilename;

	// statistics
	string aisAcceptRatioOutFilename;
	string aisAcceptRatioOutType;
	string mc3AcceptRatioOutFilename;
	string mc3AcceptRatioOutType;
	string mcmcAcceptRatioOutFilename;
	string mcmcAcceptRatioOutType;
	string samplingTimeOutFilename;
	string samplingTimeOutType;
	string logOutFilename;
	
	
	Options() :
		//logger(_logger)
		generalOpts("General options"),
		modelOpts("Model options"),
		methodOpts("Algorithmic options"),
		inputOpts("Input option"),
		outputOpts("Output option"),
		allOpts("All options")
	{
		generalOpts.add_options()
			("help,h",                 opts::bool_switch(&showHelp),
			                           "print a help message and exit")
			("version",                opts::bool_switch(&showVersion),
			                           "print the version and exit")
			("verbose,v",              value<int>(&verbosity)
			                           ->default_value(0)->implicit_value(1),
				                       "set the verbosity level")
			("quiet,q",                opts::bool_switch(&quiet),
			                           "use a quiet mode (do not print anything unnecessary)")
			("config-file",            value<string>(&configFilename)
			                           ->default_value(""),
				                       "specify a configuration file to read options from")
			("seed",                   value<unsigned int>(&rngSeed)
				                       ->default_value_text("random")
			                           ->validator([&] (const string& s) -> unsigned int { 
			                           		if (s == "random") return std::random_device()();
			                           		else return validate<unsigned int>(s);
			                           }),
				                       "set a seed for the random number generator"
				                       ", possible values: random, <non-negative integer>")
			;

		inputOpts.add_options()
			("data-in-file,D",         value<string>(&dataInFilename),
				                       "specify a data input file")
			("data-vars",              value<IntRangeList>(&dataVariables),
				                       "select a subset of data variables"
				                       ", example values: '2,5,7' or '1,3,5-9' or '2-8'")
			("data-rows",              value<IntRangeList>(&dataSamples),
				                       "select a subset of data rows"
				                       ", example values: '1-100' or '10-20,40-80' or '1,3,5,7,9-99'")
	//		("data-variables",         value<string>(&dataVariablesStr)->default_value("-"),
	//		                           "select the data variables")
	//		("data-rows",              value<string>(&dataSamplesStr)->default_value("-"),
	//		                           "select the data rows")
			("score-in-file",          value<string>(&scoreInFilename),
				                       "specify a score input file")
			("sample-in-file",         value<string>(&sampleInFilename),
				                       "specify a sample input file")
			;

		modelOpts.add_options()
			("prior-type",             value<string>(&priorType)
			                           ->default_value("ordermodular"),
				                       "set structure prior type"
				                       ", possible values: ordermodular, modular")
			("structure-prior",        value<string>(&structurePriorFunction)
			                           ->default_value("uniform"),
			                           "set structure prior structure factor"
			                           ", possible values: uniform, indegree-uniform")
			("order-prior",            value<string>(&orderPriorFunction)
			                           ->default_value("uniform"),
			                           "set structure prior order factor"
			                           ", possible values: uniform")
			("max-indegree,m",         value<NoneOrInt>(&maxIndegree)
			                           ->default_value_text("none"),
				                       "set maximum indegree"
				                       ", possible values: none, <non-negative integer>")
			("score-type",             value<string>(&scoreType)
			                           ->default_value("bdeu"),
				                       "set score function"
				                       ", possible values: k2, bdeu, ll, mdl, aic")
			("score-param",            value<double>(&equivalentSampleSize)
			                           ->default_value(1),
				                       "set score parameter (equivalent sample size for BDeu score)")
			;
		
		methodOpts.add_options()
//			("task,T",                 value<std::vector<std::string>>(&taskType),
			("task,T",                 value<std::vector<std::string>>(&taskType)
			                           ->validator([&] (const string& s) {
			                           		return validateSubsetOf(s,
			                           			split("scores,samples,normconst,normconst-harmonic"
			                           				",arcprobs,sampling-time", ","));
			                           }),
				                       "set the computational task(s)"
				                       ", possible values: arcprobs, normconst, normconst-harmonic, samples")
			//("arcs",                   value<string>(&arcsStr),
			//("arcs",                   value<std::vector<Arc>>(&arcs),
			//	                       "set arc(s) for which to compute the probability")
			("method,M",               value<string>(&computationMethod)
			                           ->validator([&] (const string& s) {
			                           		return validateOneOf(s, split("exact,mcmc,mc3,ais,mc", ","));
			                           }),
				                       "set computation method"
				                       ", possible values: exact, mcmc, mc3, ais, mc"
				                       " (use mc if the samples are read from a file)")
			("scoring-pre-gammas",     value<Flag>(&precomputeGammas)
			                           ->default_value(NO),
			                           "precompute values of gamma function for computation of scores")
			("adtree",                 value<Flag>(&useADTree)
			                           ->default_value(NO),
			                           "use ADTree in score computation")
			("adtree-min-count",       value<int>(&adTreeMinCount)
			                           ->default_value(0),
				                       "set minimum count for ADNode before resorting to leaf lists")
			("adtree-max-depth",       value<int>(&adTreeMaxDepth)
			                           ->default_value(0),
				                       "set maximum depth for ADNode before resorting to leaf lists")
			("burnin-fraction",        value<double>(&burninFraction)
			                           ->default_value(0),
				                       "set a constant relative fraction of samples to discard as burn-in")
			("mc-burnin-samples",      value<int>(&mcBurninSteps)
			                           ->default_value(0),
				                       "set a constant number of samples to discard as burn-in")
			("mcmc-step-proposal",     value<string>(&mcmcProposalDist)
			                           ->default_value("swap"),
				                       "set the MCMC step proposal distribution"
				                       ", possible values: swap")
			("mcmc-sample-spacing",     value<int>(&mcmcSampleSpacing)
			                           ->default_value(1),
				                       "set the number of MCMC steps per sample")
			("mcmc-burnin-steps",       value<int>(&mcmcBurninSteps)
			                           ->default_value(0),
				                       "set the number of burn-in steps (ran before starting actual sampling)")
			("mc3-levels",             value<int>(&mc3NLevels)
			                           ->default_value(0),
				                       "set the number of intermediate annealing levels")
			("mc3-level-scheme",       value<string>(&mc3LevelScheme)
			                           ->default_value("linear"),
				                       "set the scheme of level exponents"
				                       ", possible values: linear, autogeometric, geometric")
			("mc3-level-param",        value<double>(&mc3LevelParam)
			                           ->default_value(0),
				                       "set the level scheme parameter")
			("mc3-step-proposal",      value<string>(&mc3ProposalDist)
			                           ->default_value("swap"),
				                       "set the MC step proposal distribution"
				                       ", possible values: swap")
			("mc3-sample-spacing",     value<int>(&mc3SampleSpacing)
			                           ->default_value(1),
				                       "set the number of steps (on all chains) per sample")
			("mc3-level-swaps",        value<int>(&mc3LevelSwapsPerStep)
			                           ->default_value(1),
				                       "set the number of level swaps to try per step")
			("mc3-burnin-steps",       value<int>(&mc3BurninSteps)
			                           ->default_value(0),
				                       "set the number of burn-in steps (ran before starting actual sampling)")
			("ais-levels",             value<int>(&aisNLevels)
			                           ->default_value(0),
				                       "set the number of annealing levels")
				                       //", examples: 100, 2D, 0.5D")
			("ais-level-scheme",       value<string>(&aisLevelScheme)
			                           ->default_value("linear"),
				                       "set the scheme of level exponents"
				                       ", possible values: linear, autogeometric, geometric")
			("ais-level-param",        value<double>(&aisLevelParam)
			                           ->default_value(0),
				                       "set the level scheme parameter")
			("ais-target-samples",     value<int>(&aisTargetDistSamples)
			                           ->default_value(1),
				                       "set the number of samples to draw per annealing")
			("ais-sample-spacing",       value<int>(&aisSampleSpacing)
			                           ->default_value(1),
				                       "set the number of steps per sample (in the target distribution)")
			("ais-step-proposal",      value<string>(&aisProposalDist)
			                           ->default_value("swap"),
				                       "set the MC step proposal distribution"
				                       ", possible values: swap")
			("sample-type",             value<string>(&sampleType)
			                           ->default_value("bucketorder"),
				                       "set partial order type"
				                       ", possible values: bucketorder, order")
			("sample-bucket-size",     value<AutoInt>(&maxBucketSize)
			                           ->default_value_text("auto"),
				                       "set (maximum) bucket size")
			("samples,s",              value<int>(&nSamples)
			                           ->default_value(0),
				                       "set the number of samples to draw")
			("sampling-time,t",        value<double>(&samplingTimeTarget)
			                           ->default_value_text("")
			                           ->validator([&] (const string& s) -> double { 
			                           		if (s.empty()) return 0;
			                           		else return validateDuration(s);
			                           }),
				                       "set a total sampling time")
			("sample-allow-interrupt", value<Flag>(&interruptMidSample)
			                           ->default_value(YES),
				                       "allow interrupting in the middle of current sample")
			("subsample-dags",         value<AutoFlag>(&subSampleDagsExplicit)
			                           ->default_value_text("auto"),
				                       "sample DAG subsamples from orders/bucketorders"
				                       ", possible values: auto, yes, no")
			("subsamples",             value<int>(&nSubSamples)
			                           ->default_value(1),
				                       "set the number of subsamples per sample")
			("subsampling-time",       value<double>(&subSamplingTimeTarget)
			                           ->default_value_text("")
			                           ->validator([&] (const string& s) -> double { 
			                           		if (s.empty()) return 0;
			                           		else return validateDuration(s);
			                           }),
				                       "set a total subsampling time per sample")
			("subsampling-timing-fraction",value<double>(&subSamplingTimingFraction)
			                           ->default_value(0),
				                       "set a fraction subsampling time to use for deciding the number of subsamples")
			("subsampling-correction", value<string>(&subSamplingType)
			                           ->default_value("weighted"),
				                       "subsampling type: weighted, ellis-wong")
			("subsampling-ew-eps",       value<double>(&subSamplingTotalProbFractionEps)
			                           ->default_value(0.05),
				                       "set the epsilon parameter of Ellis-Wong method")
			("subsampling-ew-max-dags", value<size_t>(&maxTotalEWSubSamples)
			                           ->default_value_text("")
			                           ->validator([&] (const string& s) -> double { 
			                           		if (s.empty()) return 0;
			                           		else return validateSize(s);
			                           }),
				                       "set an upper limit for the total number of unique DAGs (Ellis-Wong method)")
			("subsampling-ew-max-mem",  value<size_t>(&maxTotalEWSubSampleMem)
			                           ->default_value_text("")
			                           ->validator([&] (const string& s) -> double { 
			                           		if (s.empty()) return 0;
			                           		else return validateSize(s);
			                           }),
				                       "set an upper limit for the memory used by stored unique DAGs"
				                       " (the actual consumption can be higher due to the overhead of the hash table)"
				                       " (Ellis-Wong method)"
				                       ", for example: 500M, 3.5G or 0.8T")
			;

		outputOpts.add_options()
			("score-out-file",         value<string>(&scoreOutFilename)
			                           ->default_value(""),
				                       "set score input file")

			("output-spacing",         value<int>(&estimateOutputSpacing)
			                           ->default_value(1),
				                       "set an interval to output the estimates"
				                       "(i.e. output the estimates only on every nth sample)")
			("normconst-out-file,o",   value<string>(&normConstOutFilename)
			                           ->default_value("-"),
				                       "specify an output file for the normalizing constant (estimate)")
			("normconst-out-type",     value<string>(&normConstOutType)
			                           ->default_value("final"),
				                       "specify how to output the normalizing constant estimate"
				                       ", possible values: final, update, cumulate")
			("normconst-h-out-file,o", value<string>(&normConstHarmonicOutFilename)
			                           ->default_value("-"),
				                       "specify an output file for the harmonic normalizing constant estimate")
			("normconst-h-out-type",   value<string>(&normConstHarmonicOutType)
			                           ->default_value("final"),
				                       "specify how to output the harmonic normalizing constant estimate"
				                       ", possible values: final, update, cumulate")
			("arcprobs-out-file,o",    value<string>(&arcProbsOutFilename)
			                           ->default_value("-"),
				                       "specify an output file for arc probabilities (estimates)")
			("arcprobs-out-type",      value<string>(&arcProbsOutType)
			                           ->default_value("final"),
				                       "specify how to output the arc probability estimates"
				                       ", possible values: final, update, cumulate")
			("arcprobs-out-format",    value<string>(&arcProbsOutFormat)
			                           ->default_value("pretty"),
				                       "specify how to format the arc probability estimates"
				                       ", possible values: pretty, plain-row")
			("feature-out-file,o",     value<string>(&featureOutFilename)
			                           ->default_value("-"),
				                       "specify an output file for feature probability (estimate)")
			("feature-out-type",       value<string>(&featureOutType)
			                           ->default_value("final"),
				                       "specify how to output the feature probability estimate"
				                       ", possible values: final, update, cumulate")

			("sample-out-file",        value<string>(&sampleOutFilename)
			                           ->default_value(""),
				                       "specify an output file for samples")
			("sample-normconst-out-file",value<string>(&sampleNormConstRatioOutFilename)
			                           ->default_value(""),
				                       "specify an output file for per-sample normalizing constant ratios")
			("sample-weight-out-file", value<string>(&sampleWeightOutFilename)
			                           ->default_value(""),
				                       "specify an output file for per-sample importance weights")
			("subsample-out-file",     value<string>(&subSampleOutFilename)
			                           ->default_value(""),
				                       "specify an output file for subsamples")
			/*("mcmc-accept-ratio-out-file",  value<string>(&mcmcAcceptRatioOutFilename)
			                           ->default_value(""),
				                       "specify an output file for Metrpollis-Hastings proposal acceptance ratios")
			("mcmc-accept-ratio-out-type",  value<string>(&mcmcAcceptRatioOutType)
			                           ->default_value("final"),
				                       "specify how to output the Metrpollis-Hastings proposal acceptance ratios")
			("mc3-accept-ratio-out-file",  value<string>(&mc3AcceptRatioOutFilename)
			                           ->default_value(""),
				                       "specify an output file for Metrpollis-Hastings proposal acceptance ratios")
			("mc3-accept-ratio-out-type",  value<string>(&mc3AcceptRatioOutType)
			                           ->default_value("final"),
				                       "specify how to output the Metrpollis-Hastings proposal acceptance ratios")
			("ais-accept-ratio-out-file",  value<string>(&aisAcceptRatioOutFilename)
			                           ->default_value(""),
				                       "specify an output file for Metrpollis-Hastings proposal acceptance ratios")
			("ais-accept-ratio-out-type",  value<string>(&aisAcceptRatioOutType)
			                           ->default_value("final"),
				                       "specify how to output the Metrpollis-Hastings proposal acceptance ratios")*/
			("multilevel-weight-file", value<string>(&levelWeightOutFilename)
			                           ->default_value(""),
				                       "specify an output file for per-sample multilevel importance weights")
			("multilevel-weight-spacing",value<int>(&levelWeightOutSpacing)
			                           ->default_value(1),
				                       "specify multilevel importance weight spacing")
			("multilevel-inv-weight-file",value<string>(&levelInvWeightOutFilename)
			                           ->default_value(""),
				                       "specify an output file for per-sample multilevel harmonic importance weights")
			("multilevel-inv-weight-spacing",value<int>(&levelInvWeightOutSpacing)
			                           ->default_value(1),
				                       "specify multilevel harmonic importance weight spacing")
			("sample-prob-out-file",   value<string>(&sampleProbOutFilename)
			                           ->default_value(""),
				                       "specify an output file for unnormalized sample probabilities")
			("subsample-count-out-file",value<string>(&subSampleCountOutFilename)
			                           ->default_value(""),
				                       "specify an output file for per-sample subsample counts")
			("sampling-time-out-file", value<string>(&samplingTimeOutFilename)
			                           ->default_value("-"),
				                       "specify an output file for total sampling time")
			("sampling-time-out-type", value<string>(&samplingTimeOutType)
			                           ->default_value("final"),
				                       "specify how to output the total sampling time"
				                       ", possible values: final, update, cumulate")
			("log-out-file",           value<string>(&logOutFilename)
			                           ->default_value(""),
				                       "specify an output file for logging")
			;

		//opts::positional_options_description pdesc;
		//pdesc.add("data-in-file", 1);
		//pdesc.add("estimate-out-file", 1);

		allOpts.add(generalOpts);
		allOpts.add(inputOpts);
		allOpts.add(modelOpts);
		allOpts.add(methodOpts);
		allOpts.add(outputOpts);
	}
	
	void parse(int argc, char** argv) {
		try {
			//opts::store(opts::command_line_parser(argc, argv).options(desc).positional(pdesc).run(), vm);
			opts::store(opts::command_line_parser(argc, argv).options(allOpts).run(), vm);
			opts::notify(vm);
			if (!configFilename.empty()) {
				std::ifstream configFile(configFilename.c_str());
				opts::store(opts::parse_config_file(configFile, allOpts), vm);
				opts::notify(vm);
			}
		}
		catch (opts::error& err) {
			throw Exception("Error: %s", err.what());
		}

		//showHelp = vm.count("help");
		//showVersion = vm.count("version");
		
		if (showHelp || showVersion)
			return;

		//quiet = vm.count("quiet");
		//precomputeGammas = vm.count("scoring-pre-gammas");
		//useADTree = vm.count("adtree");
		
		// task
		/*for (auto tt : split(vm["task"].as<std::string>(), ",")) {
			if (!contains(split("scores,samples,normconst,normconst-harmonic,arcprobs,sampling-time", ","), tt))
				throw Exception("Invalid task type '%s'", tt);
			taskType.push_back(tt);
		}*/

		// subsampling
		if (subSampleDagsExplicit.isSpecial()) {
			subSampleDags = false;
			if (sampleType == "bucketorder" || sampleType == "order")
				if (priorType == "modular") {
					subSampleDags = true;
				}
		}
		else {
			//subSampleDags = toBoolean(subSampleDagsExplicit);
			subSampleDags = subSampleDagsExplicit.value();
		}

		/*// sampling time
		if (samplingTimeTargetStr.empty()) {
			samplingTimeTarget = 0;
		}
		else {
			try {
				samplingTimeTarget = parseDuration(samplingTimeTargetStr);
			} catch (Exception& e) {
				throw Exception("Failed to parse sampling time: %s", e.what());
			}
		}

		// subsampling time
		if (subSamplingTimeTargetStr.empty()) {
			subSamplingTimeTarget = 0;
		}
		else {
			try {
				subSamplingTimeTarget = parseDuration(subSamplingTimeTargetStr);
			} catch (Exception& e) {
				throw Exception("Failed to parse subsampling time: %s", e.what());
			}
		}

		// EW subsampling max dags
		if (maxTotalEWSubSamplesStr.empty()) {
			maxTotalEWSubSamples = 0;
		}
		else {
			try {
				maxTotalEWSubSamples = parseSize(maxTotalEWSubSamplesStr);
			} catch (Exception& e) {
				throw Exception("Failed to parse sample-subsamples-ew-max-dags: %s", e.what());
			}
		}

		// EW subsampling max mem
		if (maxTotalEWSubSampleMemStr.empty()) {
			maxTotalEWSubSampleMem = 0;
		}
		else {
			try {
				maxTotalEWSubSampleMem = parseSize(maxTotalEWSubSampleMemStr);
			} catch (Exception& e) {
				throw Exception("Failed to parse sample-subsamples-ew-max-mem: %s", e.what());
			}
		}*/
	}
	
};

#endif
