#include "TSpearmanGWAS.h"

#include "coretools/Files/TOutputFile.h"
#include "coretools/Files/TInputFile.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Math/mathFunctions.h"
#include "coretools/Main/TRandomGenerator.h"

#include "genometools/DefaultColumnNames.h"
#include "genometools/TSampleLikelihoods.h"


#include <map>

namespace PopulationTools{

using coretools::TOutputFile;
using coretools::TInputFile;
using coretools::instances::logfile;
using coretools::instances::parameters;

//------------------------------------------------
//TSpearmanGWASPopulation
//------------------------------------------------

void TSpearmanGWASPopulation::resize(size_t Size){
	clear();
	_data.reserve(Size);
	_dosage.reserve(Size);
}

void TSpearmanGWASPopulation::clear(){
	_data.clear();
	_dosage.clear();
}

void TSpearmanGWASPopulation::addSample(double Data) {
	_data.emplace_back(Data);
}

void TSpearmanGWASPopulation::finalizeDataReading(){
	_ranksData = coretools::ranks(_data);

	// standardize, calculate RSS of total model and prepare XT_X_inv_X for regression
	coretools::standardizeZeroMeanUnitVar(_ranksData);
	_RSStotal = coretools::sumOfSquares(_ranksData);
	_XT_X_inv_X.resize(_ranksData.size());
	for(size_t i = 0; i < _ranksData.size(); ++i){
		_XT_X_inv_X[i] = _ranksData[i] / _RSStotal;
	}

	// prepare storage for dosage
	_dosage.resize(_data.size());
}

void TSpearmanGWASPopulation::prepareBootstraps(const size_t NumBootstraps){
	_bootstraps.resize(NumBootstraps);
	for(size_t i = 0; i < NumBootstraps; ++i){
		_bootstraps[i].resize(_data.size());
		coretools::instances::randomGenerator().fillRandomOrderOfIndexes(_bootstraps[i]);
	}
}

void TSpearmanGWASPopulation::updateDosage(genometools::TSampleLikelihoods<genometools::HighPrecisionPhredIntProbability> *GenotypeLikelihoods){
	// calculate dosage for each sample 
	for(size_t i = 0; i < _data.size(); ++i){
		_dosage[i] = GenotypeLikelihoods[i].meanPosteriorGenotype();
	}
}

double TSpearmanGWASPopulation::RSS_nullModel() { 
	return _RSStotal;
}

double TSpearmanGWASPopulation::_regressAndCalulateRSSModel(const std::vector<double> & ranksDosage){
	// run regression
	double beta = coretools::sumPairwiseProduct(_XT_X_inv_X, ranksDosage);
	
	// calculate RSS for model
	double RSS_model = 0.0;
	for(size_t i = 0; i < _ranksData.size(); ++i){
		double tmp = _ranksData[i] - beta * ranksDosage[i];
		RSS_model += tmp * tmp;
	}
	return RSS_model;
}		

double TSpearmanGWASPopulation::_regressAndCalulateRSSModelBootstrap(size_t bootstrapIdx, const std::vector<double> & ranksDosage){
	// run regression
	double beta = 0.0;
	for(size_t i = 0; i < _ranksData.size(); ++i){
		beta += _XT_X_inv_X[i] * ranksDosage[ _bootstraps[bootstrapIdx][i] ];
	}

	// calculate RSS for model
	double RSS_model = 0.0;
	for(size_t i = 0; i < _ranksData.size(); ++i){
		double tmp = _ranksData[i] - beta * ranksDosage[ _bootstraps[bootstrapIdx][i] ];
		RSS_model += tmp * tmp;
	}
	return RSS_model;
}	

double TSpearmanGWASPopulation::regressSpearmanAndBootstrap(std::vector<double> & bootstrapsRSSModel){
	std::vector<double> ranksDosage = coretools::ranks(_dosage);
	coretools::standardizeZeroMeanUnitVar(ranksDosage);

	// fill bootstraps
	for(size_t b = 0; b < _bootstraps.size(); ++b){
		bootstrapsRSSModel[b] += _regressAndCalulateRSSModelBootstrap(b, ranksDosage);	
	}

	// return model RSS
	return _regressAndCalulateRSSModel(ranksDosage);
}


//------------------------------------------------
// TSpearmanGWAS
//------------------------------------------------

// read data
std::map<std::string, double> TSpearmanGWAS::_readDataIntoMap(std::string_view Filename){	
	TInputFile in(Filename, coretools::FileType::Header);

	// get sample and data columns
	size_t sampleIdx = in.indexOfFirstMatch(genometools::defaultColumnNames_sample);
	int dataIdx = -1;
	if(parameters().exists("dataCol")){
		// user provided column name
		auto dataColName = parameters().get<std::string>("dataCol");
		logfile().list("Reading data from column '", dataColName, "' of file '", Filename, "'. (parameter 'dataCol')");		
		dataIdx = in.index(dataColName);
	} else {
		// use first column in file that is not sample nor population
		for(size_t i = 0; i < in.numCols(); ++i){
			if(std::find(genometools::defaultColumnNames_sample.begin(), genometools::defaultColumnNames_sample.end(), in.header()[i]) == genometools::defaultColumnNames_sample.end()
			&& std::find(genometools::defaultColumnNames_population.begin(), genometools::defaultColumnNames_population.end(), in.header()[i]) == genometools::defaultColumnNames_population.end()){
				//column is not a sample or population
				dataIdx = i;								
				logfile().list("Reading data from column '", in.header()[i], "' of file '", Filename, "'. (specify column with 'dataCol') ");
				break;				
			}			
		}

		// if no column was found: error!
		if(dataIdx < 0){
			UERROR("Sample file '", Filename, "' lacks a data column!");
		}
	}
	
	// read into map
	std::map<std::string, double> data;
	for(; !in.empty(); in.popFront()){
		data.emplace(in.get(sampleIdx), in.get<double>(dataIdx));
	}

	return data;
}


TSpearmanGWAS::TSpearmanGWAS(){	
	//read samples
	std::string sampleFile = parameters().get<std::string>("samples");
	_samples.readSamples(sampleFile);	

	//open VCF
	_vcfFilename = parameters().get<std::string>("vcf");
	logfile().list("Reading vcf from file '" + _vcfFilename + "'.");
	_vcfFile.openStream(_vcfFilename);

	//enable parsers
	_vcfFile.enablePositionParsing();
	_vcfFile.enableVariantParsing();
	_vcfFile.enableVariantQualityParsing();
	_vcfFile.enableFormatParsing();
	_vcfFile.enableSampleParsing();

	//Match samples
	_samples.fillVCFOrder(_vcfFile.parser.samples);		

	// read data into map
	auto data = _readDataIntoMap(sampleFile);
		
	// create populations and add samples in correct order
	_populations.resize(_samples.numPopulations());		
	for(size_t s = 0; s < _samples.numSamples(); ++s){
		_populations[_samples.populationIndex(s)].addSample(data[_samples.sampleName(s)]);
	}

	const auto tmp     = coretools::str::readBeforeLast(_vcfFilename, ".vcf");
	auto tag = parameters().get("out", tmp);
	_outname = tag + "_SpearmanGWAS.txt.gz";
	logfile().list("Writing Spearman GWAS results to file '" + _outname + "'. (parameter 'out')");

	//limit lines?
	_limitLines = parameters().exists("limitLines");
	if(_limitLines){
		_maxNumLines = parameters().get<long>("limitLines");
	} else {
		_maxNumLines = 0;
	}

	// num bootstraps
	_numBootstraps = parameters().get<int>("bootstraps", 1000);
	logfile().listFlush("Preparing ", _numBootstraps, " bootstraps ... ");
	for(auto& p: _populations){
		p.finalizeDataReading();
		p.prepareBootstraps(_numBootstraps);
	}
	logfile().write("done! (parameter 'bootstraps')");
};

double TSpearmanGWAS::_calculateF(const double RSS_null, const double RSS_model, const double F_scale){
	return (RSS_null - RSS_model) / RSS_model * F_scale;
}

void TSpearmanGWAS::run(){
	// open VCF reader	
	genometools::TPopulationLikelihoodReaderLocus reader(false);
	reader.openVCF(_vcfFilename);
	logfile().endIndent();

	//open output file	
	std::vector<std::string> header = {genometools::defaultColumnNames_chromosome[0], genometools::defaultColumnNames_position[0], genometools::defaultColumnNames_reference[0], genometools::defaultColumnNames_alternative[0], "F", "p"};	
	TOutputFile out(_outname, header);	

	// calculate RSS_null and F_scale
	double RSS_null = 0.0;
	for(size_t p = 0; p < _populations.size(); ++p){
		RSS_null += _populations[p].RSS_nullModel();
	}
	double F_scale = (double) (_samples.numSamples() - _samples.numPopulations()) / (double) _samples.numPopulations();

	//progress
	coretools::TTimer timer;	

	// likelihood and tmp storage
	genometools::TPopulationLikehoodLocus< genometools::TSampleLikelihoods<genometools::HighPrecisionPhredIntProbability > > data(_samples.numSamples());
	std::vector<double> bootstrappedRSS_model(_numBootstraps);

	// run through VCF file
	logfile().startIndent("Parsing VCF file and perform Spearman GWAS:");
	while (reader.readDataFromVCF(data, _samples)) {
		// update dosage for each sample in populations, calculate rho and bootstrap
		double RSS_model = 0.0;
		std::fill(bootstrappedRSS_model.begin(), bootstrappedRSS_model.end(), 0.0);
	
		for(size_t p = 0; p < _populations.size(); ++p){
			_populations[p].updateDosage(&data[_samples.startIndex(p)]);
			RSS_model += _populations[p].regressSpearmanAndBootstrap(bootstrappedRSS_model);			
		}

		// calculate F
		double F = _calculateF(RSS_null, RSS_model, F_scale);

		// determine p value from bootstraps
		size_t nBelow = 0;
		for(size_t b = 0; b < _numBootstraps; ++b){
			double F_bootstrap = _calculateF(RSS_null, bootstrappedRSS_model[b], F_scale);
			if(F_bootstrap <= F){
				++nBelow;
			}
		}
		double p = (double) nBelow / (double) _numBootstraps;

		// output
		reader.writePosition(out);
		out.write(F, p);
		out.endln();		
		if(_limitLines && out.curLine() >= _maxNumLines){
			break;
		}		
	}

	// report final status
	logfile().endIndent();
	reader.concludeFilters();
	logfile().endIndent();
};



}; // namespace PopulationTools
