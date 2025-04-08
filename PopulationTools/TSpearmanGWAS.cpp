#include "TSpearmanGWAS.h"

#include "coretools/Files/TOutputFile.h"
#include "coretools/Files/TInputFile.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Math/mathFunctions.h"
#include "coretools/Main/TRandomGenerator.h"
#include "coretools/Math/mathFunctions.h"

#include "coretools/algorithms.h"
#include "genometools/DefaultColumnNames.h"
#include "genometools/VCF/TPopulationLikelihoods.h"


#include <map>

namespace PopulationTools{

using coretools::TOutputFile;
using coretools::TInputFile;
using coretools::instances::logfile;
using coretools::instances::parameters;

namespace SpearmanGWASimpl{

//--------------------------------------------
// TSpearmanGWASBootstraps
//--------------------------------------------
TSpearmanGWASBootstraps::TSpearmanGWASBootstraps(size_t SampleSize, size_t NumBootstraps){
	_bootstraps.resize(NumBootstraps);
	for(size_t i = 0; i < NumBootstraps; ++i){
		_bootstraps[i].resize(SampleSize);
		coretools::instances::randomGenerator().fillRandomOrderOfIndexes(_bootstraps[i]);
	}
}

TSpearmanGWASBootstrapLibrary::TSpearmanGWASBootstrapLibrary(size_t NumBootstraps):
	_numBootstraps(NumBootstraps){};

const TSpearmanGWASBootstraps& TSpearmanGWASBootstrapLibrary::get(size_t SampleSize) {
	// returnbootstrap object, create if missing
	auto lib = _bootstraps.find(SampleSize);
	if(lib == _bootstraps.end()){
		_bootstraps.emplace(std::pair<size_t, TSpearmanGWASBootstraps>(
			std::piecewise_construct,
              std::forward_as_tuple(SampleSize),
              std::forward_as_tuple(SampleSize, _numBootstraps)));
		lib = _bootstraps.find(SampleSize);
	}
	return lib->second;
}

//------------------------------------------------
//TSpearmanGWASPopulation
//------------------------------------------------

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
	_dataWithData.resize(_data.size());
}

void TSpearmanGWASPopulation::updateDosage(genometools::TSampleLikelihoods<coretools::HPPhredInt> *GenotypeLikelihoods){

	// calculate dosage for each sample 
	_sampleSize = 0;
	for(size_t i = 0; i < _data.size(); ++i){
		if(!GenotypeLikelihoods[i].isMissing()){
			_dosage[_sampleSize] = GenotypeLikelihoods[i].meanPosteriorGenotype();
			_dataWithData[_sampleSize] = _data[i];
			++_sampleSize;
		}
	}

	//resize vectors, does not actually free storage so it can be reused for next locus.
	_dataWithData.resize(_sampleSize);
	_dosage.resize(_sampleSize);
}	

void TSpearmanGWASPopulation::regressSpearmanAndAddRSS(double& RSS_null, double& RSS_model){
	// rank, standardize, calculate RSS of total model and prepare XT_X_inv_X for regression
	_ranksData = coretools::ranks(_dataWithData);
	coretools::standardizeZeroMeanUnitVar(_ranksData);
	_XT_X_inv_X.resize(_ranksData.size());
	for(size_t i = 0; i < _ranksData.size(); ++i){
		_XT_X_inv_X[i] = _ranksData[i] / _RSStotal;
	}

	// rank and standardize dosage
	_ranksDosage = coretools::ranks(_dosage);
	coretools::standardizeZeroMeanUnitVar(_ranksDosage);

	// perform regression
	double beta = coretools::sumPairwiseProduct(_XT_X_inv_X, _ranksDosage);
	
	// calculate RSS for null and model
	RSS_null += coretools::sumOfSquares(_ranksData);
	for(size_t i = 0; i < _ranksData.size(); ++i){
		double tmp = _ranksData[i] - beta * _ranksDosage[i];
		RSS_model += tmp * tmp;
	}
}

void TSpearmanGWASPopulation::_regressSpearmanAndAddRSSBootstrap(double& RSS_null, double& RSS_model, const std::vector<size_t>& Bootstraps){
	// perform regression: sum pairwise product with bootstrap
	double beta = 0.0;
	for(size_t i = 0; i < _sampleSize; ++i){
		beta += _XT_X_inv_X[i] * _ranksDosage[ Bootstraps[i] ];
	}
	
	// calculate RSS for null and model
	RSS_null += coretools::sumOfSquares(_ranksData);
	for(size_t i = 0; i < _ranksData.size(); ++i){
		double tmp = _ranksData[i] - beta * _ranksDosage[ Bootstraps[i] ];
		RSS_model += tmp * tmp;
	}
}

void TSpearmanGWASPopulation::bootstrap(std::vector<double>& RSS_null, std::vector<double>& RSS_model, TSpearmanGWASBootstrapLibrary& BootstrapLib){
	// get boostraps
	const TSpearmanGWASBootstraps& bootstraps = BootstrapLib.get(_sampleSize);

	// conduct bootstraps
	for(size_t b = 0; b < BootstrapLib.numBootstraps(); ++b){
		_regressSpearmanAndAddRSSBootstrap(RSS_null[b], RSS_model[b], bootstraps.get(b));
	}
}

} // end namespace SpearmanGWASimpl

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

	// finalize data reading
	for(auto& p : _populations){
		p.finalizeDataReading();
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
	_bootstrapPThreshold = parameters().get<double>("bootstrapThreshold", 0.01);
	if(_bootstrapPThreshold <= 0.0){
		logfile().list("Will not conduct any bootstraps. (set 'bootstrapThreshold' > 0.0 to conduct bootstraps)");
	} else {
		_numBootstraps = parameters().get<int>("bootstraps", 1000);
		if(_numBootstraps == 0){
			logfile().list("Will not conduct any bootstraps. (set 'bootstraps' > 0 to conduct bootstraps)");
		} else  {
			logfile().list("Will assess p-values using  ", _numBootstraps, " bootstraps if the F-based p-value is < ", _bootstrapPThreshold, ". (parameters 'bootstraps' and 'bootstrapThreshold')");
		}
	}
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
	std::vector<std::string> header = {genometools::defaultColumnNames_chromosome[0], genometools::defaultColumnNames_position[0], genometools::defaultColumnNames_reference[0], genometools::defaultColumnNames_alternative[0], "F", "p", "p_bootstrap"};	
	TOutputFile out(_outname, header);	

	// prepare bootstrap library
	SpearmanGWASimpl::TSpearmanGWASBootstrapLibrary bootstrapLib(_numBootstraps);	

	//progress
	coretools::TTimer timer;	

	// likelihood and tmp storage
	genometools::TPopulationLikehoodLocus< genometools::TSampleLikelihoods<coretools::HPPhredInt > > data(_samples.numSamples());
	std::vector<double> bootstrappedRSS_null(_numBootstraps);
	std::vector<double> bootstrappedRSS_model(_numBootstraps);

	// run through VCF file
	logfile().startIndent("Parsing VCF file and perform Spearman GWAS:");
	while (reader.readDataFromVCF(data, _samples)) {
		// update dosage for each sample in populations, calculate rho and bootstrap
		double RSS_null = 0.0;
		double RSS_model = 0.0;
		size_t totSamples = 0;
		
		for(size_t p = 0; p < _populations.size(); ++p){
			_populations[p].updateDosage(&data[_samples.startIndex(p)]);
			_populations[p].regressSpearmanAndAddRSS(RSS_null, RSS_model);		
			totSamples += _populations[p].sampleSize();	
		}

		// calculate F and p-value
		double d1 = _populations.size();
		double d2 = totSamples - _populations.size();
		double Fscale = d2 / d1;
		double F = _calculateF(RSS_null, RSS_model, Fscale);
		coretools::Probability x = coretools::Probability((d1*F) / (d1*F + d2));
		double p = coretools::TIncompleteBeta::incompleteBeta(d1 / 2.0, d2 / 2.0, x);

		// output
		reader.writePosition(out);
		out.write(F, p);

		// calculate p-value via bootstrap
		if(p < _bootstrapPThreshold && _numBootstraps > 0){
			std::fill(bootstrappedRSS_null.begin(), bootstrappedRSS_null.end(), 0.0);
			std::fill(bootstrappedRSS_model.begin(), bootstrappedRSS_model.end(), 0.0);
			for(size_t p = 0; p < _populations.size(); ++p){
				_populations[p].bootstrap(bootstrappedRSS_null, bootstrappedRSS_model, bootstrapLib);
			}

			// determine p value from bootstraps
			size_t nBelow = 0;
			for(size_t b = 0; b < _numBootstraps; ++b){
				double F_bootstrap = _calculateF(bootstrappedRSS_null[b], bootstrappedRSS_model[b], Fscale);
				if(F_bootstrap <= F){
					++nBelow;
				}
			}
			double pBootstrap = (double) nBelow / (double) _numBootstraps;
			out.write(pBootstrap);
		} else {
			out.write("-");
		}

		// end line and decide to break if limit lines was reached
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
