/*
 * TBAmDownsampler.cpp
 *
 *  Created on: May 31, 2020
 *      Author: phaentu
 */

#include "TBamDownsampler.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Main/TRandomGenerator.h"
#include "coretools/Strings/concatenateString.h"


namespace GenomeTasks{

using coretools::Probability;
using coretools::instances::parameters;
using coretools::instances::logfile;
using coretools::instances::randomGenerator;
using coretools::str::toString;

//-----------------------------------------
// TBamSample
//-----------------------------------------
TBamSample::TBamSample(const Probability &Prob, const std::string &OutName, BAM::TBamFile &bamFile)
	: _prob(Prob), _outName(OutName), _out(_outName, bamFile){};

void TBamSample::sample(BAM::TBamFile & bamFile){
	const auto& read = bamFile.curRead();
	if(_discard.isInBlacklist(read.name())){
		_discard.remove(read.name());
	} else if(_keep.isInBlacklist(read.name())){
		bamFile.writeCurAlignment(_out);
		_keep.remove(read.name());
	} else if(randomGenerator().getRand() < _prob){
		bamFile.writeCurAlignment(_out);
		if(read.isProperPair()){
			_keep.add(read.name());
		}
	} else {
		//filtered out
		if(read.isProperPair()){
			_discard.add(read.name());
		}
	}
};

void TBamSample::downsampleRead(BAM::TAlignment & alignment){
	//parse again to get original bases and qualities
	alignment.parse();

	//downsample
	alignment.downsampleAlignment(_prob);

	//write
	_out.writeAlignment(alignment);
};

//-----------------------------------------
// TBamDownsampler_base
//-----------------------------------------

void TBamDownsampler::_readVectorOfDownsamplingProbabilities(){
    //read downsampling rates
    if(parameters().exists("prob")) {
        parameters().fill("prob", _probs);
    } else if(parameters().exists("depth")){
        std::vector<double> depths;
		parameters().fill("depth", depths);
		double averageDepth;
		if (parameters().exists("averageDepth")) {
			averageDepth = parameters().get<double>("averageDepth");
		} else {
			logfile().list("No averageDepth given, will calculate it. Use 'averageDepth' to safe time!");
			averageDepth = _genome.bamFile().averageDepth();
			logfile().list("Average depth estimated to ", averageDepth);
		}
		for (auto &it : depths) {
			if (averageDepth >= it) {
				_probs.emplace_back(it / averageDepth);
			} else {
				throw coretools::TUserError("Average Depth must be equal or bigger than provided lists of depths");
			}
		}
	} else {
        throw coretools::TUserError("Either argument 'prob' or 'depth' must be provided!");
    }
	//get unique names
	std::map <Probability, int> fracNames;
	for(size_t i=0; i<_probs.size(); ++i){
		std::map<Probability, int>::iterator it = fracNames.find(_probs[i]);
		if(it == fracNames.end()){
			fracNames.emplace(_probs[i],1);
			_names.push_back(toString(_probs[i]));
		} else {
			++(it->second);
			_names.push_back(toString(_probs[i]) + "_" + toString(it->second));
		}
	}
};

//-----------------------------------------
// TBamDownsampler
//-----------------------------------------
TBamDownsampler::TBamDownsampler() {
	//read downsampling rates
	_readVectorOfDownsamplingProbabilities();

	if (*_probs.begin() == 1.0) logfile().warning("Probability of 1 will result in identical file!");
	std::string filePrefix;
	if (parameters().exists("downsampleBases")) {
		logfile().list("Will downsample by removing bases (i.e. replacing random bases with Ns). (parameter 'downsampleBases')");
		_writeN = true;
	} else {
		logfile().list("Will downsample by removing reads. (use 'downsampleBases' to remove bases)");
		_writeN = false;
	}

	if (!_writeN && parameters().exists("separateReads")) {
		filePrefix = _genome.outputName() + "_separated_";
		// report
		logfile().list("Will separate reads with probabilities (parameter 'prob'): " +
					   coretools::str::concatenateString(_probs, ", "));

		separateReads = true;
		// check that sum <= 1.0
		double sum    = 0.0;
		for (auto &d : _probs) {
			sum += d.get();
			_cumulProbs.push_back(sum);
		}
		_cumulProbs.push_back(1.0); // always add an extra at end to ease search

		coretools::user_assert(sum <= 1.0, "Separation probabilities must sum to <= 1.0, not ", sum, "!");
	} else {
		filePrefix = _genome.outputName() + "_downsampled_";
		// report
		logfile().list("Will accept reads with probabilities (parameter 'prob'): " +
					   coretools::str::concatenateString(_probs, ", "));
	}

	// create downsampling objects
	_bamSamples.reserve(_probs.size());
	for (size_t i = 0; i < _probs.size(); ++i) {
		std::string filename = filePrefix + _names[i] + ".bam";
		_bamSamples.emplace_back(_probs[i], filename, _genome.bamFile());
	}
};

void TBamDownsampler::run() {
	// traverse BAM and downsample
	_genome.bamFile().startProgressReporting();
	while (_genome.bamFile().readNextAlignment()) {
		if (separateReads) {
			sample();
		} else if (_writeN) {
			BAM::TAlignment alignment;
			_genome.bamFile().fill(alignment);
			for (auto &s : _bamSamples) { s.downsampleRead(alignment); }
		} else {
			for (auto &s : _bamSamples) { s.sample(_genome.bamFile()); }
		}
		_genome.bamFile().printProgress();
	}
	_genome.bamFile().printEndWithSummary(_genome.outputName());
};

void TBamDownsampler::sample(){
	const auto& read = _genome.bamFile().curRead();
	if(_discard.isInBlacklist(read.name())){
		_discard.remove(read.name());
	} else {
		auto mate = _mateWasWritten.find(read.name());
		if(mate != _mateWasWritten.end()){
			_genome.bamFile().writeCurAlignment(_bamSamples[mate->second]._out);
			_mateWasWritten.erase(mate);
		} else {
			//assing to a bam file
			double r = randomGenerator().getRand();

			size_t index = 0;
			while(r > _cumulProbs[index]){
				++index;
			}
			if(index < _bamSamples.size()){
				//write
				_genome.bamFile().writeCurAlignment(_bamSamples[index]._out);
				if(read.isProperPair()){
					_mateWasWritten.emplace(read.name(), index);
				}
			} else {
				//discard read
				if(read.isProperPair()){
					_discard.add(read.name());
				}
			}
		}
	}
}

}; // end namespace
