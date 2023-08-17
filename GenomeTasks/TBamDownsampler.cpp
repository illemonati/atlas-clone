/*
 * TBAmDownsampler.cpp
 *
 *  Created on: May 31, 2020
 *      Author: phaentu
 */

#include "TBamDownsampler.h"

#include <algorithm>
#include <cstdint>
#include <exception>
#include <map>
#include <stddef.h>
#include <utility>

#include "TAlignment.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Main/TRandomGenerator.h"
#include "coretools/Strings/stringFunctions.h"

namespace GenomeTasks{

using coretools::Probability;
using coretools::instances::parameters;
using coretools::instances::logfile;
using coretools::instances::randomGenerator;
using coretools::str::toString;

//-----------------------------------------
// TBamSample
//-----------------------------------------
TBamSample::TBamSample(const Probability & Prob, const std::string & OutName){
	_prob = Prob;
	_outName = OutName;
};

void TBamSample::open(BAM::TBamFile & bamFile){
	_out.open(_outName, bamFile);
};

void TBamSample::close(){
	_out.close();
};

void TBamSample::sample(BAM::TBamFile & bamFile){
	if(_discard.isInBlacklist(bamFile.curName())){
		_discard.remove(bamFile.curName());
	} else if(_keep.isInBlacklist(bamFile.curName())){
		bamFile.writeCurAlignment(_out);
		_keep.remove(bamFile.curName());
	} else if(randomGenerator().getRand() < _prob){
		bamFile.writeCurAlignment(_out);
		if(bamFile.curIsProperPair()){
			_keep.add(bamFile.curName());
		}
	} else {
		//filtered out
		if(bamFile.curIsProperPair()){
			_discard.add(bamFile.curName());
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

void TBamDownsampler_base::_readVectorOfDownsamplingProbabilities(){
    //read downsampling rates
    if(parameters().parameterExists("prob")) {
        parameters().fillParameterIntoContainerRepeatIndexes("prob", _probs, ',');
    } else if(parameters().parameterExists("depth")){
        std::vector<double> depths;
        parameters().fillParameterIntoContainerRepeatIndexes("depth", depths, ',');
        double averageDepth = parameters().getParameter<double>("averageDepth");
        for(auto& it : depths){
            if(averageDepth >= it){
                _probs.push_back(it / averageDepth);
            } else{
                UERROR("Average Depth must be equal or bigger than provided lists of depths");
            }
        }
    } else {
        UERROR("Either argument 'prob' or 'depth' must be provided!");
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
TBamDownsampler::TBamDownsampler() : TBamDownsampler_base(){
	//read downsampling rates
	_readVectorOfDownsamplingProbabilities();

	if (*_probs.begin() == 1.0) logfile().warning("Probability of 1 will result in identical file!");
	std::string filePrefix;
	if (parameters().parameterExists("downsampleBases")) {
		logfile().list("Will downsample by removing bases (i.e. replacing random bases with Ns). (parameter 'downsampleBases')");
		_writeN = true;
	} else {
		logfile().list("Will downsample by removing reads. (use 'downsampleBases' to remove bases)");
		_writeN = false;
	}

	if (!_writeN && parameters().parameterExists("separateReads")) {
		filePrefix = _outputName + "_separated_";
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

		if (sum > 1.0) { UERROR("Separation probabilities must sum to <= 1.0, not ", sum, "!"); }
	}
	else {
		filePrefix = _outputName + "_downsampled_";
		// report
		logfile().list("Will accept reads with probabilities (parameter 'prob'): " +
					   coretools::str::concatenateString(_probs, ", "));
	}

	// create downsampling objects
	for (size_t i = 0; i < _probs.size(); ++i) {
		std::string filename = filePrefix + _names[i] + ".bam";
		_bamSamples.emplace_back(_probs[i], filename);
	}

	// open bam files for writing
	for (auto &s : _bamSamples) { s.open(_bamFile); }
};

void TBamDownsampler::run() {
	// traverse BAM and downsample
	_bamFile.startProgressReporting();
	while (_bamFile.readNextAlignment()) {
		if (separateReads) {
			sample();
		} else if (_writeN) {
			BAM::TAlignment alignment;
			_bamFile.fill(alignment);
			for (auto &s : _bamSamples) { s.downsampleRead(alignment); }
		} else {
			for (auto &s : _bamSamples) { s.sample(_bamFile); }
		}
		_bamFile.printProgress();
	}
	_bamFile.printEndWithSummary(_outputName);

	// close
	for (auto &s : _bamSamples) { s.close(); }
};

void TBamDownsampler::sample(){
	if(_discard.isInBlacklist(_bamFile.curName())){
		_discard.remove(_bamFile.curName());
	} else {
		auto mate = _mateWasWritten.find(_bamFile.curName());
		if(mate != _mateWasWritten.end()){
			_bamFile.writeCurAlignment(_bamSamples[mate->second]._out);
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
				_bamFile.writeCurAlignment(_bamSamples[index]._out);
				if(_bamFile.curIsProperPair()){
					_mateWasWritten.emplace(_bamFile.curName(), index);
				}
			} else {
				//discard read
				if(_bamFile.curIsProperPair()){
					_discard.add(_bamFile.curName());
				}
			}
		}
	}
}

//-----------------------------------------
// TBamSeparator
//-----------------------------------------
TBamSeparator::TBamSeparator() : TBamDownsampler_base(){
	//read downsampling rates
	_readVectorOfDownsamplingProbabilities();

	//check that sum <= 1.0
	double sum = 0.0;
	for(auto& d : _probs){
		sum += d.get();
		_cumulProbs.push_back(sum);
	}
	_cumulProbs.push_back(1.0); //always add an extra at end to ease search

	if(sum > 1.0){
		UERROR("Separation probabilities must sum to <= 1.0, not ", sum, "!");
	}

	//report
	logfile().list("Will separate reads with probabilities (parameter 'prob'): " + coretools::str::concatenateString(_probs, ", "));
	if(*_probs.begin() == 1.0) logfile().warning("Probability of 1 will result in identical file!");
};

void TBamSeparator::run(){
	//open bam files for writing
	std::vector<BAM::TOutputBamFile> out;
	for(auto& n : _names){
		std::string filename = _outputName + "_downsampled_" + n + ".bam";
		out.emplace_back(filename, _bamFile);
	}

	//Prepare lists to keep track of mates
	std::map<std::string, uint16_t> _mateWasWritten;
	BAM::TAlignmentList _discard;


	//traverse BAM and separate
	_bamFile.startProgressReporting();
	while(_bamFile.readNextAlignment()){
		//check if mate was already written
		if(_discard.isInBlacklist(_bamFile.curName())){
			_discard.remove(_bamFile.curName());
		} else {
			auto mate = _mateWasWritten.find(_bamFile.curName());
			if(mate != _mateWasWritten.end()){
				_bamFile.writeCurAlignment(out[mate->second]);
				_mateWasWritten.erase(mate);
			} else {
				//assing to a bam file
				double r = randomGenerator().getRand();

				size_t index = 0;
				while(r < _cumulProbs[index]){
					++index;
				}
				if(index < out.size()){
					//write
					_bamFile.writeCurAlignment(out[index]);
					if(_bamFile.curIsProperPair()){
						_mateWasWritten.emplace(_bamFile.curName(), index);
					}
				} else {
					//discard read
					if(_bamFile.curIsProperPair()){
						_discard.add(_bamFile.curName());
					}
				}
			}
		}
	}

	//close
	for(auto& s : out){
		s.close();
	}
};

}; // end namespace
