/*
 * TBAmDownsampler.cpp
 *
 *  Created on: May 31, 2020
 *      Author: phaentu
 */

#include "TBamDownsampler.h"

namespace GenomeTasks{

//-----------------------------------------
// TBamSample
//-----------------------------------------
TBamSample::TBamSample(const double Prob, const std::string OutName){
	_prob = Prob;
	_outName = OutName;
};

void TBamSample::open(BAM::TBamFile & bamFile){
	_out.open(_outName, bamFile);
};

void TBamSample::close(TLog* logfile){
	_out.close(logfile);
};

void TBamSample::sample(BAM::TBamFile & bamFile, TRandomGenerator & randomGenerator){
	if(_discard.isInBlacklist(bamFile.curName())){
		_discard.remove(bamFile.curName());
	} else if(_keep.isInBlacklist(bamFile.curName())){
		bamFile.writeCurAlignment(_out);
		_keep.remove(bamFile.curName());
	} else if(randomGenerator.getRand() < _prob){
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

void TBamSample::downsampleRead(BAM::TAlignment & alignment, TRandomGenerator & randomGenerator, const TGenotypeMap & genoMap, const TQualityMap & qualMap){
	//parse again to get original bases and qualities
	alignment.parse(genoMap, qualMap);

	//downsample
	alignment.downsampleAlignment(_prob, randomGenerator);

	//write
	_out.writeAlignment(alignment, genoMap, qualMap);
};

//-----------------------------------------
// TBamDownsampler_base
//-----------------------------------------
TBamDownsampler_base::TBamDownsampler_base(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator):TGenome_basic(Parameters, Logfile, RandomGenerator){

};

void TBamDownsampler_base::_readVectorOfDownsamplingProbabilities(TParameters & Params){
	//read downsampling rates
	Parameters.fillParameterIntoProbabilityVector("prob", _probs, ',');

	//get unique names
	std::map <double, int> fracNames;
	for(size_t i=0; i<_probs.size(); ++i){
		std::map<double, int>::iterator it = fracNames.find(_probs[i]);
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
TBamDownsampler::TBamDownsampler(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator):TBamDownsampler_base(Parameters, Logfile, RandomGenerator){
	//read downsampling rates
	_readVectorOfDownsamplingProbabilities(Params);

	//report
	_logfile->list("Will accept reads with probabilities (parameter 'prob'): " + concatenateString(_probs, ", "));
	if(*_probs.begin() == 1.0) _logfile->warning("Probability of 1 will result in identical file!");

	//create downsampling objects
	for(size_t i=0; i<_probs.size(); ++i){
		std::string filename = _outputName + "_downsampled_" + _names[i] + ".bam";
		_bamSamples.emplace_back(_probs[i], filename);
	}

	//open bam files for writing
	for(auto& s : _bamSamples){
		s.open(_bamFile);
	}
};

void TBamDownsampler::downsample(){
	//traverse BAM and downsample
	_bamFile.startProgressReporting();
	while(_bamFile.readNextAlignment()){
		for(auto& s : _bamSamples){
			s.sample(_bamFile, *_randomGenerator);
		}
		_bamFile.printProgress();
	}
	_bamFile.printEndWithSummary();

	//close
	for(auto& s : _bamSamples){
		s.close(_logfile);
	}
};

//-----------------------------------------
// TBamReadDownsampler
//-----------------------------------------
TBamReadDownsampler::TBamReadDownsampler(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator):TBamDownsampler(Parameters, Logfile, RandomGenerator){

};

void TBamReadDownsampler::downsample(){
	BAM::TAlignment alignment;

	//traverse BAM and downsample
	_bamFile.startProgressReporting();
	while(_bamFile.readNextAlignment()){
		_bamFile.fill(alignment);
		for(auto& s : _bamSamples){
			s.downsampleRead(alignment, *_randomGenerator, _genoMap, _qualMap);
		}
		_bamFile.printProgress();
	}
	_bamFile.printEndWithSummary();

	//close
	for(auto& s : _bamSamples){
		s.close(_logfile);
	}
};

//-----------------------------------------
// TBamSeparator
//-----------------------------------------
TBamSeparator::TBamSeparator(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator):TBamDownsampler_base(Parameters, Logfile, RandomGenerator){
	//read downsampling rates
	_readVectorOfDownsamplingProbabilities(Params);

	//check that sum <= 1.0
	double sum = 0.0;
	for(auto& d : _probs){
		sum += d;
		_cumulProbs.push_back(sum);
	}
	_cumulProbs.push_back(1.0); //always add an extra at end to ease search

	if(sum > 1.0){
		throw "Separation probabilities must sum to <= 1.0, not " + toString(sum) + "!";
	}

	//report
	_logfile->list("Will separate reads with probabilities (parameter 'prob'): " + concatenateString(_probs, ", "));
	if(*_probs.begin() == 1.0) _logfile->warning("Probability of 1 will result in identical file!");
};

void TBamSeparator::separate(){
	//open bam files for writing
	std::vector<BAM::TOutputBamFile> out;
	for(auto& n : _names){
		std::string filename = _outputName + "_downsampled_" + n + ".bam";
		out.emplace_back(filename, _bamFile);
	}

	//preapre lists to keep track of mates
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
				double r = _randomGenerator->getRand();

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
		s.close(_logfile);
	}
};

}; // end namespace
