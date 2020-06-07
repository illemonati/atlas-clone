/*
 * TQualityDistribution.cpp
 *
 *  Created on: Jun 2, 2020
 *      Author: phaentu
 */

#include "TQualityDistribution.h"

namespace GenomeTasks{

//-----------------------------------
// TQualityDistribution
//-----------------------------------
TQualityDistribution::TQualityDistribution(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator):TGenome_parsed(Parameters, Logfile, RandomGenerator){

};

void TQualityDistribution::_handleAlignment(){
	for(auto& b : _alignment){
		_qualDist.add(b.readGroupID, b.recalibratedQualityAsPhredInt);
	}
};

void TQualityDistribution::compileQualityDistribution(){
	//initialize counts
	_qualDist.clear();
	_qualDist.resize(_bamFile.readGroups.size());

	//traverseBAM
	_traverseBAMPassedQC();

	//print distribution
	std::string filename = _outputName + "_qualityDistribution.h";
	_logfile->listFlush("Writing quality distribution to '" + filename + "' ...");
	TOutputFile out(filename, {"readGroup", "quality", "counts"});

	//get read group names
	std::vector<std::string> readGroupNames;
	_bamFile.readGroups.fillVectorWithNames(readGroupNames);

	//write combined
	_qualDist.writeCombined(out, "allReadGroups");
	_qualDist.write(out, readGroupNames);
	_logfile->done();
};

//-----------------------------------
// TQualityTransformation
//-----------------------------------
TQualityTransformation::TQualityTransformation(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator):TGenome_parsed(Parameters, Logfile, RandomGenerator){
	//check what we compare
	if(Parameters.parameterExists("recal2")){
		std::string recalstring = Parameters.getParameterString("recal2");
		_logfile->startIndent("Comparing recalibrated qualities to those recalibrated with altermative parameters:");
		_otherSeqErrors.createModels(recalstring, _bamFile.readGroups, _logfile);

		_compareToOtherSeqErrors = true;
		_label1 = "recalibratedQuality";
		_label2 = "recalbratedQuality2";
	} else {
		_logfile->startIndent("Comparing original to recalibrated qualities.");
		_compareToOtherSeqErrors = false;
		_label1 = "originalQuality";
		_label2 = "recalbratedQuality";
	}
};

void TQualityTransformation::_handleAlignment(){
	if(_compareToOtherSeqErrors){
		for(auto& b : _alignment){
			if(b.base != N){
				_transformations[_alignment.readGroupId()].add(b.recalibratedQualityAsPhredInt, _otherSeqErrors.getPhredInt(b));
			}
		}
	} else {
		for(auto& b : _alignment){
			if(b.base != N){
				_transformations[_alignment.readGroupId()].add(b.originalQuality_phredInt, b.recalibratedQualityAsPhredInt);
			}
		}
	}
};

void TQualityTransformation::compileQualityTransformation(){
	//initialize transformations
	_transformations.resize(_bamFile.readGroups.size());

	//traverseBAM
	_traverseBAMPassedQC();

	//write read group specific files
	_logfile->startIndent("Writing quality transformation for each read group:");
	for(uint16_t rg=0; rg<_bamFile.readGroups.size(); ++rg){
		std::string filename = _outputName + _bamFile.readGroups.getName(rg) + "_qualityTransformation.txt";
		_logfile->listFlush("Writing '" + filename + "' ...");
		_transformations[rg].writeAsMatrix(filename, _label1, _label2);
		_logfile->done();
		_logfile->conclude("R squared for read group " + _bamFile.readGroups.getName(rg) + " is " + toString(_transformations[rg].RSquared()) + ".");
	}

	//write combined distribution
	TCountDistributionVector combined;
	for(auto& t : _transformations){
		combined.add(t);
	}

	std::string filename = _outputName + "_qualityTransformation.txt";
	_logfile->listFlush("Writing quality transformation of total data to '" + filename + "' ...");
	combined.writeAsMatrix(filename, _label1, _label2);
	_logfile->done();
	_logfile->conclude("R squared for total data is " + toString(combined.RSquared()) + ".");

};

}; // end namespace
