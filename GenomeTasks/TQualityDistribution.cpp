/*
 * TQualityDistribution.cpp
 *
 *  Created on: Jun 2, 2020
 *      Author: phaentu
 */

#include "TQualityDistribution.h"

#include <stdint.h>
#include <algorithm>
#include <exception>
#include <memory>

#include "genometools/GenotypeTypes.h"
#include "genometools/PhredProbabilityTypes.h"
#include "TAlignment.h"
#include "TBamFile.h"
#include "coretools/Files/TFile.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "TReadGroups.h"
#include "TSequencedBase.h"

namespace GenomeTasks{
using coretools::instances::logfile;
using coretools::instances::parameters;

//-----------------------------------
// TQualityDistribution
//-----------------------------------
void TQualityDistribution::_handleAlignment(){
	for(auto& b : _alignment){
		if(b.base != genometools::Base::N){
			_qualDist.add(b.readGroupID, b.recalibratedQualityAsPhredInt.get());
		}
	}
};

void TQualityDistribution::compileQualityDistribution(){
	//initialize counts
	_qualDist.clear();
	_qualDist.resize(_bamFile.numReadGroups());

	//traverseBAM
	_traverseBAMPassedQC();

	//print distribution
	std::string filename = _outputName + "_qualityDistribution.h";
	logfile().listFlush("Writing quality distribution to '" + filename + "' ...");
	coretools::TOutputFile out(filename, {"readGroup", "quality", "counts"});

	//get read group names
	std::vector<std::string> readGroupNames;
	_bamFile.readGroups().fillVectorWithNames(readGroupNames);

	//write combined
	_qualDist.writeCombined(out, "allReadGroups");
	_qualDist.write(out, readGroupNames);
	logfile().done();
};

//-----------------------------------
// TQualityTransformation
//-----------------------------------
TQualityTransformation::TQualityTransformation():TGenome_parsed(){
	//check what we compare
	if(parameters().parameterExists("RGInfo2")){
		BAM::RGInfo::TReadGroupInfo RGInfo2(_bamFile.readGroupsMutable(), parameters().getParameter("RGInfo2"));
		_otherSeqErrors.initialize(RGInfo2);

		_compareToOtherSeqErrors = true;
		_label1 = "recalibratedQuality";
		_label2 = "recalbratedQuality2";
	} else {
		logfile().startIndent("Comparing original to recalibrated qualities.");
		_compareToOtherSeqErrors = false;
		_label1 = "originalQuality";
		_label2 = "recalbratedQuality";
	}
};

void TQualityTransformation::_handleAlignment(){
	if(_compareToOtherSeqErrors){
		for(auto& b : _alignment){
			if(b.base != genometools::Base::N){
				_transformations[_alignment.readGroupId()].add(b.recalibratedQualityAsPhredInt.get(), _otherSeqErrors.phredInt(b).get());
			}
		}
	} else {
		for(auto& b : _alignment){
			if(b.base != genometools::Base::N){
				_transformations[_alignment.readGroupId()].add(b.originalQuality_phredInt.get(), b.recalibratedQualityAsPhredInt.get());
			}
		}
	}
};

void TQualityTransformation::run(){
	//initialize transformations
	_transformations.resize(_bamFile.numReadGroups());

	//traverseBAM
	_traverseBAMPassedQC();

	//write read group specific files
	logfile().startIndent("Writing quality transformation for each read group:");
	for (size_t rg=0; rg<_bamFile.numReadGroups(); ++rg){
		std::string filename = _outputName + _bamFile.readGroups().getName(rg) + "_qualityTransformation.txt";
		logfile().listFlush("Writing '" + filename + "' ...");
		_transformations[rg].writeAsMatrix(filename, _label1, _label2);
		logfile().done();
		logfile().conclude("R squared for read group " + _bamFile.readGroups().getName(rg) + " is ", _transformations[rg].RSquared(), ".");
	}

	//write combined distribution
	coretools::TCountDistributionVector combined;
	for(auto& t : _transformations){
		combined.add(t);
	}

	std::string filename = _outputName + "_qualityTransformation.txt";
	logfile().listFlush("Writing quality transformation of total data to '" + filename + "' ...");
	combined.writeAsMatrix(filename, _label1, _label2);
	logfile().done();
	logfile().conclude("R squared for total data is ", combined.RSquared(), ".");

};

}; // end namespace
