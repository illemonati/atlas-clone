/*
 * TQualityDistribution.cpp
 *
 *  Created on: Jun 2, 2020
 *      Author: phaentu
 */

#include "TQualityDistribution.h"
#include "coretools/Main/TParameters.h"

namespace GenomeTasks{
using coretools::instances::logfile;
using coretools::instances::parameters;

TQualityTransformation::TQualityTransformation() {
	//check what we compare
	if(parameters().exists("RGInfo2")){
		BAM::RGInfo::TReadGroupInfo RGInfo2(_alnTraverser.bamFile().readGroups(), parameters().get("RGInfo2"));
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
}

void TQualityTransformation::_traverseAlignments() {
	for (; !_alnTraverser.endOfAlignments(); _alnTraverser.nextAlignment()) {
		const auto& alignment = _alnTraverser.alignment();
		if (_compareToOtherSeqErrors) {
			for (auto &b : alignment) {
				if (b.base != genometools::Base::N) {
					_transformations[alignment.readGroupId()].add(b.recalQuality.get(),
																  _otherSeqErrors.phredInt(b).get());
				}
			}
		} else {
			for (auto &b : alignment) {
				if (b.base != genometools::Base::N) {
					_transformations[alignment.readGroupId()].add(b.originalQuality.get(), b.recalQuality.get());
				}
			}
		}
	}
}

void TQualityTransformation::run(){
	//initialize transformations
	_transformations.resize(_alnTraverser.bamFile().numReadGroups());

	//traverseBAM
	_traverseAlignments();

	//write read group specific files
	logfile().startIndent("Writing quality transformation for each read group:");
	for (size_t rg=0; rg<_alnTraverser.bamFile().numReadGroups(); ++rg){
		std::string filename = _alnTraverser.outputName() + _alnTraverser.bamFile().readGroups().getName(rg) + "_qualityTransformation.txt";
		logfile().listFlush("Writing '" + filename + "' ...");
		_transformations[rg].writeAsMatrix(filename, _label1, _label2);
		logfile().done();
		logfile().conclude("R squared for read group " + _alnTraverser.bamFile().readGroups().getName(rg) + " is ", _transformations[rg].RSquared(), ".");
	}

	//write combined distribution
	coretools::TCountDistributionVector combined;
	for(auto& t : _transformations){
		combined.add(t);
	}

	std::string filename = _alnTraverser.outputName() + "_qualityTransformation.txt";
	logfile().listFlush("Writing quality transformation of total data to '" + filename + "' ...");
	combined.writeAsMatrix(filename, _label1, _label2);
	logfile().done();
	logfile().conclude("R squared for total data is ", combined.RSquared(), ".");

};

}; // end namespace
