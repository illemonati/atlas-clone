/*
 * TPileup.cpp
 *
 *  Created on: Jun 4, 2020
 *      Author: wegmannd
 */

#include "TPileup.h"

namespace GenomeTasks{

//---------------------------------
// TPileup
//---------------------------------
TPileup::TPileup(TParameters & Parameters, TLog* Logfile, TRandomGenerator* RandomGenerator):TGenome_windows(Parameters, Logfile, RandomGenerator){
	//open output file
	std::string filename = _outputName + "_pileup.txt.gz";
	_logfile->list("Writing pileup to file '" + filename + "'.");
	out.open(filename);

	//add header
	std::vector<std::string> header = {"chromosome", "position", "reference", "depth", "refDepth", "bases", "numA", "numC", "numG", "numT", "numAlleles", "numFirstMate", "numSecondMate", "numFwd", "numRev"};
	_genoLik.addNames(header, _genoMap);
	out.writeHeader(header);

	if(Parameters.parameterExists("printAll")){
		printOnlySitesWithData = true;
		_logfile->list("Will all sites, including those without data. (parameter 'printAll')");
	} else {
		printOnlySitesWithData = true;
		_logfile->list("Will print only sites with data. (use 'printAll' to print all)");
	}
};

void TPileup::_handleWindow(){
	_logfile->listFlushTime("Writing pileup ...");

	uint32_t pos = 0;
	for(auto it = _window.begin(); it != _window.endPos(); ++it, ++pos){
		TSite& site = *it;

		out << _window._chrName;
		out << _window.posInRef(pos);

		out << site.referenceBase;
		out << site.depth();
		out << site.refDepth();
		out << site.getBases(_genoMap);

		site.countAlleles(_alleleCounts);
		out << _alleleCounts[A] << _alleleCounts[C] << _alleleCounts[G] << _alleleCounts[T];
		out << _alleleCounts.numAlleles();

		site.countMates(_counts);
		out << _counts[0] << _counts[1];

		site.countFwdRev(_counts);
		out << _counts[0] << _counts[1];

		_genotypeLikelihoodCalculator.calculateGenotypeLikelihoods(site.bases, _genoLik);
		_genoLik.write(out);
		out << std::endl;
	}

	_logfile->doneTime();
};

void TPileup::printPileup(){
	_traverseBAMWindows();
};

}; // end namespace
