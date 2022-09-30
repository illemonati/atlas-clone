/*
 * TPileup.cpp
 *
 *  Created on: Jun 4, 2020
 *      Author: wegmannd
 */

#include "TPileup.h"

#include <stdint.h>
#include <algorithm>
#include <ostream>
#include <vector>

#include "GenotypeTypes.h"
#include "TFastaBuffer.h"
#include "TGenotypeData.h"
#include "TGenotypeLikelihoodCalculator.h"
#include "TLog.h"
#include "TParameters.h"
#include "TSite.h"
#include "TWindow.h"

namespace GenomeTasks{

using coretools::instances::logfile;
using coretools::instances::parameters;

namespace impl {

void addNames(std::vector<std::string> &vec) {
	using GT = genometools::Genotype;
	for (auto g = GT::min; g < GT::max; ++g) vec.push_back("P(D|" + genometools::toString(g) + ")");
};

void write(const GenotypeLikelihoods::TGenotypeLikelihoods& lhs, coretools::TOutputFile & out) {
	for (const auto l: lhs)
		out << l;
};

bool parseField(std::set<std::string> &fields, const std::string &tag, const std::string &explanation) {
	if (fields.find(tag) != fields.end()) {
		logfile().list(explanation + " (" + tag + ")");
		fields.erase(fields.find(tag));
		return true;
	}
	return false;
}

} // namespace

//---------------------------------
// TPileup
//---------------------------------
TPileup::TPileup():TGenome_windows(){
	//open output file
	const std::string filename = _outputName + "_pileup.txt.gz";
	logfile().list("Writing pileup to file '" + filename + "'.");
	_out.open(filename);

	//parse output fields
	logfile().startIndent("Will print the following pileup fields (parameter 'fields'):");
	std::set<std::string> fields;
	parameters().fillParameterIntoContainerWithDefault("fields", fields, ',', {"depth", "bases", "qualities", "alleles", "mates", "strands", "likelihoods"});

	_printSettings.set<Depth>(impl::parseField(fields, "depth", "Sequencing depth"));
	_printSettings.set<Bases>(impl::parseField(fields, "bases", "Pileup bases"));
	_printSettings.set<Qualities>(impl::parseField(fields, "qualities", "Pileup qualities"));
	_printSettings.set<Alleles>(impl::parseField(fields, "alleles", "Allele counts"));
	_printSettings.set<Mates>(impl::parseField(fields, "mates", "Mate information"));
	_printSettings.set<Strand>(impl::parseField(fields, "strands", "Strand information"));
	_printSettings.set<Likelihoods>(impl::parseField(fields, "likelihoods", "Genotype likelihoods"));
	logfile().endIndent();

	//check if unknown fields were given
	if(!fields.empty()){
		if(fields.size() == 1){
			throw "Unknown field '" + *fields.begin() + "'! Valid fields are 'depth', 'bases', 'qualities', 'alleles', 'mates' and 'strands'.";
		} else {
			std::string f;
			for(auto i : fields){
				f += '\'' + i + "', ";
			}
			throw "Unknown fields: " + f.substr(0, f.size() - 2) + "! Valid fields are 'depth', 'bases', 'qualities', 'alleles', 'mates',  'strands' and 'likelihoods'.";
		}
	}

	//compile header
	std::vector<std::string> header = {"chromosome", "position"};
	if(_reference){ header.push_back("reference"); }
	if(_printSettings.get<Depth>()){
		header.push_back("depth");
		if(_reference){
			header.push_back("refDepth");
		}
	}
	if(_printSettings.get<Bases>()){ header.push_back("bases"); }
	if(_printSettings.get<Qualities>()){ header.push_back("qualities"); }
	if(_printSettings.get<Alleles>()){
		header.push_back("numA");
		header.push_back("numC");
		header.push_back("numG");
		header.push_back("numT");
		if(_reference){
			header.push_back("numRef");
			header.push_back("numNonRef");
		}
		header.push_back("numAlleles");
	}
	if(_printSettings.get<Mates>()){
		header.push_back("numFirstMate");
		header.push_back("numSeondMate");
	}
	if(_printSettings.get<Strand>()){
		header.push_back("numForwardStrand");
		header.push_back("numReverseStrand");
	}
	if(_printSettings.get<Likelihoods>()){
		impl::addNames(header);
	}

	_out.writeHeader(header);

	//print all sites, also those without data?
	if(parameters().parameterExists("printAll")){
		_printSettings.set<OnlySitesWithData>(false);
		logfile().list("Will print all sites that pass filters, including those without data. (parameter 'printAll')");
	} else {
		_printSettings.set<OnlySitesWithData>(true);
		logfile().list("Will print only sites with data. (use 'printAll' to print all)");
	}
};


void TPileup::_handleWindow(){
	using genometools::Base;
	logfile().listFlushTime("Writing pileup ...");

	for (size_t pos = 0; pos < _window.size(); ++pos) {
		const auto &site = _window[pos];
		if (_printSettings.get<OnlySitesWithData>() && site.empty()) continue;
		_out << _window.chrName();
		_out << _window.positionOnChr(pos) + 1; // positions are zero-based internally

		if(_reference){
			_out << site.refBase;
		}
		if(_printSettings.get<Depth>()){
			_out << site.depth();
			if(_reference){
				_out << site.refDepth();
			}
		}
		if(_printSettings.get<Bases>()){
			_out << site.getBases();
		}
		if(_printSettings.get<Qualities>()){
			_out << site.getQualities();
		}
		if(_printSettings.get<Alleles>()){
			GenotypeLikelihoods::TBaseCounts alleleCounts;
			site.countAlleles(alleleCounts);
			_out << alleleCounts[Base::A] << alleleCounts[Base::C] << alleleCounts[Base::G] << alleleCounts[Base::T];
			if(_reference){
				_out << alleleCounts[site.refBase] << alleleCounts.size() - alleleCounts[site.refBase];
			}
			_out << (int) coretools::numNonZero(alleleCounts);
		}
		if(_printSettings.get<Mates>()){
			std::array<int, 2> mateCounts;
			site.countMates(mateCounts);
			_out << mateCounts[0] << mateCounts[1];
		}
		if(_printSettings.get<Strand>()){
			std::array<int, 2> strandCounts;
			site.countFwdRev(strandCounts);
			_out << strandCounts[0] << strandCounts[1];
		}
		if(_printSettings.get<Likelihoods>()){
			_genoLik = _genotypeLikelihoodCalculator.calculateGenotypeLikelihoods(site);
			impl::write(_genoLik, _out);
		}
		_out.endln();
	}

	logfile().doneTime();
};

void TPileup::printPileup(){
	_traverseBAMWindows();
};

}; // end namespace
