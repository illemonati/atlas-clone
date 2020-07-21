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

	//parse output fields
	_logfile->startIndent("Will print the following fields (parameter 'fields'):");
	std::set<std::string> fields;
	Parameters.fillParameterIntoSetWithDefault("fields", fields, ',', "depth,bases,qualities,alleles,mates,strands,likelihoods");

	_parseField(fields, "depth", _printDepth, "sequencing depth");
	_parseField(fields, "bases", _printBases, "pileup bases");
	_parseField(fields, "qualities", _printQualities, "pileup qualities");
	_parseField(fields, "alleles", _printAlleles, "allele counts");
	_parseField(fields, "mates", _printMates, "mate information");
	_parseField(fields, "strands", _printStrand, "strand information");
	_parseField(fields, "likelihoods", _printLikelihoods, "genotype likelihoods (likelihoods)");

	//check if unknown fields were given
	if(!fields.empty()){
		if(fields.size() == 1){
			throw "Unknown field '" + fields[0] + "'! Valid fields are 'depth', 'bases', 'qualities', 'alleles', 'mates' and 'strands'.";
		} else {
			bool first = true;
			std::string f;
			for(auto i : fields){
				if(first){
					first = false;
					f = "'";
				} else {
					f += ", '";
				}
				f += i + "'";
			}
			throw "Unknown fields: " + f + "! Valid fields are 'depth', 'bases', 'qualities', 'alleles', 'mates' and 'strands'.";
		}
	}

	//compile header
	std::vector<std::string> header = {"chromosome", "position"};
	if(_reference.hasReference()){ header.push_back("reference"); }
	if(_printDepth){
		header.push_back("depth");
		if(_reference.hasReference()){
			header.push_back("refDepth");
		}
	}
	if(_printBases){ header.push_back("bases"); }
	if(_printQualities){ header.push_back("qualities"); }
	if(_printAlleles){
		header.push_back("numA");
		header.push_back("numC");
		header.push_back("numG");
		header.push_back("numT");
		header.push_back("numRef");
		header.push_back("numNonRef");
		header.push_back("numAlleles");
	}
	if(_printMates){
		header.push_back("numFirstMate");
		header.push_back("numSeondMate");
	}
	if(_printStrand){
		header.push_back("numForwardStrand");
		header.push_back("numReverseStrand");
	}
	if(_printLikelihoods){
		_genoLik.addNames(header, _genoMap);
	}

	out.writeHeader(header);

	//print all sites, also those without data?
	if(Parameters.parameterExists("printAll")){
		printOnlySitesWithData = true;
		_logfile->list("Will all sites that pass filters, including those without data. (parameter 'printAll')");
	} else {
		printOnlySitesWithData = true;
		_logfile->list("Will print only sites with data. (use 'printAll' to print all)");
	}
};

void TPileup::_parseField(std::set<std::string> & fields, const std::string tag, bool & flag, const std::string explanation){
	if(fields.find(tag)){
		flag = true;
		_logfile->list(explanation + " (" + tag + ")");
			fields.erase(fields.find(tag));
	} else {
		flag = false;
	}
};

void TPileup::_handleWindow(){
	_logfile->listFlushTime("Writing pileup ...");

	uint32_t pos = 0;
	for(auto it = _window.begin(); it != _window.endPos(); ++it, ++pos){
		GenotypeLikelihoods::TSite& site = *it;

		out << _window._chrName;
		out << _window.posInRef(pos);

		//reference
		if(_reference.hasReference()){
			out << site._referenceBase;
		}

		//depth
		if(_printDepth){
			out << site.depth();
			if(_reference.hasReference()){
				out << site.refDepth();
			}
		}

		if(_printBases){
			out << site.getBases(_genoMap);
		}
		if(_printQualities){
			out << site.getQualities(_qualMap);
		}

		if(_printAlleles){
			site.countAlleles(_alleleCounts);
			out << _alleleCounts[A] << _alleleCounts[C] << _alleleCounts[G] << _alleleCounts[T];
			out << _alleleCounts.numAlleles();
		}

		if(_printMates){
			site.countMates(_counts);
			out << _counts[0] << _counts[1];
		}

		if(_printStrand){
			site.countFwdRev(_counts);
			out << _counts[0] << _counts[1];
		}

		if(_printLikelihoods){
			_genotypeLikelihoodCalculator.calculateGenotypeLikelihoods(site._bases, _genoLik);
			_genoLik.write(out);
		}

		out << std::endl;
	}

	_logfile->doneTime();
};

void TPileup::printPileup(){
	_traverseBAMWindows();
};

}; // end namespace
