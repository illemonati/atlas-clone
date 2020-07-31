/*
 * TChromosomes.cpp
 *
 *  Created on: Apr 25, 2019
 *      Author: linkv
 */

#include "TChromosomes.h"
#include "TLog.h"

namespace BAM{

//---------------------------------------------------------
// TChromosome
//---------------------------------------------------------
TChromosome::TChromosome(const uint32_t RefID, const std::string Name, const uint32_t Length){
	name = Name;
	length = Length;
	ploidy = 2; //default: diploid
	inUse = true;

	//set TGenomePosition
	chrStart.move(RefID, 0);
	chrEnd.move(RefID, length); //end is not included
};

std::string TChromosome::compileSamHeader() const{
	std::string header = "@SQ\tSN:" + name + "\tLN:" + toString(length);

	if(!alternateLocus.empty()){
		header += "\tAH:" + alternateLocus;
	}
	if(!alternativeNames.empty()){
		header += "\tAN:" + alternativeNames;
	}
	if(!assemblyIdentifier.empty()){
		header += "\tAS:" + assemblyIdentifier;
	}
	if(!description.empty()){
		header += "\tDS:" + description;
	}
	if(!MD5.empty()){
		header += "\tM5:" + MD5;
	}
	if(!species.empty()){
		header += "\tSP:" + species;
	}
	if(!topology.empty()){
		header += "\tTP:" + topology;
	}
	if(!uri.empty()){
		header += "\tUR:" + uri;
	}

	return header + "\n";
};

//---------------------------------------------------------
// TChromosomes
//---------------------------------------------------------
void TChromosomes::clear(){
	_chromosomes.clear();
	_curChr = _chromosomes.end();
};

void TChromosomes::appendChromosome(const std::string name, const uint32_t length){
	_chromosomes.emplace_back(_chromosomes.size(), name, length);
};

const TChromosome& TChromosomes::_find(const std::string chrName) const{
	for(auto& c : _chromosomes){
		if(c.name == chrName)
			return c;
	}
	throw "Chromosome '" + chrName + "' not found in BAM header!";
};

void TChromosomes::_setAllNotInUse(){
	for(auto& c : _chromosomes){
		c.inUse = false;
	}
};

void TChromosomes::limitAndSetPloidy(TParameters & params, TLog* logfile){
	limitChr(params, logfile);
	setPloidy(params, logfile);
};

void TChromosomes::limitChr(TParameters & params, TLog* logfile){
	if(params.parameterExists("chr")){
		logfile->startIndent("Will limit analysis to the following chromosomes (parameter 'chr':");
		std::vector<std::string> vec;
		fillVectorFromString(params.getParameterString("chr"), vec, ',');
		_useSpecifiedChr(vec);
		writeUsedChromosomes(logfile);
	} else {
		if(params.parameterExists("limitChr")){
			std::string limitName = params.getParameterString("limitChr");
			logfile->list("Will limit analysis to all chromosomes up to and including " + limitName + ". (parameter 'limitChr')");
			_useUpTo(limitName);
		}
	}
};

void TChromosomes::_useSpecifiedChr(const std::vector<std::string> & chrNames){
	//set all chromosomes to false
	_setAllNotInUse();

	//set matching to to true
	for(auto& it : chrNames){
		//find chromosome
		auto c = _find(it);
		c.inUse = true;
	}
};

void TChromosomes::_useUpTo(const std::string & name){
	_find(name);

	//set in use
	bool use = true;
	for(auto& it : _chromosomes){
		if(it.name == name){
			use = false;
		}
		it.inUse = use;
	}
};

void TChromosomes::writeUsedChromosomes(TLog* logfile){
	for(auto& it : _chromosomes){
		if(it.inUse){
			logfile->list(it.name);
		}
	}

	logfile->endIndent();
};

void TChromosomes::setPloidy(TParameters & params, TLog* logfile){
	if(params.parameterExists("ploidy")){
		_specifyPloidy(params.getParameterString("ploidy"), logfile);
	} else if(params.parameterExists("haploid")){
		std::vector<std::string> vec;
		fillVectorFromString(params.getParameterString("haploid"), vec, ',');
		if(vec.size() == 1 && vec[0] == "all"){
			logfile->list("Assuming all chromosomes are haploid. (parameter 'haploid')");
			for(auto& c : _chromosomes){
				c.ploidy = 1;
			}
		} else {
			_setToHaploid(vec, logfile);
		}
	} else {
		logfile->list("Will assume that all chromosomes are diploid (use 'ploidy' or 'haploid' to change).");
	}
};

void TChromosomes::_specifyPloidy(const std::string & ploidyFileName, TLog* logfile){
	std::ifstream ploidyFile(ploidyFileName.c_str());
	logfile->list("Reading ploidy specification per chromosome from file '" + ploidyFileName + "'.");
	if(!ploidyFile)
		throw "Failed to open file '" + ploidyFileName + "'!";

	logfile->startIndent("Setting ploidy for following chromosomes:");
	uint32_t numChrInFile = 0;
	while(ploidyFile.good() && !ploidyFile.eof()){
		//read line
		std::string line;
		std::getline(ploidyFile, line);
		std::vector<std::string> vec;
		fillVectorFromStringWhiteSpace(line, vec, true);

		//skip empty lines
		if(vec.size() > 0){
			auto c = _find(vec[0]);
			c.ploidy = convertString<int>(vec[1]);

			//report
			if(c.ploidy == 1){
				logfile->list(c.name + ": haploid");
			} else if(c.ploidy == 2){
				logfile->list(c.name + ": diploid");
			} else {
				throw "Unsupported ploidy '" + toString(c.ploidy) + "'! Currently only support haploid (1) and diploid (2).";
			}
			++numChrInFile;
		}
	}
	ploidyFile.close();
	if(numChrInFile < _chromosomes.size()){
		logfile->list("Remaining " + toString(_chromosomes.size() - numChrInFile) + " chromosomes: diploid");
	}
	logfile->endIndent();
};


void TChromosomes::_setToHaploid(const std::vector<std::string> & chrNames, TLog* logfile){
	logfile->startIndent("Setting the following chromosomes to be haploid:");
	for(auto& it : chrNames){
		auto c = _find(it);
		c.ploidy = 1;
		logfile->list(c.name);
	}

	if(chrNames.size() < _chromosomes.size()){
		logfile->list("Will assume that the remaining " + toString(_chromosomes.size() - chrNames.size()) + " chromosomes are diploid.");
	}

	logfile->endIndent();
};

//getters
uint32_t TChromosomes::referenceLength() const{
	long totLength = 0;
	for(auto& c : _chromosomes){
		if(c.inUse)
			totLength += c.length;
	}
	return totLength;
};

bool TChromosomes::exists(const std::string chrName) const{
	for(auto& c : _chromosomes){
		if(c.name == chrName)
			return true;
	}
	return false;
};

const TChromosome& TChromosomes::getChromosome(const std::string chrName) const{
	return _find(chrName);
};

uint32_t TChromosomes::refID(const std::string chrName) const{
	return _find(chrName).refID();
};

uint32_t TChromosomes::length(const uint32_t & RefID) const{
	return _chromosomes[RefID].length;
};

std::string TChromosomes::name(const uint32_t & RefID) const{
	return _chromosomes[RefID].name;
};

bool TChromosomes::inUse(const uint32_t & RefID) const{
	return _chromosomes[RefID].inUse;
};

uint8_t TChromosomes::ploidy(const uint32_t & RefID) const{
	return _chromosomes[RefID].ploidy;
};

const TGenomePosition& TChromosomes::chrStart(const uint32_t & RefID) const{
	return _chromosomes[RefID].chrStart;
};

const TGenomePosition& TChromosomes::chrEnd(const uint32_t & RefID) const{
	return _chromosomes[RefID].chrEnd;
};

std::string TChromosomes::compileSamHeader() const{
	std::string header;
	for(auto& c : _chromosomes){
		header += c.compileSamHeader();
	}
	return header;
};

}; //end namespace
