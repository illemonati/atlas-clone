/*
 * TSimulatorQuality.cpp
 *
 *  Created on: Oct 5, 2017
 *      Author: vivian
 */

#include "TSimulatorRead.h"
#include "TParameters.h"


//----------------------------------
//TSimulatorRead
//----------------------------------
TSimulatorRead::TSimulatorRead(std::string readGroupName, int MaxPrintPhredInt, TRandomGenerator* RandomGenerator){
	//set variables
	randomGenerator = RandomGenerator;
	maxPrintPhredInt = MaxPrintPhredInt;

	readLengthDist = NULL;
	readLengthInitialized = false;

	qualityDist = NULL;
	qualityDistInitialized = false;
	qualDistType = "";

	qualityTransform = NULL;
	qualityTransformInitialized = false;

	hasPMD = false;
	isInitialized = false;

	isContaminated = false;
	contaminationRate = 0.0;
	contaminationSource = NULL;

	//initialize bamAlignment
	_name = readGroupName;
	bamAlignment.AddTag("RG", "Z", _name);
	bamAlignment.MapQuality = 50;
	bamAlignment.Name = "*";
	bases = NULL;
	phredIntQualities = NULL;
	fragmentLength = 0;

	//tmp variables
	p = 0;
	previousBase = N;
	tmp_qual = 0;
};

bool TSimulatorRead::checkInitialization(){
	isInitialized = readLengthInitialized && qualityDistInitialized && qualityTransformInitialized;
	return isInitialized;
}

void TSimulatorRead::setReadLengthDistribution(std::string s){
	size_t pos = s.find("(");
	std::string tmp;

	if(pos == std::string::npos)
		throw "Unable to understand read length distribution '" + s + "'!";

	//initialize appropriate function
	std::string type = s.substr(0, pos);
	s.erase(0, pos);

	if(type == "gamma")
		readLengthDist = new TSimulatorReadLengthGamma(s, randomGenerator);
	else if(type == "gammaMode"){
		readLengthDist = new TSimulatorReadLengthGammaMode(s, randomGenerator);

	}
	else if(type == "fixed")
		readLengthDist = new TSimulatorReadLength(s, randomGenerator);
	else throw "Unknown read length distribution '" + type + "'!";

	//initialize bases and qualities
	bases = new Base[readLengthDist->max()];
	phredIntQualities = new int[readLengthDist->max()];

	//update status
	readLengthInitialized = true;
	checkInitialization();
};

void TSimulatorRead::setQualityDistribution(std::string s){
	size_t pos = s.find("(");
	std::string tmp;

	if(pos == std::string::npos)
		throw "Unable to understand read length distribution '" + s + "'!";

	//initialize appropriate function
	std::string type = s.substr(0, pos);
	s.erase(0, pos);
	if(type == "fixed")
		qualityDist = new TSimulatorQualityDist(s);
	else if(type == "normal")
		qualityDist = new TSimulatorQualityDistNormal(s, randomGenerator);
	else if(type == "binned")
		qualityDist = new TSimulatorQualityDistBinned(s, randomGenerator);
	else throw "Unknown read quality distribution '" + type + "'!";

	qualDistType = type;

	qualityDistInitialized = true;
	checkInitialization();
};

void TSimulatorRead::setQualityTransformation(const std::string & type, const std::string & args, TLog* logfile){
	if(!qualityDistInitialized)
		throw "Can not initialize quality transformation in TSimulatorRead: quality distribution not initialized!";

	if(type == "none")
		qualityTransform = new TSimulatorQualityTransformation(qualityDist, randomGenerator);
	else if(type == "recal")
		qualityTransform = new TSimulatorQualityTransformationRecal(args, readLengthDist->max(), qualityDist, randomGenerator);
	else if(type == "BQSR"){
		if(qualDistType != "normal") throw "Cannot apply BQSR transformation to any quality distribution besides 'normal'!";
		qualityTransform = new TSimulatorQualityTransformationBQSR(args, readLengthDist, logfile, qualityDist, randomGenerator);
	}
	else throw "Unknown quality transformation type '" + type + "'!";

	qualityTransformInitialized = true;
	checkInitialization();
};

void TSimulatorRead::setPMD(const std::string & pmdStringCT, const std::string & pmdStringGA){
	pmdObject.initializeFunction(pmdStringCT, pmdCT);
	pmdObject.initializeFunction(pmdStringGA, pmdGA);

	hasPMD = pmdObject.hasDamage();
	checkInitialization();
};

void TSimulatorRead::applyPMD(Base* _bases, int & len, int & fragmentLength){
	if(hasPMD){
		for(p=0; p<len; ++p){
			if(_bases[p] == C ){
				if(randomGenerator->getRand() < pmdObject.getProbCT(p)){
					_bases[p] = T;
				}
			} else if(_bases[p] == G){
				if(randomGenerator->getRand() < pmdObject.getProbGA(fragmentLength - p - 1)){
					_bases[p] = A;
				}
			}
		}
	}
};

void TSimulatorRead::setContamination(double rate, TSimulatorReference* source){
	isContaminated = true;
	contaminationRate = rate;
	contaminationSource = source;

	//check
	if(contaminationRate < 0.0)
		throw "Contamination rate must be >= 0.0!";
	if(contaminationRate > 1.0)
		throw "Contamination rate must be <= 0.0!";
	if(contaminationRate == 0.0)
		isContaminated = false;
};

void TSimulatorRead::simulate(Base* haplotype, const long & pos, TSimulatorBamFile & bamFile){
	//pick a fragment and read length
	readLengthDist->sample(bamAlignment.Length, fragmentLength);

	//fill bam alignment
	bamAlignment.CigarData.clear();
	bamAlignment.CigarData.push_back(BamTools::CigarOp('M', bamAlignment.Length));
	bamAlignment.Position = pos;
	bamAlignment.QueryBases = "";
	bamAlignment.Qualities = "";

	//copy bases
	if(isContaminated && randomGenerator->getRand() < contaminationRate)
		memcpy(contaminationSource->getPointerToRef(), contaminationSource->getPointerToRef() + pos, bamAlignment.Length);
	else
		memcpy(bases, haplotype + pos, bamAlignment.Length);

	//apply PMD
	applyPMD(bases, bamAlignment.Length, fragmentLength);

	//simulate qualities and errors
	qualityTransform->simulateQualitiesAndErrors(bases, phredIntQualities, bamAlignment.Length);

	//add to alignment
	bamAlignment.QueryBases = "";
	bamAlignment.Qualities = "";
	for(p=0; p<bamAlignment.Length; ++p){
		bamAlignment.QueryBases += genoMap.baseToChar[bases[p]];
		bamAlignment.Qualities += (char) (std::min(phredIntQualities[p], maxPrintPhredInt) + 33);
	}

	//save
	bamFile.saveAlignment(bamAlignment);
};

void TSimulatorRead::printDetails(TLog* logfile){
	if(readLengthInitialized)
		readLengthDist->printDetails(logfile);
	else throw "Read length distribution not initialized!";

	if(qualityDistInitialized)
		qualityDist->printDetails(logfile);
	else throw "Read quality distribution not initialized!";

	if(qualityTransformInitialized)
		qualityTransform->printDetails(logfile);
	else throw "Quality transformation not initialized!";

	if(hasPMD){
		logfile->list("C->T PMD as " + pmdObject.getFunctionString(pmdCT));
		logfile->list("G->A PMD as " + pmdObject.getFunctionString(pmdGA));
	} else {
		logfile->list("No PMD.");
	}

	if(isContaminated)
		logfile->list("Contaminated with rate " + toString(contaminationRate) + ".");
	else
		logfile->list("Read group is not contaminated.");
};




