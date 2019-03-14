/*
 * TSimulatorQuality.cpp
 *
 *  Created on: Oct 5, 2017
 *      Author: vivian
 */

#include "TSimulatorRead.h"

//----------------------------------
//TSimulatorSingleEndRead
//----------------------------------
TSimulatorSingleEndRead::TSimulatorSingleEndRead(std::string readGroupName, int readGroupNumber, int MaxPrintPhredInt, TRandomGenerator* RandomGenerator){
	type = "single-end";

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
	initializeAlignment(bamAlignment);
	readNamePrefix = "ATL:0:A:1:" + toString(readGroupNumber) + ":"; //"<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:"  Still need to add "<x-pos>:<y-pos>"
	readXPos = 1;
	readYPos = 1;

	bases = NULL;
	phredIntQualities = NULL;
};

TSimulatorSingleEndRead::~TSimulatorSingleEndRead(){
	if(readLengthInitialized){
		delete readLengthDist;
		delete[] bases;
		delete[] phredIntQualities;
	}
	if(qualityDistInitialized)
		delete qualityDist;
	if(qualityTransformInitialized)
		delete qualityTransform;
};

void TSimulatorSingleEndRead::initializeAlignment(BamTools::BamAlignment & alignment){
	alignment.AddTag("RG", "Z", _name);
	alignment.MapQuality = 50;
};

bool TSimulatorSingleEndRead::checkInitialization(){
	isInitialized = readLengthInitialized && qualityDistInitialized && qualityTransformInitialized;
	return isInitialized;
}

void TSimulatorSingleEndRead::setReadLengthDistribution(std::string s, TLog* logfile){
	size_t pos = s.find("(");
	std::string tmp;

	if(pos == std::string::npos)
		throw "Unable to understand read length distribution '" + s + "'!";

	//initialize appropriate function
	std::string type = s.substr(0, pos);
	s.erase(0, pos);

	if(type == "gamma")
		readLengthDist = new TSimulatorReadLengthGamma(s, randomGenerator, logfile);
	else if(type == "gammaMode"){
		readLengthDist = new TSimulatorReadLengthGammaMode(s, randomGenerator, logfile);

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

void TSimulatorSingleEndRead::setQualityDistribution(std::string s){
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

void TSimulatorSingleEndRead::setQualityTransformation(const std::string & type, const std::string & args, TLog* logfile){
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

void TSimulatorSingleEndRead::setPMD(const std::string & pmdStringCT, const std::string & pmdStringGA){
	pmdObject.initializeFunction(pmdStringCT, pmdCT);
	pmdObject.initializeFunction(pmdStringGA, pmdGA);

	hasPMD = pmdObject.hasDamage();
	checkInitialization();
};

void TSimulatorSingleEndRead::applyPMD(Base* _bases, BamTools::BamAlignment & alignment, int & fragmentLength){
	if(hasPMD){
		if(alignment.IsReverseStrand()){
			for(int p=0; p<alignment.Length; ++p){
				if(_bases[p] == C ){
					if(randomGenerator->getRand() < pmdObject.getProbGA(fragmentLength - alignment.Length + p))
						_bases[p] = T;
				} else if(_bases[p] == G){
					if(randomGenerator->getRand() < pmdObject.getProbGA(alignment.Length - p))
						_bases[p] = A;
				}
			}
		} else {
			for(int p=0; p<alignment.Length; ++p){
				if(_bases[p] == C ){
					if(randomGenerator->getRand() < pmdObject.getProbCT(p))
						_bases[p] = T;
				} else if(_bases[p] == G){
					if(randomGenerator->getRand() < pmdObject.getProbGA(fragmentLength - p - 1))
						_bases[p] = A;
				}
			}
		}
	}
};

void TSimulatorSingleEndRead::setContamination(double rate, TSimulatorReference* source){
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

void TSimulatorSingleEndRead::fillAlignmentDetails(BamTools::BamAlignment & alignment, const Base* theBases, const int* thePhredIntQualities){
	bamAlignment.CigarData.clear();
	bamAlignment.CigarData.push_back(BamTools::CigarOp('M', bamAlignment.Length));
	bamAlignment.QueryBases = "";
	bamAlignment.Qualities = "";

	//set read name
	++readXPos;
	if(readXPos == 65536){
		++readYPos;
		readXPos = 1;
	}
	bamAlignment.Name = readNamePrefix + toString(readXPos) + ":" + toString(readYPos);

	//copy bases and qualities
	for(int p=0; p<bamAlignment.Length; ++p){
		bamAlignment.QueryBases += genoMap.baseToChar[theBases[p]];
		bamAlignment.Qualities += (char) (std::min(thePhredIntQualities[p], maxPrintPhredInt) + 33);
	}
};

void TSimulatorSingleEndRead::simulate(Base* haplotype, const long & pos, TSimulatorBamFile & bamFile){
	//pick a fragment and read length
	int fragmentLength;
	readLengthDist->sample(bamAlignment.Length, fragmentLength);

	//fill bam alignment
	bamAlignment.Position = pos;
	bamAlignment.SetIsReverseStrand( randomGenerator->getRand() < 0.5);

	//copy bases
	if(isContaminated && randomGenerator->getRand() < contaminationRate)
		memcpy(bases, contaminationSource->getPointerToRef() + pos, bamAlignment.Length);
	else
		memcpy(bases, haplotype + pos, bamAlignment.Length);

	//apply PMD
	applyPMD(bases, bamAlignment, fragmentLength);

	//simulate qualities and errors
	qualityTransform->simulateQualitiesAndErrors(bases, phredIntQualities, bamAlignment.Length);

	//add to alignment
	fillAlignmentDetails(bamAlignment, bases, phredIntQualities);

	//save
	bamFile.saveAlignment(bamAlignment);
};

void TSimulatorSingleEndRead::printDetails(TLog* logfile){
	logfile->list(type + ".");
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

//----------------------------------
// TSimulatorPairedEndReads
//----------------------------------
TSimulatorPairedEndReads::TSimulatorPairedEndReads(std::string readGroupName, int readGroupNumber, int MaxPrintQual, TRandomGenerator* RandomGenerator)
:TSimulatorSingleEndRead(readGroupName, readGroupNumber, MaxPrintQual, RandomGenerator){
	type = "paired-end";
};

TSimulatorPairedEndReads::~TSimulatorPairedEndReads(){};


void TSimulatorPairedEndReads::simulate(Base* haplotype, const long & pos, TSimulatorBamFile & bamFile){
	//pick a fragment and read length
	int readLength; int fragmentLength;
	readLengthDist->sample(readLength, fragmentLength);

	//fill first mate
	bamAlignment.Position = pos;
	bamAlignment.SetIsReverseStrand(false);

	//fill second mate: reuse from idle pile, or create anew
	BamTools::BamAlignment* secondMate;
	if(bamAlignmentsSecondMate_idle.size() > 0){

	} else { //create one anew

	}




	//copy bases for both mates
	if(isContaminated && randomGenerator->getRand() < contaminationRate)
		memcpy(bases, contaminationSource->getPointerToRef() + pos, bamAlignment.Length);
	else
		memcpy(bases, haplotype + pos, bamAlignment.Length);

	//apply PMD
	applyPMD(bases, bamAlignment, fragmentLength);

	//simulate qualities and errors
	qualityTransform->simulateQualitiesAndErrors(bases, phredIntQualities, bamAlignment.Length);

	//add to alignment
	fillAlignmentDetails(bamAlignment, bases, phredIntQualities);

	//save
	bamFile.saveAlignment(bamAlignment);



};

void TSimulatorPairedEndReads::writeUnwrittenAlignments(const long & pos){

}


