/*
 * TAtlasTestRecalibration.cpp
 *
 *  Created on: Jan 22, 2018
 *      Author: linkv
 */

#include "TAtlasTestRecalibration.h"


//------------------------------------------
//TAtlasTest_recalSimulation
//------------------------------------------
TAtlasTest_recalSimulation::TAtlasTest_recalSimulation(TParameters & params, TLog* logfile):TAtlasTest(params, logfile){
	_name = "recalSimulation";
	filenameTag = _testingPrefix + _name;
	bamFileName = filenameTag + ".bam";
	meanQual = params.getParameterIntWithDefault("recal_meanQual", 25);
	sdphredInt = params.getParameterDoubleWithDefault("recal_sdQual", 10);
	minPhredInt = params.getParameterIntWithDefault("recal_minQual", 0);
	maxPhredInt = params.getParameterIntWithDefault("recal_maxQual", 42);
	qualityDist = params.getParameterStringWithDefault("recal_qualityDist", "normal(" + toString(meanQual) + "," + toString(sdphredInt) + ")[" + toString(minPhredInt) + "," + toString(maxPhredInt) + "]");
//	recalParamString = params.getParameterStringWithDefault("recal_recalParams", "2,0,0.1,0.001,1{20}");
	recalParamString = params.getParameterStringWithDefault("recal_recalParams","0.908163,0.0022877,-0.0160425,0.000170256,0.120439,1.50259,1.55807,0.607032,0.775844,1.1983,3.52317,-0.0538213,0.392298,1.07254,1.41819,-0.387901,0.949369,1.17807,1.3996,0.0631075,0.834644,1.08996,2.29066,-0.102391");
	fillVectorFromStringAnySkipEmpty(recalParamString, tmpVec, ",");
	repeatIndexes(tmpVec, trueParams);
	recalParamsFileName = filenameTag + "_true_recalibrationEM.txt";
	poolRGFileName = filenameTag + "_poolThese.txt";
}

bool TAtlasTest_recalSimulation::run(){
	//1) Write recal params to file
	//-----------------------------
	std::vector<std::string> paramVector;
	fillVectorFromStringAny(recalParamString, paramVector, ",");
	outRecalParams.open(recalParamsFileName.c_str());
	if(!outRecalParams) throw "Failed to open file '" + recalParamsFileName + "'!";

	outRecalParams << "readGroup\tmate\tmodel\tquality\tposition\tcontext\n";
	for(int i=0; i<3; ++i){
		std::string RGName = "RG" + toString(i);
		outRecalParams << RGName << "\tfirst\tqualFuncPosFuncContext\t";
		//quality
		outRecalParams << trueParams[0] << "," << trueParams[1] << "\t";
		//position
		outRecalParams << trueParams[2] << "," << trueParams[3] << "\t";
		//context
		outRecalParams << trueParams[4];
		for(unsigned int p=5; p<paramVector.size(); ++p)
			outRecalParams << "," << paramVector[p] ;
		//add likelihood
		outRecalParams << "\n";
	}
	outRecalParams.close();

	//1) Run ATLAS to simulate BAM file
	//-----------------------------
	//TODO: find minimal data necessary to run test in order to speed up

	_testParams.addParameter("out", filenameTag);
	_testParams.addParameter("chrLength", "2000000");
	_testParams.addParameter("ploidy", "1");
	_testParams.addParameter("depth", "1");
	_testParams.addParameter("qualityDist", qualityDist);
	_testParams.addParameter("recal", recalParamsFileName);
//	_testParams.addParameter("recal", "recal[" + recalParamString + "]");
//	_testParams.addParameter("readLength", "gamma(" + toString(alpha) + "," + toString(beta)+ ")[" + toString(minReadLen) + "," + toString(maxReadLen));
	_testParams.addParameter("readLength", "single:fixed(70)");

	if(!runTGenomeFromInputfile("simulate"))
		return false;

	logfile->newLine();

	//2) Run recal
	//-----------------------------
	//open pool read group file
	outRecalPool.open(poolRGFileName.c_str());
	if(!outRecalPool) throw "Failed to open file '" + poolRGFileName + "'!";
	outRecalPool << "RG0 RG1\n";
	outRecalPool.close();

	_testParams.clear();
	_testParams.addParameter("bam", bamFileName);
	_testParams.addParameter("poolReadGroups", poolRGFileName);

	if(!runTGenomeFromInputfile("recal"))
		return false;

	//3) check if results are OK
	//--------------------------
	if(checkRecalFile() == true)
		return true;
	else return false;
};

bool TAtlasTest_recalSimulation::checkRecalFile(){
	logfile->startIndent("Checking recal file:");

	//open quality file
	std::string filename = filenameTag + "_recalibrationEM.txt";
	logfile->listFlush("Opening file '" + filename + "' for reading ...");
	std::ifstream in(filename.c_str());
	if(!in)
		throw "Failed to open file '" + filename + "'!";
	logfile->done();

	//some variables
	std::string tmp, tmp2;
	std::vector<double> estimatedParams, estimatedParams2;

	//skip header
	getline(in, tmp);

	//read estimated params for RG0 and RG1
	getline(in, tmp);
	getline(in, tmp2);
	fillVectorFromStringAny(tmp, estimatedParams, "\t");
	fillVectorFromStringAny(tmp2, estimatedParams2, "\t");

	//parse true params
	logfile->startIndent("Checking parameter values for pooled read groups RG0 and RG1");
	for(unsigned int i=3; i<estimatedParams.size(); ++i){ //first three are model information
		if(estimatedParams[i] != estimatedParams2[i]){
			logfile->newLine();
			logfile->conclude("esimated value for parameter number " + toString(i) + " in RG0: " + toString(estimatedParams[i]) + " is not the same as in RG1: " + toString(estimatedParams2[i]));
			return false;
		}
		if(estimatedParams[i] != trueParams[i]){
			logfile->newLine();
			logfile->conclude("esimated value for parameter number " + toString(i) + ": " + toString(estimatedParams[i]) + " and true value: " + toString(trueParams[i-1]));
		}
	}
	logfile->done();

	//read estimated params for RG2
	logfile->startIndent("Checking parameter values for RG2");
	getline(in, tmp);
	fillVectorFromStringAny(tmp, estimatedParams, "\t");

	for(unsigned int i=1; i<estimatedParams.size()-1; ++i){ //first one is read group name, last one LL
		if(estimatedParams[i] != trueParams[i]){
			logfile->newLine();
			logfile->conclude("esimated value for parameter number " + toString(i) + ": " + toString(estimatedParams[i]) + " and true value: " + toString(trueParams[i-1]));
		}
	}
	logfile->done();
	return true;
}

//------------------------------------------
//TAtlasTest_BQSRSimulation
//------------------------------------------

TAtlasTest_BQSRSimulation::TAtlasTest_BQSRSimulation(TParameters & params, TLog* logfile):TAtlasTest(params, logfile){
	_name = "BQSRSimulation";
	filenameTag = _testingPrefix + _name;
	bamFileName = filenameTag + ".bam";
	fastaFileName = filenameTag + ".fasta";
	meanQual = params.getParameterIntWithDefault("BQSR_meanQual", 25);
	sdphredInt = params.getParameterDoubleWithDefault("BQSR_sdQual", 10);
	minPhredInt = params.getParameterIntWithDefault("BQSR_minQual", 0);
	maxPhredInt = params.getParameterIntWithDefault("BQSR_maxQual", 42);
	qualityDist = params.getParameterStringWithDefault("BQSR_qualityDist", "normal(" + toString(meanQual) + "," + toString(sdphredInt) + ")[" + toString(minPhredInt) + "," + toString(maxPhredInt) + "]");
//	alpha = params.getParameterDoubleWithDefault("alpha", 10.0);
//	beta = params.getParameterDoubleWithDefault("beta", 0.2);
	minReadLen = params.getParameterIntWithDefault("BQSR_minReadLen", 30);
	maxReadLen = params.getParameterIntWithDefault("BQSR_maxReadLen", 100);
//	readLengthDist = params.getParameterStringWithDefault("readLength", "gamma(alpha,beta)[min,max]");
	positionEffectSlope = params.getParameterDoubleWithDefault("BQSR_positionEffectSlope", 0.0144928);
	positionEffectIntercept = params.getParameterDoubleWithDefault("BQSR_positionEffectIntercept", 0.485507);
	phi1 = params.getParameterIntWithDefault("BQSR_phi1", 35);
	phi2 = params.getParameterDoubleWithDefault("BQSR_phi2", 1.2);
	revIntercept = params.getParameterDoubleWithDefault("BQSR_revIntercept", 1.5);
	acceptedDelta = params.getParameterDoubleWithDefault("BQSR_acceptedDelta", 1);

}

bool TAtlasTest_BQSRSimulation::run(){

	//TODO: find minimal data necessary to run test in order to speed up

	//1) Run ATLAS to simulate BAM file
	//-----------------------------
	_testParams.addParameter("out", filenameTag);
	_testParams.addParameter("qualityDist", qualityDist);
	_testParams.addParameter("chrLength", "5000000");
	_testParams.addParameter("refDiv", "0.0");
	_testParams.addParameter("ploidy", "1");
	_testParams.addParameter("BQSRTransformation", "[" + toString(phi1) + "," + toString(phi2) + "," + toString(revIntercept) + "]");
//	_testParams.addParameter("readLength", "gamma(" + toString(alpha) + "," + toString(beta)+ ")[" + toString(minReadLen) + "," + toString(maxReadLen));
	_testParams.addParameter("readLength", "fixed(70)");


	if(!runTGenomeFromInputfile("simulate"))
		return false;

	logfile->newLine();

	//1) Run BQSR
	//-----------------------------
	_testParams.addParameter("bam", bamFileName);
	_testParams.addParameter("fasta", fastaFileName);
	_testParams.addParameter("storeInMemory", "");
	_testParams.addParameter("estimateBQSRPosition", "");
	_testParams.addParameter("maxPos", "110");

	if(!runTGenomeFromInputfile("BQSR"))
		return false;


	//3) check if results are OK
	//--------------------------
	if(checkBQSRQualityFile() == true && checkBQSRPositionFile() == true) return true;
	else return false;
};

double TAtlasTest_BQSRSimulation::trueQual(int & phi1, double & phi2, int & fakeQual){
	double tmpPhi1 = (double) phi1;
	double tmpFakeQual = (double) fakeQual;
	double exp1, exp2;
	exp1 = pow(10.0,-1.0/10.0*phi2*tmpFakeQual);
	exp2 = pow(10.0, -tmpPhi1/10.0);

	double trueQual = -10.0 * log10(exp1 + exp2);
	return trueQual;
}

bool TAtlasTest_BQSRSimulation::checkBQSRQualityFile(){
	logfile->startIndent("Checking BQSR Quality table:");

	//open quality file
	std::string filename = filenameTag + "_BQSR_ReadGroup_Quality_Table.txt";
	logfile->listFlush("Opening file '" + filename + "' for reading ...");
	std::ifstream in(filename.c_str());
	if(!in)
		throw "Failed to open file '" + filename + "'!";
	logfile->done();

	//skip header and quality 0
	std::string tmp;
	getline(in, tmp);
	getline(in, tmp);

	//some variables
	std::vector<std::string> line;
	int numLines = 0;
	int QualityScoreAsPhredInt;
	double EmpiricalQuality;
	double Log10Observations;
	int unacceptablesCount = 0;
	double maxEmpiricQual = 0;

	//parse file line by line check contents
	logfile->listFlush("Parsing file ...");
	while(in.good() && !in.eof()){
		//read line into vector
		++numLines;
		fillVectorFromLineWhiteSpaceSkipEmpty(in, line);
		QualityScoreAsPhredInt = stringToInt(line[1]);
		EmpiricalQuality = stringToDouble(line[3]);
		Log10Observations = stringToDouble(line[4]);
		if(Log10Observations >= 5.5 && fabs(EmpiricalQuality - trueQual(phi1, phi2, QualityScoreAsPhredInt)) > acceptedDelta){
			++unacceptablesCount;
		}
		if(Log10Observations >= 4.5 && (EmpiricalQuality > maxEmpiricQual)){
			maxEmpiricQual = EmpiricalQuality;
		}
	}
	if(unacceptablesCount > 0){
		logfile->newLine();
		logfile->conclude("There were " + toString(unacceptablesCount) + " empirical quality scores that did not match.");
		return false;
	}
	if(fabs(maxEmpiricQual - phi1) > acceptedDelta){
		logfile->newLine();
		logfile->conclude("There is at least one empirical quality scores that was estimated to be larger than phi1.");
		return false;
	}
	logfile->done();
	logfile->endIndent();

	return true;
}

//void TAtlasTest_BQSRSimulation::calculateSlopeIntercept(){
//	std::map<std::string, std::string> readLengthMap;
//
//	double sum = 0.0;
//	//gamma density starts at 0 but p at 1!
//	for(int p=1; p<(maxReadLen + 1) ; ++p)
//		sum += (double) p * readLengthDist->positionProbs[p-1];
//
//	m = (1.0 - revIntercept) / (sum - maxReadLength);
//	intercept = revIntercept - m * maxReadLength;
//
//	if(intercept < 0) throw "The value given for the reverse intercept results in a negative intercept!";
//}

double TAtlasTest_BQSRSimulation::trueScaling(int & pos){
	double trueScaling = positionEffectIntercept + positionEffectSlope * pos;
	return trueScaling;
}

bool TAtlasTest_BQSRSimulation::checkBQSRPositionFile(){
	logfile->startIndent("Checking BQSR Position table:");

	//open quality file
	std::string filename = filenameTag + "_BQSR_ReadGroup_Position_Table.txt";
	logfile->listFlush("Opening file '" + filename + "' for reading ...");
	std::ifstream in(filename.c_str());
	if(!in)
		throw "Failed to open file '" + filename + "'!";
	logfile->done();

	//skip header and quality 0
	std::string tmp;
	getline(in, tmp);
	getline(in, tmp);

	//some variables
	std::vector<std::string> line;
	int numLines = 0;
	int Position;
	double Scaling;
	double Log10Observations;
	int unacceptablesCount = 0;

	//parse file line by line check contents
	logfile->listFlush("Parsing file ...");
	while(in.good() && !in.eof()){
		//read line into vector
		++numLines;
		fillVectorFromLineWhiteSpaceSkipEmpty(in, line);
		Position = stringToInt(line[1]);
		Scaling = stringToDouble(line[3]);
		Log10Observations = stringToDouble(line[4]);
		if(Log10Observations > 4.5 && fabs(trueScaling(Position) - Scaling) >= 0.1) ++unacceptablesCount;
	}
	if(unacceptablesCount > 0){
		logfile->newLine();
		logfile->conclude("There were " + toString(unacceptablesCount) + " scaling factor estimates that did not match.");
		return false;
	}
	logfile->done();
	logfile->endIndent();

	return true;
}

//------------------------------------------
//TAtlasTest_qualityTransformationRecal
//------------------------------------------

TAtlasTest_qualityTransformationRecalPlain::TAtlasTest_qualityTransformationRecalPlain(TParameters & params, TLog* logfile):TAtlasTest(params, logfile){
	_name = "qualityTransformation";
	filenameTag = _testingPrefix + _name;
	bamFileName = filenameTag + ".bam";
	std::string onlyRecalParams = "1,0;0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0";
	recalParamString = params.getParameterStringWithDefault("recal_recalParams","qualFuncPosFuncContext[" + onlyRecalParams +"]");
	maxReadLength = 70;
	randomGenerator = new TRandomGenerator();
	qualDistString = "10";
	qualityDist =  new TSimulatorQualityDist(qualDistString);
	recalObject = new TSimulatorQualityTransformationRecal("qualFuncPosFuncContext[" + onlyRecalParams + "]", maxReadLength, qualityDist, randomGenerator);
}

bool TAtlasTest_qualityTransformationRecalPlain::run(){
	//1) Run ATLAS to simulate BAM file
	//-----------------------------
	_testParams.addParameter("out", filenameTag);
	_testParams.addParameter("chrLength", "2000000");
	_testParams.addParameter("depth", "4");
	_testParams.addParameter("recal", recalParamString);
	_testParams.addParameter("readLength", "single:fixed("+toString(maxReadLength) + ")");
	_testParams.addParameter("qualityDist", "fixed(" + qualDistString + ")");


	if(!runTGenomeFromInputfile("simulate"))
		return false;

	logfile->newLine();

	//1) Run qualityTransformation
	//-----------------------------
	_testParams.addParameter("bam", bamFileName);
	_testParams.addParameter("recal", recalParamString);


	if(!runTGenomeFromInputfile("qualityTransformation"))
		return false;


	//3) check if results are OK
	//--------------------------
	if(readTransformationFile() == true){
		if(checkTransformation(qualDistVec) == true) return true;
	}
//	else return false;
	return false;
};

bool TAtlasTest_qualityTransformationRecalPlain::readTransformationFile(){
	//open quality file
	std::string filename = filenameTag + "_total_qualityTransformation.txt";
	logfile->listFlush("Opening file '" + filename + "' for parsing ...");
	std::ifstream in(filename.c_str());
	if(!in){
		throw "Failed to open file '" + filename + "'!";
	}
	//parse file line by line check contents
	std::vector<double> tmp;
	while(in.good() && !in.eof()){
		fillVectorFromLineAny(in, tmp, "\t");
		qualTransTable.push_back(tmp);
	}
	logfile->done();
	return true;
}

bool TAtlasTest_qualityTransformationRecalPlain::checkTransformation(std::vector<int> trueQualScores){
	//find true quality scores
	std::vector<int> transformedQualScores;
	int numQualScores = trueQualScores.size();

	transformedQualScores.push_back(19);
	transformedQualScores.push_back(30);
	transformedQualScores.push_back(40);
	transformedQualScores.push_back(60);

	//is the rest = 0?
	double s = 0.0;
	for(unsigned int i=1; i<qualTransTable.size(); ++i){
		for(unsigned int j=1; j<qualTransTable[i].size(); ++j){
			s += qualTransTable[i][j];
		}
	}
	if(s > 1.0001 || s < 0.9999){
		logfile->newLine();
		logfile->conclude("Proportions in qualityTransformation table don't sum to one!");
		return false;
	}

	//are the qualities transformed correctly? lines=transformed, columns=true
	double fracObservationsFound;
	double fracObservationsExpected = (double) (1.0 / (double) numQualScores);
	for(int qI=0; qI<(int) numQualScores; ++qI){
		fracObservationsFound = qualTransTable[trueQualScores[qI]+1][transformedQualScores[qI]+1];
		if( fracObservationsFound < (fracObservationsExpected - 0.0009) || fracObservationsFound > (fracObservationsExpected + 0.0009)){ //+1 for header and line names
			logfile->newLine();
			logfile->conclude("Wrong transformation of " + toString(trueQualScores[qI]+1) + "! Found " + toString(fracObservationsFound) + " observations in [" + toString(transformedQualScores[qI]+1) + "][" + toString(trueQualScores[qI]) + "] instead of " + toString(fracObservationsExpected) + "!");
			return false;
		}
	}
	return true;
}

//---------------------------------------
TAtlasTest_qualityTransformationRecalBinned::TAtlasTest_qualityTransformationRecalBinned(TParameters & params, TLog* logfile):TAtlasTest_qualityTransformationRecalPlain(params, logfile){
	recalParamString = params.getParameterStringWithDefault("recal_recalParams", "qualFuncPosFuncContext[2,0;0,0;0{20}]");
	qualDistString = "(10,15,20,30)";
	qualityDist =  new TSimulatorQualityDistBinned(qualDistString, randomGenerator);
	recalObject = new TSimulatorQualityTransformationRecal(recalParamString, maxReadLength, qualityDist, randomGenerator);
	fillVectorFromStringAnySkipEmpty(qualDistString, qualDistVec, ",");
}

bool TAtlasTest_qualityTransformationRecalBinned::run(){
	//1) Run ATLAS to simulate BAM file
	//-----------------------------
	_testParams.addParameter("out", filenameTag);
	_testParams.addParameter("chrLength", "2000000");
	_testParams.addParameter("depth", "2");
	_testParams.addParameter("recalTransformation", recalParamString);
	_testParams.addParameter("readLength", "single:fixed("+toString(maxReadLength) + ")");
	_testParams.addParameter("qualityDist", "binned(" + qualDistString + ")");

	if(!runTGenomeFromInputfile("simulate"))
		return false;

	logfile->newLine();

	//1) Run qualityTransformation
	//-----------------------------
	_testParams.addParameter("bam", bamFileName);
	_testParams.addParameter("recal", recalParamString);

	if(!runTGenomeFromInputfile("qualityTransformation"))
		return false;

	//3) check if results are OK
	//--------------------------
	if(readTransformationFile() == true){
		if(checkTransformation(qualDistVec) == true)
			return true;
	}
	return false;
}
