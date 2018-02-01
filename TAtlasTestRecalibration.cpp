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
	recalParamString = params.getParameterStringWithDefault("recal_recalParams", "2,0,0.1,0.001,1{20}");
}

bool TAtlasTest_recalSimulation::run(){

	//TODO: find minimal data necessary to run test in order to speed up

	//1) Run ATLAS to simulate BAM file
	//-----------------------------
	_testParams.addParameter("out", filenameTag);
	_testParams.addParameter("qualityDist", qualityDist);
	_testParams.addParameter("chrLength", "10000000");
	_testParams.addParameter("ploidy", "1");
	_testParams.addParameter("recalTransformation", "recal[" + recalParamString + "]");
//	_testParams.addParameter("readLength", "gamma(" + toString(alpha) + "," + toString(beta)+ ")[" + toString(minReadLen) + "," + toString(maxReadLen));
	_testParams.addParameter("readLength", "fixed(70)");


	if(!runTGenomeFromInputfile("simulate"))
		return false;

	logfile->newLine();

	//1) Run recal
	//-----------------------------
	_testParams.addParameter("bam", bamFileName);
//	_testParams.addParameter("recal", filenameTag + "_recalibrationEM.txt");


	if(!runTGenomeFromInputfile("recal"))
		return false;


	//3) check if results are OK
	//--------------------------
	if(checkRecalFile() == true) return true;
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
	std::string tmp;
	std::vector<double> estimatedParams;
	std::vector<std::string> tmpVec;
	std::vector<double> trueParams;


	//skip header
	getline(in, tmp);

	//read estimated params
	getline(in, tmp);
	fillVectorFromStringAny(tmp, estimatedParams, "\t");

	//parse true params
	fillVectorFromStringAnySkipEmpty(recalParamString, tmpVec, ",");
	repeatIndexes(tmpVec, trueParams);

	for(unsigned int i=1; i<estimatedParams.size(); ++i){ //first one is read group name
		if(estimatedParams[i] != trueParams[i]){
			logfile->newLine();
			logfile->conclude("esimated value for parameter number " + toString(i) + ": " + toString(estimatedParams[i]) + " and true value: " + toString(trueParams[i-1]));
		}

	}
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
	int QualityScore;
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
		QualityScore = stringToInt(line[1]);
		EmpiricalQuality = stringToDouble(line[3]);
		Log10Observations = stringToDouble(line[4]);
		if(Log10Observations >= 5.5 && fabs(EmpiricalQuality - trueQual(phi1, phi2, QualityScore)) > acceptedDelta){
			std::cout << QualityScore << " "<<EmpiricalQuality << " " << trueQual(phi1, phi2, QualityScore) << std::endl;
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
		std::cout << maxEmpiricQual << std::endl;
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
//TAtlasTest_qualityTransformation
//------------------------------------------

TAtlasTest_qualityTransformationRecalPlain::TAtlasTest_qualityTransformationRecalPlain(TParameters & params, TLog* logfile):TAtlasTest(params, logfile){
	_name = "qualityTransformation";
	filenameTag = _testingPrefix + _name;
	bamFileName = filenameTag + ".bam";
	recalParamString = params.getParameterStringWithDefault("recal_recalParams", "1,0{23}");
	maxReadLength = 70;
	randomGenerator = TRandomGenerator();
	binnedQualDist = "10";
	qualityDist = TSimulatorQualityDistBinned(binnedQualDist, & randomGenerator);
	recalObject = TSimulatorQualityTransformationRecal(recalParamString, maxReadLength, & qualityDist, & randomGenerator);

	//parse true params
	std::vector<std::string> tmpVec;
	fillVectorFromStringAnySkipEmpty(recalParamString, tmpVec, ",");
	repeatIndexes(tmpVec, trueParams);

}

bool TAtlasTest_qualityTransformationRecalPlain::run(){

	//TODO: find minimal data necessary to run test in order to speed up

	//1) Run ATLAS to simulate BAM file
	//-----------------------------
	_testParams.addParameter("out", filenameTag);
	_testParams.addParameter("chrLength", "2000000");
	_testParams.addParameter("depth", "2");
	_testParams.addParameter("ploidy", "2");
	_testParams.addParameter("recalTransformation", "recal[" + recalParamString + "]");
	_testParams.addParameter("readLength", "fixed("+toString(maxReadLength) + ")");
	_testParams.addParameter("qualityDist", "fixed(10)");


	if(!runTGenomeFromInputfile("simulate"))
		return false;

	logfile->newLine();

	//1) Run qualityTransformation
	//-----------------------------
	_testParams.addParameter("bam", bamFileName);
	_testParams.addParameter("recal", "recal[" + recalParamString + "]");


	if(!runTGenomeFromInputfile("qualityTransformation"))
		return false;


	//3) check if results are OK
	//--------------------------
	if(readTransformationFile() == true){

	}
//	else return false;
	return false;
};

bool TAtlasTest_qualityTransformationRecalPlain::readTransformationFile(){
	//open quality file
	std::string filename = filenameTag + "_total_qualityTransformation.txt";
	logfile->listFlush("Opening file '" + filename + "' for reading ...");
	std::ifstream in(filename.c_str());
	if(!in){
		throw "Failed to open file '" + filename + "'!";
	}
	//parse file line by line check contents
	std::vector<double> tmp;
	logfile->listFlush("Parsing file ...");
	while(in.good() && !in.eof()){
		fillVectorFromLineAny(in, tmp, "\t");
		qualTransTable.push_back(tmp);
	}
	logfile->done();
	return true;
}

bool TAtlasTest_qualityTransformationRecalPlain::checkTransformationBinned(std::vector<int> binnedQualScores){
	//find true quality scores
	std::vector<int> transformedQualScores;
	for(int i=0; binnedQualScores.size(); ++i){
		transformedQualScores.push_back(recalObject.getTransformedQuality(binnedQualScores[i],0,0));
/*
		fakeError = pow(10, (double) binnedQualScores[i]/-10.0);
		transFakeError = log(fakeError /(1.0 -fakeError));
		transFakeError = transFakeError * trueParams[0] + transFakeError * trueParams[0];
		for(unsigned int p=4; p<trueParams.size(); ++p){
			transFakeError += trueParams[p];
		}
		trueQualScores[i] = round(-10.0 * log10(transFakeError));
*/
	}

	//are the qualities transformed correctly?
	double fracObservations = 1.0 / binnedQualScores.size();
	for(unsigned int qI=0; qI<binnedQualScores.size(); ++qI){
		if(qualTransTable[binnedQualScores[qI]][transformedQualScores[qI]] != fracObservations){
			logfile->newLine();
			logfile->conclude("Wrong transformation of " + toString(binnedQualScores[qI]) + "!");
		}
	}

	//is the rest = 0?
	double s = 0.0;
	for(unsigned int i=0; i<qualTransTable.size(); ++i){
		for(unsigned int j=0; j<qualTransTable[0].size(); ++j){
			s += qualTransTable[i][j];
		}
	}
	if(s != 1.0){
		logfile->newLine();
		logfile->conclude("Proportions in qualityTransformation table don't sum to one!");
		return false;
	}

	return true;
}

//---------------------------------------
TAtlasTest_qualityTransformationRecalBinned::TAtlasTest_qualityTransformationRecalBinned(TParameters & params, TLog* logfile):TAtlasTest_qualityTransformationRecalPlain(params, logfile){
	recalParamString = params.getParameterStringWithDefault("recal_recalParams", "2,0{23}");
}

bool TAtlasTest_qualityTransformationRecalBinned::run(){
	return true;
}
