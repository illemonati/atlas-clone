/*
 * TAtlasTestPMD.cpp
 *
 *  Created on: Feb 22, 2018
 *      Author: vivian
 */

#include "TAtlasTestPMD.h"


TAtlasTest_PMDEmpiric::TAtlasTest_PMDEmpiric(TParameters & params, TLog* logfile):TAtlasTest(params, logfile){
	_name = "recalSimulation";
	filenameTag = _testingPrefix + _name;
	bamFileName = filenameTag + ".bam";
	fastaFileName = filenameTag + ".fasta";

	//open PMD param file
	pmdEmpiricFileName = filenameTag + "_true_PMD_input_Empiric.txt";
	poolRGFileName = filenameTag + "_poolThese.txt";

	alpha = params.getParameterDoubleWithDefault("recal_alpha", 10);
	beta = params.getParameterDoubleWithDefault("recal_beta", 0.2);
	minReadLength = params.getParameterDoubleWithDefault("recal_minReadLength", 30);
	maxReadLength = params.getParameterDoubleWithDefault("recal_maxReadLength", 100);

	firstPMDStringCT = params.getParameterStringWithDefault("PMD_firstPMDStringCT", "0.0570425,0.0136436,0.00557868,0.00369679,0.00217176,0.00165969,0.0013597,0.00128248,0.00115739,0.0010344,0.00100565,0.000947369,0.000914736,0.000946235,0.000931833,0.00090844,0.000981253,0.000969278,0.000952248,0.000983697,0.000967508,0.000912304,0.000872001,0.00106018,0.000734042,0.000768929,0.000640332,0.000614326,0.0005853,0.000602209,0.000544814,0.000542236,0.000562432,0.00060454,0.000559798,0.000579571,0.000613039,0.000606546,0.000512908,0.000409839,0.000504727,0.000540894,0.000660901,0.000449767,0.000563391,0.000545252,0.00056032,0.000537928,0.000578507,0.000568476");
	firstPMDStringGA = params.getParameterStringWithDefault("PMD_firstPMDStringGA", "0.0566675,0.0122868,0.00522836,0.00318553,0.00207032,0.00171248,0.00123768,0.00117464,0.000974231,0.000900184,0.000890755,0.000791377,0.000755095,0.000789783,0.000779903,0.000752922,0.000769629,0.000748781,0.000842946,0.000809977,0.000810973,0.000702441,0.000638453,0.000634849,0.000614554,0.000579583,0.000500633,0.000477177,0.000415798,0.00038988,0.000362353,0.000378398,0.000336733,0.000377509,0.00030131,0.000326029,0.000301191,0.00032936,0.00034176,0.000249148,0.000303415,0.000273776,0.000259635,0.000280269,0.000305902,0.000215535,0.000184103,0.000240266,0.000237668,0.000261064");
	secondPMDStringCT = params.getParameterStringWithDefault("PMD_secondPMDStringCT", "0.0738236,0.0169882,0.0068439,0.00512189,0.00298767,0.00235187,0.00226778,0.00181904,0.00188279,0.00159739,0.00144078,0.00162674,0.00151388,0.00157765,0.00148488,0.00164889,0.00164377,0.00157734,0.00156031,0.00193791,0.00207055,0.00186685,0.00188191,0.00189957,0.00183199,0.00181242,0.00183904,0.00174362,0.00160242,0.00182506,0.00181626,0.00188778,0.00192195,0.00239329,0.00153655,0.00169314,0.00173276,0.00167762,0.00166104,0.00169973,0.00192706,0.00201291,0.00179049,0.00154331,0.00164789,0.0020585,0.0017981,0.00170639,0.00193168,0.00173178");
	secondPMDStringGA = params.getParameterStringWithDefault("PMD_secondPMDStringGA", "0,0,0.00293742,0.00336413,0.0048651,0.00218997,0.00245297,0.00209862,0.00221646,0.00173175,0.00206741,0.00204128,0.00190718,0.00176191,0.00177644,0.00189132,0.00210642,0.0017206,0.00168852,0.00184077,0.00178956,0.00167079,0.00153733,0.00158146,0.0013248,0.00173105,0.00159785,0.00168088,0.00174343,0.00148742,0.00160376,0.00180732,0.00154377,0.00158591,0.0013558,0.00153654,0.00135425,0.00132407,0.00171275,0.00141973,0.00153039,0.00146885,0.00133669,0.00149182,0.00173479,0.00141099,0.00157696,0.00171272,0.00133564,0.00161568");
	thirdPMDStringCT = params.getParameterStringWithDefault("PMD_thirdPMDStringCT", "0.484202,0.293208,0.1472,0.0944159,0.0664974,0.0529067,0.0426632,0.0361048,0.0339675,0.0310974,0.0290736,0.0294713,0.029148,0.0290933,0.0294277,0.0302419,0.0309769,0.0333302,0.0339012,0.0362583,0.035577,0.035068,0.034547,0.0332448,0.0339795,0.0323574,0.0337345,0.0319812,0.0318846,0.03236,0.03354,0.0374068,0.0381569,0.0378421,0.0384532,0.03858,0.0375318,0.0383115,0.0378018,0.0377071,0.0371242,0.0382571,0.0368973,0.0373601,0.0371789,0.0361634,0.0367845,0.0363578,0.0361899,0.037026");
	thirdPMDStringGA = params.getParameterStringWithDefault("PMD_thirdPMDStringGA", "0.0521108,0.0522647,0.0571158,0.0489263,0.0444985,0.0392937,0.0392914,0.0380545,0.0373113,0.0372499,0.0371585,0.036443,0.0381607,0.0361248,0.0378769,0.0367767,0.0372836,0.0374682,0.0368751,0.0364539,0.0363257,0.0358006,0.0365954,0.0353917,0.0358281,0.0356717,0.0356513,0.0348792,0.0339361,0.034636,0.0348795,0.0338893,0.0341659,0.0333493,0.0336529,0.0324415,0.033668,0.032749,0.0324516,0.0324685,0.0312622,0.0320816,0.0322268,0.0321283,0.0313664,0.030891,0.0317229,0.0305144,0.0300429,0.0298577");

	CTpatterns[0] = firstPMDStringCT;
	CTpatterns[1] = secondPMDStringCT;
	CTpatterns[2] = thirdPMDStringCT;

	GApatterns[0] = firstPMDStringGA;
	GApatterns[1] = secondPMDStringGA;
	GApatterns[2] = thirdPMDStringGA;
}

bool TAtlasTest_PMDEmpiric::run(){
	//1) Write PMD params to file
	outPMD.open(pmdEmpiricFileName.c_str());
	if(!outPMD) throw "Failed to open file '" + pmdEmpiricFileName + "'!";

	outPMD << "HWI-ST558:341:C7R9TACXX_BAR8_UDG_less99	Empiric[" << firstPMDStringCT << "]	Empiric[" << firstPMDStringGA << "]" << std::endl;
	outPMD << "HWI-ST558:341:C7R9TACXX_BAR8_UDG_plus100 Empiric[" << secondPMDStringCT << "]	Empiric[" << secondPMDStringGA << "]" << std::endl;
	outPMD << "HWI-ST558:341:C7R9TACXX_BAR8_plus100 Empiric[" << thirdPMDStringCT << "]	Empiric[" << thirdPMDStringGA << "]" << std::endl;

	//TODO: find minimal data necessary to run test in order to speed up

	//2) Run ATLAS to simulate BAM file
	//-----------------------------

	_testParams.addParameter("out", filenameTag);
	_testParams.addParameter("chrLength", "10000000");
	_testParams.addParameter("ploidy", "2");
	_testParams.addParameter("depth", "2");
	_testParams.addParameter("readLength", "gamma(" + toString(alpha) + "," + toString(beta)+ ")[" + toString(minReadLength) + "," + toString(maxReadLength));
	_testParams.addParameter("pmdFile", pmdEmpiricFileName);

	if(!runTGenomeFromInputfile("simulate"))
		return false;

	logfile->newLine();

	//3) Run estimatePMD
	//-----------------------------
	//open pool read group file
	outPool.open(poolRGFileName.c_str());
	if(!outPool) throw "Failed to open file '" + poolRGFileName + "'!";
	outPool << "HWI-ST558:341:C7R9TACXX_BAR8_UDG_less99 HWI-ST558:341:C7R9TACXX_BAR8_UDG_plus100";

	_testParams.clear();
	_testParams.addParameter("bam", bamFileName);
	_testParams.addParameter("fasta", fastaFileName);
	_testParams.addParameter("poolReadGroups", poolRGFileName);

	if(!runTGenomeFromInputfile("estimatePMD"))
		return false;

	//4) check if results are OK
	//--------------------------
	if(checkPMDEmpiricFile() == true) return true;
	else return false;
};


bool TAtlasTest_PMDEmpiric::checkPMDEmpiricFile(){
	logfile->startIndent("Checking Empiric PMD file:");

	//open quality file
	std::string filename = filenameTag + "_PMD_input_Empiric.txt";
	logfile->listFlush("Opening file '" + filename + "' for reading ...");
	std::ifstream in(filename.c_str());
	if(!in)
		throw "Failed to open file '" + filename + "'!";
	logfile->done();

	//some variables
	std::string tmp;
	std::vector<std::string> tmpVec;
	std::string CTString, GAString;
	std::vector<double> CTestimated, GAestimated, CTtrue, GAtrue;

	for(int rg=0; rg<3; ++rg){
		//read estimated params
		logfile->list("Reading line for read group " + toString(rg));
		getline(in, tmp);

		//CT string
		std::string::size_type pos1 = tmp.find_first_of('[');
		if(pos1 == std::string::npos) throw "Can not find '[' in '" + tmp + "'!";
		std::string::size_type pos2 = tmp.find_first_of(']');
		if(pos2 == std::string::npos) throw "Can not find ']' in '" + tmp + "'!";
		CTString = tmp.substr((pos1+1),pos2-(pos1+1));
		fillVectorFromStringAny(CTString, CTestimated, ",");

		//GA string
		pos1 = tmp.find_last_of('[');
		if(pos1 == std::string::npos) throw "Can not find second '[' in '" + tmp + "'!";
		pos2 = tmp.find_last_of(']');
		if(pos2 == std::string::npos) throw "Can not find second ']' in '" + tmp + "'!";
		GAString = tmp.substr((pos1+1),pos2-(pos1+1));
		fillVectorFromStringAny(GAString, GAestimated, ",");

		//parse true params
		fillVectorFromStringAny(CTpatterns[0], CTtrue, ",");
		fillVectorFromStringAny(GApatterns[0], GAtrue, ",");

		//compare
		for(unsigned int i=0; i<CTtrue.size(); ++i){
			logfile->conclude("At pos " + toString(i) + " the true C to T damage proportion is = " + toString(CTtrue[i]) + " and it was estimated to be = " + toString(CTestimated[i]));
		}
	}

	return false;
}
