/*
 * estimHet.cpp
 *
 *  Created on: Feb 19, 2015
 *      Author: wegmannd
 */

#include "TMain.h"
#include "allTasks.h"
#include "Tests/TAtlasTest.h"

void addTaks(TMain & main) {
    // Use main.addRegularTask to add a regular task (shown in list of available tasks)
    main.addRegularTask("simulate", new TTask_simulate());

/*

	taskMap_regular["splitRGByLength"] = new TTask_splitRGbyLength();
	taskMap_regular["mergeRG"] = new TTask_mergeReadGroups();
	taskMap_regular["pileup"] = new TTask_pileup();
	taskMap_regular["BAMDiagnostics"] = new TTask_BAMDiagnostics();
	taskMap_regular["readOverlap"] = new TTask_assessReadOverlap();
	taskMap_regular["mergeReads"] = new TTask_mergeReads();
	taskMap_regular["splitMerge"] = new TTask_splitMerge();
	taskMap_regular["duplication"] = new TTask_assessDuplication();
	taskMap_regular["writeReadInfoPerSite"] = new TTask_mateInfo();
	taskMap_regular["assessSoftClipping"] = new TTask_assessSoftClipping();
	taskMap_regular["contextStats"] = new TTask_contextStats();
	taskMap_regular["removeSoftClippedBases"] = new TTask_removeSoftClippedBasesFromReads();
	taskMap_regular["filter"] = new TTask_filterBAM();
	taskMap_regular["writeDepthPerWindow"] = new TTask_writeDepthPerWindow();
	taskMap_regular["depthPerSiteDist"] = new TTask_depthPerSiteDist();
	taskMap_regular["writeDepthPerSite"] = new TTask_writeDepthPerSite();
	taskMap_regular["createDepthMask"] = new TTask_createDepthMask();
	taskMap_regular["createNonRefMask"] = new TTask_createNonRefMask();
	taskMap_regular["allelicDepth"] = new TTask_allelicDepth();
	taskMap_regular["downsample"] = new TTask_downsample();
	taskMap_regular["separateReads"] = new TTask_separateReads();
	taskMap_regular["downsampleReads"] = new TTask_downSampleReads();
	taskMap_regular["PMD"] = new TTask_estimatePMD();
	taskMap_regular["PMDS"] = new TTask_PMDS();
	taskMap_regular["PSMC"] = new TTask_PSMC();
	taskMap_regular["recal"] = new TTask_recal();
	taskMap_regular["recalLL"] = new TTask_recalLL();
	taskMap_regular["BQSR"] = new TTask_BQSR();
	taskMap_regular["BAMUpdateQualities"] = new TTask_recalBAM();
	taskMap_regular["qualityDist"] = new TTask_qualityDist();
	taskMap_regular["qualityTransformation"] = new TTask_qualityTransformation();
	taskMap_regular["binQualityScores"] = new TTask_binQualityScores();
	taskMap_regular["call"] = new TTask_call();
	taskMap_regular["theta"] = new TTask_estimateTheta();
	taskMap_regular["thetaRatio"] = new TTask_estimateThetaRatio();
	taskMap_regular["thetaQC"] = new TTask_downsamplingThetaQC();
	taskMap_regular["GLF"] = new TTask_GLF();
	taskMap_regular["printGLF"] = new TTask_printGLF();
	taskMap_regular["majorMinor"] = new TTask_majorMinor();
	taskMap_regular["geneticDist"] = new TTask_estimateDist();
	taskMap_regular["alleleCounts"] = new TTask_estimateAlleleCounts();
	taskMap_regular["transformCountFormat"] = new TTask_transformAlleleCountFormat();
	taskMap_regular["alleleFreq"] = new TTask_estimateAlleleFreq();
	taskMap_regular["inbreeding"] = new TTask_estimateInbreeding();
	taskMap_regular["VCFAssessAllelicBalance"] = new TTask_VCFDiagnostics();
	taskMap_regular["VCFToInvariantBed"] = new TTask_VCFToInvariantBed();
	taskMap_regular["VCFToBeagle"] = new TTask_VCFToBeagle();
    taskMap_regular["VCFToLFMM"] = new TTask_VCFToLFMM();
    taskMap_regular["VCFFixInt"] = new TTask_VCFFixInt();
	taskMap_regular["VCFCompare"] = new TTask_VCFCompare();
	taskMap_regular["simulate"] = new TTask_simulate();
	//taskMap_regular["filterVCF"] = new TTask_filterVCF();

	//and debug tasks
	taskMap_debug["inbreedingLikelihood"] = new TTask_inbreedingLikelihood();
	taskMap_debug["thetaLLSurface"] = new TTask_thetaLLSurface();
	taskMap_debug["alleleFrequencyLikelihoods"] = new TTask_writeAlleleFrequencyLikelihoods();
	taskMap_debug["polymorphicWindows"] = new TTask_identifyPolymorphicWindows();
	taskMap_debug["testBED"] = new TTask_testBED(); //TODO: write as test!


*/
    // Use main.addDebugTask to add a debug task (not shown in list of available tasks)
    // main.addDaddTaksebugTask("debugTask", new TTask_debugTask());
}

void addTests(TMain & main){
    // Use testing.addTest to add a single test
    main.addTest("empty", new TTest());
    main.addTest("pileup", new TAtlasTest_pileup());
    // Use testing.addTestSuite to add a test suite (= a combination of tests)
    main.addTestSuite("exampleSuite", {"example", "empty"});
}

/*
TAtlasTestList(){
	//fill map of tests
	//NOTE: order will be the order in which test will be run if initialized from TParameters
	testMap.emplace_back("empty", &createInstance<TAtlasTest>);
	testMap.emplace_back("pileup", &createInstance<TAtlasTest_pileup>);
	testMap.emplace_back("allelicDepth", &createInstance<TAtlasTest_allelicDepth>);
	testMap.emplace_back("recalSimulation", &createInstance<TAtlasTest_recalSimulation>);
	testMap.emplace_back("BQSRSimulation", &createInstance<TAtlasTest_BQSRSimulation>);
	testMap.emplace_back("qualityTransformationRecalPlain", &createInstance<TAtlasTest_qualityTransformationRecalPlain>);
	testMap.emplace_back("qualityTransformationRecalBinned", &createInstance<TAtlasTest_qualityTransformationRecalBinned>);
	testMap.emplace_back("PMDEmpiric", &createInstance<TAtlasTest_PMDEmpiric>);
	testMap.emplace_back("theta", &createInstance<TAtlasTest_theta>);
	testMap.emplace_back("invariantBed", &createInstance<TAtlasTest_invariantBed>);
	testMap.emplace_back("mergePairs", &createInstance<TAtlasTest_mergePairs>);
	testMap.emplace_back("splitMerge", &createInstance<TAtlasTest_mergeSplitPairs>);
	testMap.emplace_back("filter", &createInstance<TAtlasTest_filter>);

	//fill map of test suites
	//Note: suites and tests within suites will be initialized in this order!
	testSuites["exampleSuite"] = {"empty", "pileup"};

	//automatically create a test suite "all"
	testSuites["all"] = {};
	for(testMapIt = testMap.begin(); testMapIt != testMap.end(); ++testMapIt)
		testSuites["all"].push_back(testMapIt->first);

}
*/
//---------------------------------------------------------------------------
//Main function
//---------------------------------------------------------------------------
int main(int argc, char* argv[]){
	TMain main("Atlas", "0.9", "https://bitbucket.org/wegmannlab/atlas");

	//add existing tasks
	addTaks(main);
    addTests(main);

	//now run program
	return main.run(argc, argv);
};

