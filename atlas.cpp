/*
 * estimHet.cpp
 *
 *  Created on: Feb 19, 2015
 *      Author: wegmannd
 */

#include "TMain.h"
#include "allTasks.h"
#include "Tests/TAtlasTest.h"
#include "TAtlasTestFilter.h"
#include "TAtlasTestMergePairs.h"
#include "TAtlasTestPMD.h"
#include "TAtlasTestRecalibration.h"
#include "TVcfTest.h"

void addTaks(TMain & main) {
    // Use main.addRegularTask to add a regular task (shown in list of available tasks)
	main.addRegularTask("splitRGByLength", new TTask_splitRGbyLength());
	main.addRegularTask("mergeRG", new TTask_mergeReadGroups());
	main.addRegularTask("pileup", new TTask_pileup());
	main.addRegularTask("BAMDiagnostics", new TTask_BAMDiagnostics());
	main.addRegularTask("readOverlap", new TTask_assessReadOverlap());
	main.addRegularTask("mergeReads", new TTask_mergeReads());
	main.addRegularTask("splitMerge", new TTask_splitMerge());
	main.addRegularTask("duplication", new TTask_assessDuplication());
	main.addRegularTask("writeReadInfoPerSite", new TTask_mateInfo());
	main.addRegularTask("assessSoftClipping", new TTask_assessSoftClipping());
	main.addRegularTask("contextStats", new TTask_contextStats());
	main.addRegularTask("removeSoftClippedBases", new TTask_removeSoftClippedBasesFromReads());
	main.addRegularTask("filter", new TTask_filterBAM());
	main.addRegularTask("writeDepthPerWindow", new TTask_writeDepthPerWindow());
	main.addRegularTask("depthPerSiteDist", new TTask_depthPerSiteDist());
	main.addRegularTask("writeDepthPerSite", new TTask_writeDepthPerSite());
	main.addRegularTask("createDepthMask", new TTask_createDepthMask());
	main.addRegularTask("createNonRefMask", new TTask_createNonRefMask());
	main.addRegularTask("allelicDepth", new TTask_allelicDepth());
	main.addRegularTask("downsample", new TTask_downsample());
	main.addRegularTask("separateReads", new TTask_separateReads());
	main.addRegularTask("downsampleReads", new TTask_downSampleReads());
	main.addRegularTask("PMD", new TTask_estimatePMD());
	main.addRegularTask("PMDS", new TTask_PMDS());
	main.addRegularTask("PSMC", new TTask_PSMC());
	main.addRegularTask("recal", new TTask_recal());
	main.addRegularTask("BQSR", new TTask_BQSR());
	main.addRegularTask("BAMUpdateQualities", new TTask_recalBAM());
	main.addRegularTask("qualityDist", new TTask_qualityDist());
	main.addRegularTask("qualityTransformation", new TTask_qualityTransformation());
	main.addRegularTask("binQualityScores", new TTask_binQualityScores());
	main.addRegularTask("call", new TTask_call());
	main.addRegularTask("theta", new TTask_estimateTheta());
	main.addRegularTask("thetaRatio", new TTask_estimateThetaRatio());
	main.addRegularTask("thetaQC", new TTask_downsamplingThetaQC());
	main.addRegularTask("GLF", new TTask_GLF());
	main.addRegularTask("printGLF", new TTask_printGLF());
	main.addRegularTask("majorMinor", new TTask_majorMinor());
	main.addRegularTask("geneticDist", new TTask_estimateDist());
	main.addRegularTask("alleleCounts", new TTask_estimateAlleleCounts());
	main.addRegularTask("transformCountFormat", new TTask_transformAlleleCountFormat());
	main.addRegularTask("alleleFreq", new TTask_estimateAlleleFreq());
	main.addRegularTask("compareAlleleFreq", new TTask_compareAlleleFreq());
	main.addRegularTask("alleleFreqLikelihoods", new TTask_alleleFreqLikelihoods());
	main.addRegularTask("inbreeding", new TTask_estimateInbreeding());
	main.addRegularTask("VCFAssessAllelicBalance", new TTask_VCFDiagnostics());
	main.addRegularTask("VCFToInvariantBed", new TTask_VCFToInvariantBed());
	main.addRegularTask("VCFToBeagle", new TTask_VCFToBeagle());
    main.addRegularTask("VCFToLFMM", new TTask_VCFToLFMM());
    main.addRegularTask("VCFToPosFile", new TTask_VCFToPosFile());
    main.addRegularTask("VCFToGenotypeTruthSetFile", new TTask_VCFToGenotypeTruthSetFile());
    main.addRegularTask("VCFFixInt", new TTask_VCFFixInt());
	main.addRegularTask("VCFCompare", new TTask_VCFCompare());
	main.addRegularTask("simulate", new TTask_simulate());

	//debug tasks
	main.addDebugTask("recalLL", new TTask_recalLL());
	main.addDebugTask("inbreedingLikelihood", new TTask_inbreedingLikelihood());
	main.addDebugTask("thetaLLSurface", new TTask_thetaLLSurface());
	main.addDebugTask("alleleFrequencyLikelihoods", new TTask_writeAlleleFrequencyLikelihoods);
	main.addDebugTask("polymorphicWindows", new TTask_identifyPolymorphicWindows());
	main.addDebugTask("testBED", new TTask_testBED()); //TODO: write as test!
	main.addDebugTask("testGL", new TTask_testGL());
};

void addTests(TMain & main){
    // Use testing.addTest to add a single test
    main.addTest("empty", new TTest());
    main.addTest("pileup", new TAtlasTest_pileup());
	main.addTest("allelicDepth", new TAtlasTest_allelicDepth());
	main.addTest("recalSimulation", new TAtlasTest_recalSimulation());
	main.addTest("BQSRSimulation", new TAtlasTest_BQSRSimulation());
	main.addTest("qualityTransformationRecalPlain", new TAtlasTest_qualityTransformationRecalPlain());
	main.addTest("qualityTransformationRecalBinned", new TAtlasTest_qualityTransformationRecalBinned());
	main.addTest("PMDEmpiric", new TAtlasTest_PMDEmpiric());
	main.addTest("theta", new TAtlasTest_theta());
	main.addTest("invariantBed", new TAtlasTest_invariantBed());
	main.addTest("mergePairs", new TAtlasTest_mergePairs());
	main.addTest("splitMerge", new TAtlasTest_mergeSplitPairs());
	main.addTest("filter", new TAtlasTest_filter());

    // Use testing.addTestSuite to add a test suite (= a combination of tests)
//    main.addTestSuite("exampleSuite", {"example", "empty"});

	//test suite "all"
	main.addTestSuite("all", {"pileup", "allelicDepth", "recalSimulation", "BQSRSimulation",
			"qualityTransformationRecalPlain", "qualityTransformationRecalBinned", "PMDEmpiric",
			"theta", "invariantBed", "mergePairs", "splitMerge", "filter"
	});

};

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

