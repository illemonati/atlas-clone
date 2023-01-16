/*
 * estimHet.cpp
 *
 *  Created on: Feb 19, 2015
 *      Author: wegmannd
 */

#include <string>

#include "coretools/Main/TMain.h"

//BAM
#include "TBamDiagnoser.h"
#include "GenomeTasks/TBamFilter.h"
#include "TContextQuantifer.h"
#include "TDuplicateQuantifier.h"
#include "TPMDEstimator.h"
#include "TPMDSCalculator.h"
#include "TPileup.h"
#include "TQualityDistribution.h"
#include "TReadGroupMerger.h"
#include "TSoftClipping.h"
#include "TAlignmentMerger.h"

//window
#include "TAlleleCountEstimator.h"
#include "TAlleleFrequencyEstimator.h"
#include "TAllelicDepthCounts.h"
#include "TBamDownsampler.h"
#include "TCaller.h"
#include "TCreateBedMask.h"
#include "TDepthWriter.h"
#include "TDistanceEstimator.h"
#include "TEstimateRecalibration.h"
#include "TEstimateTheta.h"
#include "TF2Estimator.h"
#include "TGLF.h"
#include "THardyWeinbergTest.h"
#include "TInbreedingEstimator.h"
#include "TMajorMinor.h"
#include "TPSMCInput.h"
#include "TPolymorhicWindowIdentifier.h"
#include "TWriteGLF.h"
#include "TSexEstimator.h"

//VCF
#include "TVcfCompare.h"
#include "TVcfConverter.h"
#include "TVcfDiagnostics.h"

//simulations
#include "TSimulator.h"

#include "TReadGroupInfo.h"

/*
//tests
#include "TAtlasTest.h"
#include "TAtlasTestFilter.h"
#include "TAtlasTestMergePairs.h"
#include "TAtlasTestPMD.h"
#include "TAtlasTestRecalibration.h"
#include "TVcfTest.h"*/

void addTaks(coretools::TMain & main) {
    // Use main.addRegularTask to add a regular task (shown in list of available tasks)

	//BAM
	main.addRegularTask("filterBAM", new GenomeTasks::TTask_filterBAM());
	main.addRegularTask("splitMerge", new GenomeTasks::TTask_splitMerge());
	main.addRegularTask("mergeRG", new GenomeTasks::TTask_mergeReadGroups());
	main.addRegularTask("pileup", new GenomeTasks::TTask_pileup());
	main.addRegularTask("BAMDiagnostics", new GenomeTasks::TTask_BAMDiagnostics());
	main.addRegularTask("readOverlap", new GenomeTasks::TTask_overlapQuantifier());
	//main.addRegularTask("duplication", new GenomeTasks::TTask_duplicationQuantifier());
	main.addRegularTask("assessSoftClipping", new GenomeTasks::TTask_assessSoftClipping());
	main.addRegularTask("removeSoftClippedBases", new GenomeTasks::TTask_removeSoftClippedBasesFromReads());
	//main.addRegularTask("qualityDist", new GenomeTasks::TTask_qualityDist());
	main.addRegularTask("qualityTransformation", new GenomeTasks::TTask_qualityTransformation());
	//main.addRegularTask("contextStats", new GenomeTasks::TTask_quantifyContext());
	main.addRegularTask("downsample", new GenomeTasks::TTask_downsample());
	//main.addRegularTask("separateReads", new GenomeTasks::TTask_separateReads());
	main.addRegularTask("downsampleReads", new GenomeTasks::TTask_downSampleReads());
	main.addRegularTask("PMD", new GenomeTasks::TTask_estimatePMD());
	main.addRegularTask("PMDS", new GenomeTasks::TTask_PMDS());

	//window tasks
	main.addRegularTask("recal", new GenomeTasks::TTask_recal());
	//TODO: add task to calculate LL of full model
	//main.addRegularTask("writeDepth", new GenomeTasks::TTask_depthWriter());
	main.addRegularTask("createMask", new GenomeTasks::TTask_createMask());
	main.addRegularTask("allelicDepth", new GenomeTasks::TTask_allelicDepth());
	main.addRegularTask("PSMC", new GenomeTasks::TTask_PSMC());
	main.addRegularTask("call", new GenomeTasks::TTask_call());
	main.addRegularTask("theta", new GenomeTasks::TTask_estimateTheta());
	main.addRegularTask("thetaGenomeWide", new GenomeTasks::TTask_estimateThetaGenomeWide());
	main.addRegularTask("thetaRatio", new GenomeTasks::TTask_estimateThetaRatio());
	main.addRegularTask("thetaQC", new GenomeTasks::TTask_downsamplingThetaQC());
	main.addRegularTask("GLF", new GenomeTasks::TTask_writeGLF());
	main.addRegularTask("sexEstimation", new GenomeTasks::TTask_estimateSex());

	//Population tools
	main.addRegularTask("printGLF", new GLF::TTask_printGLF());
	main.addRegularTask("majorMinor", new PopulationTools::TTask_majorMinor());
	main.addRegularTask("geneticDist", new PopulationTools::TTask_estimateDist());
	main.addRegularTask("alleleCounts", new PopulationTools::TTask_estimateAlleleCounts());
	main.addRegularTask("alleleCountLikelihoods", new PopulationTools::TTask_writeAlleleCountLikelihoods);
	main.addRegularTask("transformAlleleCountFormat", new PopulationTools::TTask_transformAlleleCountFormat());
	main.addRegularTask("alleleFreq", new PopulationTools::TTask_estimateAlleleFreq());
	main.addRegularTask("compareAlleleFreq", new PopulationTools::TTask_compareAlleleFreq());
	main.addRegularTask("alleleFreqLikelihoods", new PopulationTools::TTask_alleleFreqLikelihoods());
	main.addRegularTask("inbreeding", new PopulationTools::TTask_estimateInbreeding());
	main.addRegularTask("polymorphicWindows", new PopulationTools::TTask_identifyPolymorphicWindows());
    main.addRegularTask("calculateF2", new PopulationTools::TTask_calculateF2());


    //VCF
	main.addRegularTask("VCFDiagnostics", new VCF::TTask_VCFDiagnostics());
	main.addRegularTask("VCFToInvariantBed", new VCF::TTask_VCFToInvariantBed());
	main.addRegularTask("convertVCF", new VCF::TTask_VcfConverter());
    main.addRegularTask("VCFFixInt", new VCF::TTask_VCFFixInt());
	main.addRegularTask("VCFCompare", new VCF::TTask_VcfCompare());
	main.addRegularTask("testHardyWeinberg", new PopulationTools::TTask_testHardyWeinberg);

	//simulations
	main.addRegularTask("simulate", new Simulations::TTask_simulate());

	// Use main.addDebugTask to add a debug task (not shown in list of available tasks)
	//main.addDebugTask("recalLL", new GenomeTasks::TTask_recalLL());
	main.addDebugTask("thetaLLSurface", new GenomeTasks::TTask_thetaLLSurface());
	BAM::RGInfo::TTask_testReadGroupInfo x;
	main.addDebugTask("json", new BAM::RGInfo::TTask_testReadGroupInfo());
};

void addTests(coretools::TMain & main){
    // Use testing.addTest to add a single test
	main.addTest("empty", new coretools::TTest());
	/*
    main.addTest("pileup", new TAtlasTest_pileup());
	main.addTest("allelicDepth", new TAtlasTest_allelicDepth());
	main.addTest("recalSimulation", new TAtlasTest_recalSimulation());
	main.addTest("BQSRSimulation", new TAtlasTest_BQSRSimulation());
	//main.addTest("qualityTransformationRecalPlain", new TAtlasTest_qualityTransformationRecalPlain());
	//main.addTest("qualityTransformationRecalBinned", new TAtlasTest_qualityTransformationRecalBinned());
	main.addTest("PMDEmpiric", new TAtlasTest_PMDEmpiric());
	main.addTest("theta", new TAtlasTest_theta());
	main.addTest("invariantBed", new TAtlasTest_invariantBed());
	main.addTest("mergePairs", new TAtlasTest_mergePairs());
	main.addTest("splitMerge", new TAtlasTest_mergeSplitPairs());
	main.addTest("filter", new TAtlasTest_filter());

    // Use testing.addTestSuite to add a test suite (= a combination of tests)
    // main.addTestSuite("exampleSuite", {"example", "empty"});

	//test suite "all"
	main.addTestSuite("all", {"pileup", "allelicDepth", "recalSimulation", "BQSRSimulation",
			"qualityTransformationRecalPlain", "qualityTransformationRecalBinned", "PMDEmpiric",
			"theta", "invariantBed", "mergePairs", "splitMerge", "filter"
	});

	//TODO: add suits for VCF, BAM, PopTools etc.
	*/
};

//---------------------------------------------------------------------------
//Main function
//---------------------------------------------------------------------------
int main(int argc, char* argv[]){
	coretools::TMain main("ATLAS", "0.9", "https://bitbucket.org/wegmannlab/atlas", "daniel.wegmann@unifr.ch");


	//add existing tasks
	addTaks(main);
    addTests(main);

	//now run program
	return main.run(argc, argv);
};


//---------------------------------------------------------------------------
// Missing TASKS
//---------------------------------------------------------------------------

/*

class TTask_splitRGbyLength:public TTask_atlas{ -> splitmerger
public:
	TTask_splitRGbyLength(){ _explanation = "Splitting single end read groups in a BAM file"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenomeWindows genome(logfile, parameters, _randomGenerator);
		genome.splitSingleEndReadGroups(parameters);
	};
};

class TTask_mateInfo:public TTask_atlas{ -> pileup
public:
	TTask_mateInfo(){ _explanation = "Writing read information per site"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenomeWindows genome(logfile, parameters, _randomGenerator);
		genome.printMateInformationPerSite(parameters);
	};
};

class TTask_writeDepthPerSite:public TTask_atlas{ -> pileup
public:
	TTask_writeDepthPerSite(){ _explanation = "Writing sequencing depth for each site"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenomeWindows genome(logfile, parameters, _randomGenerator);
		genome.writeDepthPerSite(parameters);
	};
};


class TTask_BQSR:public TTask_atlas{ -> suppress?
public:
	TTask_BQSR(){
		_explanation = "Estimating BQSR error re-calibration parameters";
		_citations.push_back("Hofmanova et al. (2016) PNAS");
	};

	void run(TParameters & parameters, TLog* logfile){
		TGenomeWindows genome(logfile, parameters, _randomGenerator);
		genome.BQSR(parameters);
	};
};

class TTask_binQualityScores:public TTask_atlas{ -> a feature of all BAM writing functions, e.g. BamFilter
public:
	TTask_binQualityScores(){ _explanation = "Binning quality scores"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenomeWindows genome(logfile, parameters, _randomGenerator);
		genome.binQualityScores(parameters);
	};
};


class TTask_recalBAM:public TTask_atlas{ -> feature of all BAM writing functions that parse, e.g. BamFilter
public:
	TTask_recalBAM(){ _explanation = "Recalibrating quality scores in a BAM file"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenomeWindows genome(logfile, parameters, _randomGenerator);
		genome.recalibrateBamFile(parameters);
	};
};

class TTask_testBED:public TTask_atlas{ ->should be a proper integration test
public:
	TTask_testBED(){ _explanation = "Testing BED files"; };

	void run(TParameters & parameters, TLog* logfile){
		TBed bed;
		bed.test();
	};
};


*/
