/*
 * estimHet.cpp
 *
 *  Created on: Feb 19, 2015
 *      Author: wegmannd
 */

#include <string>

#include "TEstimateErrors.h"
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
#include "TIlluminaIdentifier.h"

//window
#include "TAlleleCountEstimator.h"
#include "TAlleleFrequencyEstimator.h"
#include "TAllelicDepthCounts.h"
#include "TBamDownsampler.h"
#include "TCaller.h"
#include "TCreateBedMask.h"
#include "TDepthWriter.h"
#include "TDistanceEstimator.h"
#include "TEstimateErrors.h"
#include "TEstimateMutationLoad.h"
#include "TEstimateRecalibration.h"
#include "TEstimateTheta.h"
#include "TF2Estimator.h"
#include "TGLF.h"
#include "THardyWeinbergTest.h"
#include "TInbreedingEstimator.h"
#include "TMajorMinor.h"
#include "TPSMCInput.h"
#include "TPolymorhicWindowIdentifier.h"
#include "TSexEstimator.h"
#include "TWriteGLF.h"

//VCF
#include "TVcfCompare.h"
#include "TVcfConverter.h"
#include "TVcfDiagnostics.h"

//simulations
#include "TSimulator.h"

#include "TReadGroupInfo.h"
#include "TAncestralAlleleEstimator.h"

void addTaks(coretools::TMain & main) {
	//BAM
	main.createRegularTask<GenomeTasks::BamFilter::TBamFilter>("filterBAM", "Writing reads that pass filters to BAM file");
	main.createRegularTask<GenomeTasks::AlignmentMerger::TAlignmentSplitMerger>("splitMerge", "Splitting single-end reads and merging paired-end reads in BAM file");
	main.createRegularTask<GenomeTasks::TReadGroupMerger>("mergeRG", "Merging read groups in a BAM file");
	main.createRegularTask<GenomeTasks::TPileup>("pileup", "Printing pileup from BAM file");
	main.createRegularTask<GenomeTasks::TBamDiagnoser>("BAMDiagnostics", "Estimating approximate depth, read length frequencies and mapping quality frequencies");
	main.createRegularTask<GenomeTasks::AlignmentMerger::TOverlapQuantifier>("readOverlap", "Estimating distribution of overlap of paired reads in BAM file");
	main.createRegularTask<GenomeTasks::TAssessSoftClipping>("assessSoftClipping", "Assessing level of soft clipping in BAM file");
	main.createRegularTask<GenomeTasks::TRemoveSoftClippedBases>("removeSoftClippedBases", "Removing soft clipped bases from reads");
	main.createRegularTask<GenomeTasks::TQualityTransformation>("qualityTransformation", "Printing Quality Transformation");
	main.createRegularTask<GenomeTasks::TBamDownsampler>("downsample", "Downsampling a BAM file");
	main.createRegularTask<GenomeTasks::TPMDEstimator>("PMD", "Estimating Post-Mortem Damage (PMD) patterns");
	main.createRegularTask<GenomeTasks::TPMDSCalculator>("PMDS", "Filtering for ancient reads using PMDS", "Skoglund et al. (2014) PNAS");
	main.createRegularTask<GenomeTasks::TIlluminaIdentifier>("identifyIlluminaReadGroups", "Reassigning read groups based on the platform unit in their name");
	//main.addRegularTask("duplication", new GenomeTasks::TTask_duplicationQuantifier());
	//main.addRegularTask("qualityDist", new GenomeTasks::TTask_qualityDist());
	//main.addRegularTask("contextStats", new GenomeTasks::TTask_quantifyContext());
	//main.addRegularTask("separateReads", new GenomeTasks::TTask_separateReads());

	//window tasks
	main.createRegularTask<GenomeTasks::TEstimateErrors>("estimateErrors", "Estimating PMD pattern and Sequencing Errors", "Kousathanas et al. (2017) Genetics");
	main.createRegularTask<GenomeTasks::TEstimateRecalibration>("recal", "Estimating error re-calibration parameters", "Kousathanas et al. (2017) Genetics");
	main.createRegularTask<GenomeTasks::TMaskCreator>("createMask", "Creating a mask BED file");
	main.createRegularTask<GenomeTasks::TAllelicDepth>("allelicDepth", "Writing genotype likelihoods to a GLF file");
	main.createRegularTask<GenomeTasks::TPSMCInput>("PSMC", "Generating a PSMC Input file probabilistically");
	main.createRegularTask<GenomeTasks::TCall>("call", "Calling genotypes");
	main.createRegularTask<GenomeTasks::TEstimateTheta>("theta", "Estimating heterozygosity (theta)", "Kousathanas et al. (2017) Genetics");
	main.createRegularTask<GenomeTasks::TEstimateThetaRatio>("thetaRatio", "Estimate the ratio in heterozygosity (theta) between genomic regions", "Kousathanas et al. (2017) Genetics");
	main.createRegularTask<GenomeTasks::TWriteGLF>("GLF", "Writing genotype likelihoods to a GLF file");
	main.createRegularTask<GenomeTasks::TSexEstimator>("sexEstimation", "Estimating the distribution of depth among sites and writing depth per window");
	main.createRegularTask<GenomeTasks::TEstimateMutationLoad>("mutationLoad", "Estimating mutation load across the genome");
	//main.addRegularTask("writeDepth", new GenomeTasks::TTask_depthWriter());

	//Population tools
	main.createRegularTask<GLF::TGLFPrinter>("printGLF", "Printing a GLF file to screen");
	main.createRegularTask<PopulationTools::TMajorMinor>("majorMinor", "Estimating major and minor alles", "Skotte et al. (2012) Genetic Epidemiology");
	main.createRegularTask<PopulationTools::TDistanceEstimator>("geneticDist", "Estimating the genetic distance between individuals");
	main.createRegularTask<PopulationTools::TAlleleCounter>("alleleCounts", "Estimating population allele counts");
	main.createRegularTask<PopulationTools::TAlleleFreqEstimator>("alleleFreq", "Estimating population allele frequencies");
	main.createRegularTask<PopulationTools::TInbreedingEstimator>("inbreeding", "Estimating the inbreeding coefficient", "Burger et al. (2020) Current Biology");
	main.createRegularTask<PopulationTools::TPolymorhicWindowIdentifier>("polymorphicWindows", "Identifying windows for which samples are polymorphic");
    main.createRegularTask<PopulationTools::TF2Estimator>("calculateF2", "Calculate F2 between different samples, and within and between populations");
	main.createRegularTask<PopulationTools::TAncestralAlleleEstimator>("ancestralAlleles", "Writing FASTA-file with ancestral alleles");

    //VCF
	main.createRegularTask<VCF::TVcfDiagnostics>("VCFDiagnostics", "Diagnosing a VCF file");
	main.createRegularTask<VCF::TVCFConverter>("convertVCF", "Converting a VCF file to other formats");
	main.createRegularTask<VCF::TVcfCompare>("VCFCompare", "Comparing genotype calls in two VCF files");
	main.createRegularTask<PopulationTools::THardyWeinbergTest>("testHardyWeinberg", "Testing for Hardy-Weinberg equilibrium across multiple populations");

	//simulations
	main.createRegularTask<Simulations::TSimulationRunner>("simulate", "Generating simulations");

	// Debug tasks
	main.createDebugTask<GenomeTasks::TEstimateThetaLLSurface>("thetaLLSurface", "Calculating the theta LL surface for each window");
	main.createDebugTask<BAM::RGInfo::TReadGroupInfoTest>("json", "Testing JSON stuff");
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

	void run(){
		TGenomeWindows genome;
		genome.splitSingleEndReadGroups();
	};
};

class TTask_mateInfo:public TTask_atlas{ -> pileup
public:
	TTask_mateInfo(){ _explanation = "Writing read information per site"; };

	void run(){
		TGenomeWindows genome;
		genome.printMateInformationPerSite();
	};
};

class TTask_writeDepthPerSite:public TTask_atlas{ -> pileup
public:
	TTask_writeDepthPerSite(){ _explanation = "Writing sequencing depth for each site"; };

	void run(){
		TGenomeWindows genome;
		genome.writeDepthPerSite();
	};
};


class TTask_BQSR:public TTask_atlas{ -> suppress?
public:
	TTask_BQSR(){
		_explanation = "Estimating BQSR error re-calibration parameters";
		_citations.push_back("Hofmanova et al. (2016) PNAS");
	};

	void run(){
		TGenomeWindows genome;
		genome.BQSR();
	};
};

class TTask_binQualityScores:public TTask_atlas{ -> a feature of all BAM writing functions, e.g. BamFilter
public:
	TTask_binQualityScores(){ _explanation = "Binning quality scores"; };

	void run(){
		TGenomeWindows genome;
		genome.binQualityScores();
	};
};


class TTask_recalBAM:public TTask_atlas{ -> feature of all BAM writing functions that parse, e.g. BamFilter
public:
	TTask_recalBAM(){ _explanation = "Recalibrating quality scores in a BAM file"; };

	void run(){
		TGenomeWindows genome;
		genome.recalibrateBamFile();
	};
};

class TTask_testBED:public TTask_atlas{ ->should be a proper integration test
public:
	TTask_testBED(){ _explanation = "Testing BED files"; };

	void run(){
		TBed bed;
		bed.test();
	};
};


*/
