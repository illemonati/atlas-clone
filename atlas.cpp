/*
 * atlas.cpp
 */

#include "TDepthCalculator.h"
#include "TEstimateHKY85.h"
#include "TFromTo.h"
#include "TGLFPrinter.h"
#include "TIlluminaIdentifier.h"
#include "TOverlappingReadsMerger.h"
#include "TPairAnalyser.h"
#include "TSafEstimator.h"
#include "TTransitionTabler.h"
#include "coretools/Main/TMain.h"

//BAM
#include "TFilterBam.h"
#include "TBamDiagnoser.h"
#include "TPMDSCalculator.h"
#include "TPileup.h"
#include "TQualityDistribution.h"
#include "TReadGroupInfo.h"
#include "TReadGroupMerger.h"
#include "TSoftClipping.h"

//window
#include "TAlleleCountEstimator.h"
#include "TAlleleFrequencyEstimator.h"
#include "TAllelicDepthCounts.h"
#include "TAncestralAlleleEstimator.h"
#include "TBamDownsampler.h"
#include "TCaller.h"
#include "TCreateBedMask.h"
#include "TDistanceEstimator.h"
#include "TErrorEstimator.h"
#include "TEstimateMutationLoad.h"
#include "TEstimateTheta.h"
#include "TF2Estimator.h"
#include "THardyWeinbergTest.h"
#include "TInbreedingEstimator.h"
#include "TMajorMinor.h"
#include "TPSMCInput.h"
#include "TPolymorphicWindowIdentifier.h"
#include "TWriteGLF.h"
#include "TSpearmanGWAS.h"

//VCF
#include "TVcfCompare.h"
#include "TVcfConverter.h"
#include "TVcfDiagnostics.h"

//simulations
#include "TSimulator.h"

//other
#include "TPositionBasedLiftOver.h"

void addTaks(coretools::TMain & main) {
	//BAM
	{
	constexpr auto groupName = "Read";
	main.createGroupedTask<GenomeTasks::TFilterBam>(groupName, "filterBAM", "Writing reads that pass filters to BAM file");
	main.createGroupedTask<GenomeTasks::TOverlappingReadsMerger>(groupName, "mergeOverlappingReads", "Merging paired-end reads in BAM file");
	main.createGroupedTask<GenomeTasks::TReadGroupMerger>(groupName, "mergeRG", "Merging read groups in a BAM file");	
	main.createGroupedTask<GenomeTasks::TBamDiagnoser>(groupName, "BAMDiagnostics", "Estimating depth and read property frequencies");
	main.createGroupedTask<GenomeTasks::TAssessSoftClipping>(groupName, "assessSoftClipping", "Assessing level of soft clipping in BAM file");
	main.createGroupedTask<GenomeTasks::TSoftClipsTrimmer>(groupName, "trimSoftClips", "Removing soft clipped bases from reads");
	main.createGroupedTask<GenomeTasks::TQualityTransformation>(groupName, "qualityTransformation", "Printing Quality Transformation");
	main.createGroupedTask<GenomeTasks::TBamDownsampler>(groupName, "downsample", "Downsampling a BAM file");
	main.createGroupedTask<GenomeTasks::TPMDSCalculator>(groupName, "PMDS", "Filtering for ancient reads using PMDS", "Skoglund et al. (2014) PNAS");	
	}

	//window tasks
	{
	constexpr auto groupName = "Site";
	main.createGroupedTask<GenotypeLikelihoods::TErrorEstimator>(groupName, "estimateErrors", "Estimating PMD pattern and Sequencing Errors", "Kousathanas et al. (2017) Genetics");
	main.createGroupedTask<GenomeTasks::TMaskCreator>(groupName, "createMask", "Creating a mask BED file");
	main.createGroupedTask<GenomeTasks::TAllelicDepth>(groupName, "allelicDepth", "Writing genotype likelihoods to a GLF file");
	main.createGroupedTask<GenomeTasks::TPSMCInput>(groupName, "PSMC", "Generating a PSMC Input file probabilistically");
	main.createGroupedTask<GenomeTasks::TCall>(groupName, "call", "Calling genotypes");
	main.createGroupedTask<GenomeTasks::TEstimateTheta>(groupName, "theta", "Estimating heterozygosity (theta)", "Kousathanas et al. (2017) Genetics");
	main.createGroupedTask<GenomeTasks::TEstimateThetaRatio>(groupName, "thetaRatio", "Estimate ratio in heterozygosity (theta) between genomic regions", "Kousathanas et al. (2017) Genetics");
	main.createGroupedTask<GenomeTasks::TWriteGLF>(groupName, "GLF", "Writing genotype likelihoods to a GLF file");	
	main.createGroupedTask<GenomeTasks::TEstimateMutationLoad>(groupName, "mutationLoad", "Estimating mutation load across the genome");
	main.createGroupedTask<GenomeTasks::TEstimateHKY85>(groupName, "HKY85", "Estimating HKY85 genotype Distribution");
	main.createGroupedTask<GenomeTasks::TPileup>(groupName, "pileup", "Printing pileup from BAM file");
	}

	//Population tools
	{
	constexpr auto groupName = "Population";
	main.createGroupedTask<PopulationTools::TGLFPrinter>(groupName, "printGLF", "Printing a GLF file to screen");
	main.createGroupedTask<PopulationTools::TMajorMinor>(groupName, "majorMinor", "Estimating major and minor alles", "Skotte et al. (2012) Genetic Epidemiology");
	main.createGroupedTask<PopulationTools::TDistanceEstimator>(groupName, "geneticDist", "Estimating the genetic distance between individuals");
	main.createGroupedTask<PopulationTools::TAlleleCounter>(groupName, "alleleCounts", "Estimating population allele counts");
	main.createGroupedTask<PopulationTools::TAlleleFreqEstimator>(groupName, "alleleFreq", "Estimating population allele frequencies");
	main.createGroupedTask<PopulationTools::TInbreedingEstimator>(groupName, "inbreeding", "Estimating the inbreeding coefficient", "Burger et al. (2020) Current Biology");
	main.createGroupedTask<PopulationTools::TPolymorphicWindowIdentifier>(groupName, "polymorphicWindows", "Identifying windows for which samples are polymorphic");
	main.createGroupedTask<PopulationTools::TF2Estimator>(groupName, "calculateF2", "Calculate F2 between samples, and within/between populations");
	main.createGroupedTask<PopulationTools::TAncestralAlleleEstimator>(groupName, "ancestralAlleles", "Writing FASTA-file with ancestral alleles");
	main.createGroupedTask<PopulationTools::TSafEstimator>(groupName, "saf", "Estimating Site Allele Frequencies");	
	}

	//VCF
	{
	constexpr auto groupName = "VCF";
	main.createGroupedTask<VCF::TVcfDiagnostics>(groupName, "VCFDiagnostics", "Diagnosing a VCF file");
	main.createGroupedTask<VCF::TVCFConverter>(groupName, "convertVCF", "Converting a VCF file to other formats");
	main.createGroupedTask<VCF::TVcfCompare>(groupName, "VCFCompare", "Comparing genotype calls in two VCF files");
	main.createGroupedTask<PopulationTools::THardyWeinbergTest>(groupName, "testHardyWeinberg", "Testing for Hardy-Weinberg equilibrium across multiple populations");
	}

	//simulations
	{
	constexpr auto groupName = "Simulation";
	main.createGroupedTask<Simulations::TSimulationRunner>(groupName, "simulate", "Simulate bam- or vcf-file[s]");
	}

	/*other
	{
	constexpr auto groupName = "Other";
	main.createGroupedTask<PopulationTools::TPositionBasedLiftOver>(groupName, "liftOver", "Position-based lift over from one reference to another");
	}
	*/

	// Debug tasks
	main.createDebugTask<GenomeTasks::TEstimateThetaLLSurface>("thetaLLSurface", "Calculating the theta LL surface for each window");
	main.createDebugTask<BAM::RGInfo::TReadGroupInfoTest>("json", "Testing JSON stuff");
	main.createDebugTask<GenomeTasks::TFromTo>("fromTo", "FromTo");
	main.createDebugTask<PopulationTools::TSpearmanGWAS>("SpearmanGWAS", "GWAS with Spearman correlation, allowing for population-specific signs");
	main.createDebugTask<GenomeTasks::TIlluminaIdentifier>("identifyIllumina", "Reassigning read groups based on the platform unit in their name");
	main.createDebugTask<GenomeTasks::TTransitionTabler>("transitionTable", "Create transition table from sequencer-start");
	main.createDebugTask<GenomeTasks::TPairAnalyser>("analysePairs", "analyse Pairs");
	main.createDebugTask<PopulationTools::TPositionBasedLiftOver>("liftOver", "Position-based lift over from one reference to another");
	main.createDebugTask<GenomeTasks::TDepthCalculator>("averageDepth", "Calculate average depth of BAM file");	

};

void addTests(coretools::TMain & ){
	// Use testing.addTest to add a single test
};

//---------------------------------------------------------------------------
//Main function
//---------------------------------------------------------------------------
int main(int argc, char* argv[]){
	coretools::TMain main("ATLAS", "2.0.0-rc.10 (Release Candidate)", "https://bitbucket.org/wegmannlab/atlas", "andreas.fueglistaler@unifr.ch");

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

*/
