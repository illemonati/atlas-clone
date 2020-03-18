/*
 * TTask.h
 *
 *  Created on: Mar 31, 2019
 *      Author: phaentu
 */

#ifndef ALLTASKS_H_
#define ALLTASKS_H_

#include "TTask.h"

//---------------------------------------------------------------------------
// all includes necessary for the application
//---------------------------------------------------------------------------
#include "TGenome.h"
#include "TDistanceEstimator.h"
#include "Simulations/TSimulator.h"
#include "TVcfDiagnostics.h"
#include "TVcfConverter.h"
#include "TMajorMinor.h"
#include "TAlleleCountEstimator.h"
#include "TAlleleFrequencyEstimator.h"
#include "TInbreedingEstimator.h"
#include "TVCFCompare.h"
#include "TPolymorhicWindowIdentifier.h"

//---------------------------------------------------------------------------
// TTask class specific to this application (optional)
//---------------------------------------------------------------------------
class TTask_atlas:public TTask{
public:
	TTask_atlas(){ _citations.push_back("Link et al. (2017) bioRxiv"); };
};

//---------------------------------------------------------------------------
// Function to create map of tasks (adjust in allTasks.cpp file)
//---------------------------------------------------------------------------
void fillTaskMaps(std::map< std::string, TTask* > & taskMap_regular, std::map< std::string, TTask* > & taskMap_debug);

//---------------------------------------------------------------------------
// Create a class for each task, as shown in the example below
//---------------------------------------------------------------------------
//	class TTask_NAME:public TTask{ //or derive from class specific task
//		public:
//		TTask_NAME(){ _explanation = "SAY WHAT THIS TASK IS DOING"; };

//		void run(TParameters & parameters, TLog* logfile){
//			SPECIFY HOW TO EXECUTE THIS TASK
//		};
//	};
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
// Read groups
//---------------------------------------------------------------------------
class TTask_splitRGbyLength:public TTask_atlas{
public:
	TTask_splitRGbyLength(){ _explanation = "Splitting single end read groups in a BAM file"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters, randomGenerator);
		genome.splitSingleEndReadGroups(parameters);
	};
};

class TTask_mergeReadGroups:public TTask_atlas{
public:
	TTask_mergeReadGroups(){ _explanation = "Merging read groups in a BAM file"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters, randomGenerator);
		genome.mergeReadGroups(parameters);
	};
};

//---------------------------------------------------------------------------
// BAM file tools
//---------------------------------------------------------------------------
class TTask_pileup:public TTask_atlas{
public:
	TTask_pileup(){ _explanation = "Printing pileup from BAM file"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters, randomGenerator);
		genome.printPileup(parameters);
	};
};

class TTask_BAMDiagnostics:public TTask_atlas{
public:
	TTask_BAMDiagnostics(){ _explanation = "Estimating approximate depth, read length frequencies and mapping quality frequencies"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters, randomGenerator);
		genome.diagnoseBamFile(parameters);
	};
};

class TTask_assessReadOverlap:public TTask_atlas{
public:
	TTask_assessReadOverlap(){ _explanation = "Estimating distribution of overlap of paired reads in BAM file"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters, randomGenerator);
		genome.assessOverlap(parameters);
	};
};

class TTask_mergeReads:public TTask_atlas{
public:
	TTask_mergeReads(){ _explanation = "Merging paired-end reads in BAM file"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters, randomGenerator);
		genome.mergePairedEndReads(parameters);
	};
};

class TTask_splitMerge:public TTask_atlas{
public:
	TTask_splitMerge(){ _explanation = "Splitting single-end reads and merging paired-end reads and in BAM file"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters, randomGenerator);
		genome.splitMerge(parameters);
	};
};

class TTask_assessDuplication:public TTask_atlas{
public:
	TTask_assessDuplication(){ _explanation = "Quantifying read duplication"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters, randomGenerator);
		genome.estimateDuplicationCounts(parameters);
	};
};

class TTask_mateInfo:public TTask_atlas{
public:
	TTask_mateInfo(){ _explanation = "Writing read information per site"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters, randomGenerator);
		genome.printMateInformationPerSite(parameters);
	};
};

class TTask_assessSoftClipping:public TTask_atlas{
public:
	TTask_assessSoftClipping(){ _explanation = "Assessing level of soft clipping in BAM file"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters, randomGenerator);
		genome.assessSoftClipping(parameters);
	};
};

class TTask_removeSoftClippedBasesFromReads:public TTask_atlas{
public:
	TTask_removeSoftClippedBasesFromReads(){ _explanation = "Removing soft clipped bases from reads"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters, randomGenerator);
		genome.removeSoftClippedBasesFromReads(parameters);
	};
};

class TTask_contextStats:public TTask_atlas{
public:
	TTask_contextStats(){ _explanation = "Writing context statistics to file"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters, randomGenerator);
		genome.contextStats(parameters);
	};
};

class TTask_filterBAM:public TTask_atlas{
public:
	TTask_filterBAM(){ _explanation = "Writing reads that pass filters to BAM file"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters, randomGenerator);
		genome.filterBAM(parameters);
	};
};

//---------------------------------------------------------------------------
// Depth
//---------------------------------------------------------------------------
class TTask_writeDepthPerWindow:public TTask_atlas{
public:
	TTask_writeDepthPerWindow(){ _explanation = "Estimating depth per window"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters, randomGenerator);
		genome.estimateApproximateDepthPerWindow(parameters);
	};
};

class TTask_depthPerSiteDist:public TTask_atlas{
public:
	TTask_depthPerSiteDist(){ _explanation = "Estimating depth per site distribution"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters, randomGenerator);
		genome.estimateDepthPerSite(parameters);
	};
};

class TTask_writeDepthPerSite:public TTask_atlas{
public:
	TTask_writeDepthPerSite(){ _explanation = "Writing sequencing depth for each site"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters, randomGenerator);
		genome.writeDepthPerSite(parameters);
	};
};

class TTask_createDepthMask:public TTask_atlas{
public:
	TTask_createDepthMask(){ _explanation = "Creating depth mask BED file"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters, randomGenerator);
		genome.createDepthMask(parameters);
	};
};

class TTask_createNonRefMask:public TTask_atlas{
public:
	TTask_createNonRefMask(){ _explanation = "Creating mask BED file with sites with non-reference alleles"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters, randomGenerator);
		genome.writeNonConservedBed(parameters);
	};
};

class TTask_allelicDepth:public TTask_atlas{
public:
	TTask_allelicDepth(){ _explanation = "Estimating allelic depth distribution"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters, randomGenerator);
		genome.allelicDepth(parameters);
	};
};

class TTask_downsample:public TTask_atlas{
public:
	TTask_downsample(){ _explanation = "Downsampling a BAM file by removing reads"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters, randomGenerator);
		genome.downSampleBamFile(parameters);
	};
};

class TTask_separateReads:public TTask_atlas{
public:
	TTask_separateReads(){ _explanation = "Separating reads into different BAM files"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters, randomGenerator);
		genome.separateReads(parameters);
	};
};

class TTask_downSampleReads:public TTask_atlas{
public:
	TTask_downSampleReads(){ _explanation = "Downsampling a BAM file by setting bases to N"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters, randomGenerator);
		genome.downSampleReads(parameters);
	};
};

//---------------------------------------------------------------------------
// PMD
//---------------------------------------------------------------------------
class TTask_estimatePMD:public TTask_atlas{
public:
	TTask_estimatePMD(){ _explanation = "Estimating Post-Mortem Damage (PMD) patterns"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters, randomGenerator);
		genome.estimatePMD(parameters);
	};
};

class TTask_PMDS:public TTask_atlas{
public:
	TTask_PMDS(){
		_explanation = "Filtering for ancient reads using PMDS";
		_citations.push_back("Skoglund et al. (2014) PNAS");
	};

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters, randomGenerator);
		genome.runPMDS(parameters);
	};
};

class TTask_PSMC:public TTask_atlas{
public:
	TTask_PSMC(){ _explanation = "Generating a PSMC Input file probabilistically"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters, randomGenerator);
		genome.generatePSMCInput(parameters);
	};
};

//---------------------------------------------------------------------------
// Recalibration
//---------------------------------------------------------------------------
class TTask_recal:public TTask_atlas{
public:
	TTask_recal(){
		_explanation = "Estimating error re-calibration parameters";
		_citations.push_back("Kousathanas et al. (2017) Genetics");
	};

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters, randomGenerator);
		genome.estimateErrorCalibrationEM(parameters);
	};
};

class TTask_recalLL:public TTask_atlas{
public:
	TTask_recalLL(){ _explanation = "Calculating LL for error re-calibration function"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters, randomGenerator);
		genome.calculateLikelihoodErrorCalibrationEM(parameters);
	};
};

class TTask_BQSR:public TTask_atlas{
public:
	TTask_BQSR(){
		_explanation = "Estimating BQSR error re-calibration parameters";
		_citations.push_back("Hofmanova et al. (2016) PNAS");
	};

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters, randomGenerator);
		genome.BQSR(parameters);
	};
};

class TTask_recalBAM:public TTask_atlas{
public:
	TTask_recalBAM(){ _explanation = "Recalibrating quality scores in a BAM file"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters, randomGenerator);
		genome.recalibrateBamFile(parameters);
	};
};

class TTask_qualityDist:public TTask_atlas{
public:
	TTask_qualityDist(){ _explanation = "Printing Quality Distribution"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters, randomGenerator);
		genome.printQualityDistribution(parameters);
	};
};

class TTask_qualityTransformation:public TTask_atlas{
public:
	TTask_qualityTransformation(){ _explanation = "Printing Quality Transformation"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters, randomGenerator);
		genome.printQualityTransformation(parameters);
	};
};

class TTask_binQualityScores:public TTask_atlas{
public:
	TTask_binQualityScores(){ _explanation = "Binning quality scores"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters, randomGenerator);
		genome.binQualityScores(parameters);
	};
};

//---------------------------------------------------------------------------
// Caller
//---------------------------------------------------------------------------
class TTask_call:public TTask_atlas{
public:
	TTask_call(){ _explanation = "Calling genotypes"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters, randomGenerator);
		genome.callGenotypes(parameters);
	};
};

//---------------------------------------------------------------------------
// Diversity (theta)
//---------------------------------------------------------------------------
class TTask_estimateTheta:public TTask_atlas{
public:
	TTask_estimateTheta(){
		_explanation = "Estimating heterozygosity (theta)";
		_citations.push_back("Kousathanas et al. (2017) Genetics");
	};

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters, randomGenerator);
		genome.estimateTheta(parameters);
	};
};

class TTask_estimateThetaRatio:public TTask_atlas{
public:
	TTask_estimateThetaRatio(){ _explanation = "Estimate the ratio in heterozygosity (theta) between genomic regions"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters, randomGenerator);
		genome.estimateThetaRatio(parameters);
	};
};

class TTask_thetaLLSurface:public TTask_atlas{
public:
	TTask_thetaLLSurface(){ _explanation = "Calculating the theta LL surface for each window"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters, randomGenerator);
		genome.calcThetaLikelihoodSurfaces(parameters);
	};
};

class TTask_downsamplingThetaQC:public TTask_atlas{
public:
	TTask_downsamplingThetaQC(){ _explanation = "QC recalibration by estimating theta on downsampled data for each window"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters, randomGenerator);
		genome.performDownsamplingThetaQC(parameters);
	};
};

//---------------------------------------------------------------------------
// Population tools
//---------------------------------------------------------------------------
class TTask_GLF:public TTask_atlas{
public:
	TTask_GLF(){ _explanation = "Writing genotype likelihoods to a GLF file"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters, randomGenerator);
		genome.writeGLF(parameters);
	}
};

class TTask_printGLF:public TTask_atlas{
public:
	TTask_printGLF(){ _explanation = "Printing a GLF file to screen"; };

	void run(TParameters & parameters, TLog* logfile){
		std::string glf = parameters.getParameterString("glf");
		TGlfReader reader(glf);
		reader.printToEnd();
	};
};

class TTask_majorMinor:public TTask_atlas{
public:
	TTask_majorMinor(){
		_citations.push_back("Skotte et al. (2012) Genetic Epidemiology");
		_explanation = "Estimating major and minor alles"; };

	void run(TParameters & parameters, TLog* logfile){
		TMajorMinor majorMinor(logfile, parameters, randomGenerator);
		majorMinor.estimateMajorMinor(parameters);
	};
};

class TTask_estimateDist:public TTask_atlas{
public:
	TTask_estimateDist(){ _explanation = "Estimating the genetic distance between individuals"; };

	void run(TParameters & parameters, TLog* logfile){
		TDistanceEstimator distEst(logfile, parameters);
		distEst.estimateDistances(parameters);
	};
};

class TTask_estimateAlleleCounts:public TTask_atlas{
public:
	TTask_estimateAlleleCounts(){ _explanation = "Estimating population allele counts"; };

	void run(TParameters & parameters, TLog* logfile){
		TAlleleCountEstimator alleleCountEst(parameters, logfile);
		alleleCountEst.estimateAlleleCounts(parameters, randomGenerator);
	};
};

class TTask_writeAlleleFrequencyLikelihoods:public TTask_atlas{
public:
	TTask_writeAlleleFrequencyLikelihoods(){ _explanation = "Writing allele frequency likelihoods"; };

	void run(TParameters & parameters, TLog* logfile){
		TAlleleCountEstimator alleleCountEst(parameters, logfile);
		alleleCountEst.writeAlleleFrequencyLikelihoods(parameters);
	};
};

class TTask_transformAlleleCountFormat:public TTask_atlas{
public:
	TTask_transformAlleleCountFormat(){ _explanation = "Transforming allele counts file format"; };

	void run(TParameters & parameters, TLog* logfile){
		TAlleleCountEstimator alleleCountEst(parameters, logfile);
		alleleCountEst.transformFormat(parameters);
	};
};

class TTask_estimateAlleleFreq:public TTask_atlas{
public:
	TTask_estimateAlleleFreq(){ _explanation = "Estimating population allele frequencies"; };

	void run(TParameters & parameters, TLog* logfile){
		TAlleleFreqEstimator alleleFreqEstimator(parameters, logfile);
		alleleFreqEstimator.estimateAlleleFreq(parameters, randomGenerator);
	};
};

class TTask_estimateInbreeding:public TTask_atlas{
public:
	TTask_estimateInbreeding(){ _explanation = "Estimating the inbreeding coefficient"; };

	void run(TParameters & parameters, TLog* logfile){
		TInbreedingEstimator inbreedingEstimator(parameters, logfile);
		inbreedingEstimator.runEstimation(parameters);
	};
};

class TTask_inbreedingLikelihood:public TTask_atlas{
public:
	TTask_inbreedingLikelihood(){ _explanation = "Estimating likelihood surfaces for the inbreeding model"; };

	void run(TParameters & parameters, TLog* logfile){
		TInbreedingEstimator inbreedingEstimator(parameters, logfile);
			if(parameters.parameterExists("llGamma"))
				inbreedingEstimator.writeLikelihoodForDebuggingGamma(parameters);
	//			if(parameters->parameterExists("llBeta"))
	//				inbreedingEstimator.writeLikelihoodForDebuggingBeta(parameters);
			else if(parameters.parameterExists("llP"))
				inbreedingEstimator.writeLikelihoodForDebuggingAlleleFreq(parameters);
			else if(parameters.parameterExists("llF"))
				inbreedingEstimator.writeLikelihoodForDebuggingF(parameters);
			else
				throw "define parameter for which to calculate likelihood surface!";
	};
};

class TTask_identifyPolymorphicWindows:public TTask_atlas{
public:
	TTask_identifyPolymorphicWindows(){ _explanation = "Identifying windows for which samples are polymoprhic"; };

	void run(TParameters & parameters, TLog* logfile){
		TPolymorhicWindowIdentifier identifier(parameters, logfile);
		identifier.identifyPolymorphicWindows(parameters, randomGenerator);
	}
};
//---------------------------------------------------------------------------
// Simulations
//---------------------------------------------------------------------------
class TTask_simulate:public TTask_atlas{
public:
	TTask_simulate(){
		_explanation = "Generating simulations";
	};

	void run(TParameters & parameters, TLog* logfile){
		//initialize simulator
		TSimulator* simulator;
		std::string method = parameters.getParameterStringWithDefault("type", "one");
		if(method == "one")
			simulator = new TSimulatorOneIndividual(logfile, parameters, randomGenerator);
		else if(method == "pair")
			simulator = new TSimulatorPairOfIndividuals(logfile, parameters, randomGenerator);
		else if(method == "SFS")
			simulator = new TSimulatorSFS(logfile, parameters, randomGenerator);
		else if(method == "HW")
			simulator = new TSimulatorHardyWeinberg(logfile, parameters, randomGenerator);
		else throw "Unknown simulation method '" + method + "'!";

		//now run simulations
		simulator->runSimulations();

		//clean up
		delete simulator;
	};
};

//---------------------------------------------------------------------------
// VCF tools
//---------------------------------------------------------------------------
class TTask_VCFDiagnostics:public TTask_atlas{
public:
	TTask_VCFDiagnostics(){ _explanation = "Diagnosing a VCF file"; };

	void run(TParameters & parameters, TLog* logfile){
		VcfDiagnostics VcfDiagnoser(parameters, logfile);
		VcfDiagnoser.assessAllelicImbalance(parameters);
	};
};

class TTask_VCFToInvariantBed:public TTask_atlas{
public:
	TTask_VCFToInvariantBed(){ _explanation = "Writing a BED file from invariant sites in a VCF file"; };

	void run(TParameters & parameters, TLog* logfile){
		VcfDiagnostics VcfDiagnoser(parameters, logfile);
		VcfDiagnoser.vcfToInvariantBed();
	};
};

class TTask_VCFToBeagle:public TTask_atlas{
public:
	TTask_VCFToBeagle(){ _explanation = "Converting a VCF file to Beagle format"; };

	void run(TParameters & parameters, TLog* logfile){
		TVcfToBeagle VcfToBeagle(parameters, logfile);
        VcfToBeagle.vcfToBeagle(parameters);
	};
};

class TTask_VCFToLFMM:public TTask_atlas{
public:
    TTask_VCFToLFMM(){ _explanation = "Converting a VCF file to LFMM format"; };

    void run(TParameters & parameters, TLog* logfile){
        std::string output = parameters.getParameterString("geno");
        if (output == "postGeno"){
            TVcfToLFMMPostGeno vcfToLFMMPostGeno(parameters, logfile);
            vcfToLFMMPostGeno.vcfToLFMM(parameters);
        }
        else if (output == "calledGeno"){
            TVcfToLFMMCalledGeno vcfToLFMMCalledGeno(parameters, logfile);
            vcfToLFMMCalledGeno.vcfToLFMM(parameters);
        }
        else {
            throw std::runtime_error("Unknown LFMM format '" + output + "'! Please choose either 'postGeno' or 'calledGeno'.");
        }
    };
};

class TTask_VCFFixInt:public TTask_atlas{
public:
	TTask_VCFFixInt(){ _explanation = "Fixing integers printed as floats in VCF file"; };

	void run(TParameters & parameters, TLog* logfile){
		VcfDiagnostics VcfDiagnoser(parameters, logfile);
		VcfDiagnoser.fixIntAsFloat();
	};
};

class TTask_VCFCompare:public TTask_atlas{
public:
	TTask_VCFCompare(){ _explanation = "Comparing genotype calls in two VCF files"; };

	void run(TParameters & parameters, TLog* logfile){
		TVCFCompare compare(logfile);
		compare.compareVCFFiles(parameters);
	};
};

/*
class TTask_filterVCF:public TTask_atlas{
public:
	TTask_filterVCF(){ _explanation = "Filtering a VCF file"; };

	void run(TParameters & parameters, TLog* logfile){
		TVcfFilter vcfFilter(parameters, logfile);
		vcfFilter.filterVCF(parameters);
	};
};
*/


#endif /* ALLTASKS_H_ */
