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
#include "Simulations/runSimulations.h"
#include "TVcfDiagnostics.h"
#include "TMajorMinor.h"
#include "TAlleleCountEstimator.h"
#include "TInbreedingEstimator.h"

//---------------------------------------------------------------------------
// TTask class specific to this application (optional)
//---------------------------------------------------------------------------
class TTask_atlas:public TTask{
public:
	TTask_atlas(){ _citations.push_back("Link et al. (2019) SOMEWHERE"); };
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
		TGenome genome(logfile, parameters);
		genome.splitSingleEndReadGroups(parameters);
	};
};

class TTask_mergeReadGroups:public TTask_atlas{
public:
	TTask_mergeReadGroups(){ _explanation = "Merging read groups in a BAM file"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters);
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
		TGenome genome(logfile, parameters);
		genome.printPileup(parameters);
	};
};

class TTask_BAMDiagnostics:public TTask_atlas{
public:
	TTask_BAMDiagnostics(){ _explanation = "Estimating approximate depth, read length frequencies and mapping quality frequencies"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters);
		genome.diagnoseBamFile(parameters);
	};
};

class TTask_assessReadOverlap:public TTask_atlas{
public:
	TTask_assessReadOverlap(){ _explanation = "Estimating distribution of overlap of paired reads in BAM file"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters);
		genome.assessOverlap(parameters);
	};
};

class TTask_mergeReads:public TTask_atlas{
public:
	TTask_mergeReads(){ _explanation = "Merging paired-end reads in BAM file"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters);
		genome.mergePairedEndReadsNoOrder(parameters);
	};
};

class TTask_assessDuplication:public TTask_atlas{
public:
	TTask_assessDuplication(){ _explanation = "Quantifying read duplication"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters);
		genome.estimateDuplicationCounts(parameters);
	};
};

class TTask_mateInfo:public TTask_atlas{
public:
	TTask_mateInfo(){ _explanation = "Writing mate information per site"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters);
		genome.printMateInformationPerSite(parameters);
	};
};

class TTask_assessSoftClipping:public TTask_atlas{
public:
	TTask_assessSoftClipping(){ _explanation = "Assessing level of soft clipping in BAM file"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters);
		genome.assessSoftClipping(parameters);
	};
};

class TTask_removeSoftClippedBasesFromReads:public TTask_atlas{
public:
	TTask_removeSoftClippedBasesFromReads(){ _explanation = "Removing soft clipped bases from reads"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters);
		genome.removeSoftClippedBasesFromReads(parameters);
	};
};

//---------------------------------------------------------------------------
// Depth
//---------------------------------------------------------------------------
class TTask_writeDepthPerWindow:public TTask_atlas{
public:
	TTask_writeDepthPerWindow(){ _explanation = "Estimating depth per window"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters);
		genome.estimateApproximateDepthPerWindow(parameters);
	};
};

class TTask_depthPerSiteDist:public TTask_atlas{
public:
	TTask_depthPerSiteDist(){ _explanation = "Estimating depth per site"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters);
		genome.estimateDepthPerSite(parameters);
	};
};

class TTask_writeDepthPerSite:public TTask_atlas{
public:
	TTask_writeDepthPerSite(){ _explanation = "Writing sequencing depth for each site"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters);
		genome.writeDepthPerSite(parameters);
	};
};

class TTask_createDepthMask:public TTask_atlas{
public:
	TTask_createDepthMask(){ _explanation = "Creating depth mask BED file"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters);
		genome.createDepthMask(parameters);
	};
};

class TTask_allelicDepth:public TTask_atlas{
public:
	TTask_allelicDepth(){ _explanation = "Estimating allelic depth distribution"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters);
		genome.allelicDepth(parameters);
	};
};

class TTask_downsample:public TTask_atlas{
public:
	TTask_downsample(){ _explanation = "Downsampling a BAM file by removing reads"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters);
		genome.downSampleBamFile(parameters);
	};
};

class TTask_downSampleReads:public TTask_atlas{
public:
	TTask_downSampleReads(){ _explanation = "Downsampling a BAM file by setting bases to N"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters);
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
		TGenome genome(logfile, parameters);
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
		TGenome genome(logfile, parameters);
		genome.runPMDS(parameters);
	};
};

class TTask_PSMC:public TTask_atlas{
public:
	TTask_PSMC(){ _explanation = "Generating a PSMC Input file probabilistically"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters);
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
		TGenome genome(logfile, parameters);
		genome.estimateErrorCalibrationEM(parameters);
	};
};

class TTask_recalLL:public TTask_atlas{
public:
	TTask_recalLL(){ _explanation = "Calculating LL for error re-calibration function"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters);
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
		TGenome genome(logfile, parameters);
		genome.BQSR(parameters);
	};
};

class TTask_recalBAM:public TTask_atlas{
public:
	TTask_recalBAM(){ _explanation = "Recalibrating quality scores in a BAM file"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters);
		genome.recalibrateBamFile(parameters);
	};
};

class TTask_qualityDist:public TTask_atlas{
public:
	TTask_qualityDist(){ _explanation = "Printing Quality Distribution"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters);
		genome.printQualityDistribution(parameters);
	};
};

class TTask_qualityTransformation:public TTask_atlas{
public:
	TTask_qualityTransformation(){ _explanation = "Printing Quality Transformation"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters);
		genome.printQualityTransformation(parameters);
	};
};

class TTask_binQualityScores:public TTask_atlas{
public:
	TTask_binQualityScores(){ _explanation = "Binning quality scores"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters);
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
		TGenome genome(logfile, parameters);
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
		TGenome genome(logfile, parameters);
		genome.estimateTheta(parameters);
	};
};

class TTask_estimateThetaRatio:public TTask_atlas{
public:
	TTask_estimateThetaRatio(){ _explanation = "Estimate the ratio in heterozygosity (theta) between genomic regions"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters);
		genome.estimateThetaRatio(parameters);
	};
};

class TTask_thetaLLSurface:public TTask_atlas{
public:
	TTask_thetaLLSurface(){ _explanation = "Calculating the theta LL surface for each window"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters);
		genome.calcThetaLikelihoodSurfaces(parameters);
	};
};


//---------------------------------------------------------------------------
// Population tools
//---------------------------------------------------------------------------
class TTask_GLF:public TTask_atlas{
public:
	TTask_GLF(){ _explanation = "Writing genotype likelihoods to a GLF file"; };

	void run(TParameters & parameters, TLog* logfile){
		TGenome genome(logfile, parameters);
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
	TTask_majorMinor(){ _explanation = "Estimating major and minor alles"; };

	void run(TParameters & parameters, TLog* logfile){
		TMajorMinor majorMinor(parameters, logfile);
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
		alleleCountEst.estimateAlleleCounts(parameters);
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

class TTask_estimateAlleleFreq:public TTask_atlas{
public:
	TTask_estimateAlleleFreq(){ _explanation = "Estimating population allele frequencies"; };

	void run(TParameters & parameters, TLog* logfile){
		TAlleleFreqEstimator alleleFreqEstimator(parameters, logfile);
		alleleFreqEstimator.estimateAlleleFreq(parameters);
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

//---------------------------------------------------------------------------
// Simulations
//---------------------------------------------------------------------------
class TTask_simulate:public TTask_atlas{
public:
	TTask_simulate(){
		_explanation = "Generating simulations";
	};

	void run(TParameters & parameters, TLog* logfile){
		runSimulations(parameters, logfile);
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
	TTask_VCFToInvariantBed(){ _explanation = "Diagnosing a VCF file"; };

	void run(TParameters & parameters, TLog* logfile){
		VcfDiagnostics VcfDiagnoser(parameters, logfile);
		VcfDiagnoser.vcfToInvariantBed();
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
		throw "run not defined!";
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
