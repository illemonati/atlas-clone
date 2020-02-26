/*
 * allTasks.cpp
 *
 *  Created on: Mar 31, 2019
 *      Author: phaentu
 */

#include "allTasks.h"

void fillTaskMaps(std::map< std::string, TTask* > & taskMap_regular, std::map< std::string, TTask* > & taskMap_debug){
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
	taskMap_regular["allelicDepth"] = new TTask_allelicDepth();
	taskMap_regular["downsample"] = new TTask_downsample();
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
};


