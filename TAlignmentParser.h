/*
 * TAlignmentParser.h
 *
 *  Created on: Nov 14, 2017
 *      Author: wegmannd
 */

#ifndef TALIGNMENTPARSER_H_
#define TALIGNMENTPARSER_H_

#include <TBedReaderWindows.h>
#include "TBamFile.h"
#include "bamtools/api/BamWriter.h"
#include "bamtools/utils/bamtools_fasta.h"
#include "TGenotypeMap.h"
#include "TLog.h"
#include "TPostMortemDamage.h"
#include "TAlignment.h"
#include "TWindow.h"
#include "TBed.h"
#include "TChromosomes.h"
#include "TFastaBuffer.h"
#include <vector>
#include "progressTools.h"
#include "TSiteSubset.h"

//-----------------------------------------------------
//TAlignmentParser
//-----------------------------------------------------
class TAlignmentParser{
protected:
	TLog* logfile;
	BAM::TBamFile& bamFile;
	GenotypeLikelihoods::TGenotypeLikelihoodCalculator& genotypeLikelihoodCalculator;

	//maps
	TGenotypeMap genoMap;
	TQualityMap qualMap;

 	//reference
	bool hasReference;
	bool chrChangedWindow;
	BAM::TFastaBuffer* fastaBuffer;

	//read trimming
	bool trimReads;
	int trimmingLength3Prime;
	int trimmingLength5Prime;

	//filters
	bool applyQualityFilter;
	TQualityFilter qualityFilter;
	bool applyContextFilter;
	std::map<BaseContext,int> ignoreTheseContexts;

	//functions for initialization
	void _setReadTrimming(TParameters & params);
	void _setQualityFilter(TParameters & params);
	void _setContextFilter(TParameters & params);
	void _setQualityRangeForPrinting(TParameters & params);

	//functions for reading
	void _fillAlignment(BAM::TAlignment & alignment);

public:
	TAlignmentParser(TParameters & params, TLog* logfile, BAM::TBamFile & BamFile, GenotypeLikelihoods::TGenotypeLikelihoodCalculator& GenotypeLikelihoodCalculator);


};

//-----------------------------------------------------
//TAlignmentParserWindows
//-----------------------------------------------------
class TAlignmentParserWindows:public TAlignmentParser{
private:
	BAM::TAlignment* oldAlignment;
	bool oldAlignmentInitialized;
	bool oldAlignmentMustBeConsidered;

	//counters
	bool hasWindowIndent;

	//window params
	unsigned int windowSize;
	unsigned int numWindowsOnChr;
	unsigned int windowNumber;
	double maxMissing;
	double maxRefN;

	//masks
	BAM::TBedReaderWindows* mask;

	//filters
	bool applyDepthFilter;
	uint32_t readUpToDepth, minDepth, maxDepth;

	//blacklist
	bool _updateBlacklist;
	BAM::TMateFinder blacklist;

	//limit windows
	long limitWindows;
	int skipWindows;

	//contructor functions
	void _setWindowParameters(TParameters & params);
	void _setWindowFilters(TParameters & params);
	void _setSiteFilters(TParameters & params);



	void _setQualityRangeForPrinting(int minQual, int maxQual);

	void _setMasks(TParameters & params);
//	void initializeSiteSubset(TParameters & params);
	void _setParsingLimits(TParameters & params);
	void _setReadTrimming(TParameters & params);
	void _setChrPloidy(TParameters & params);

	//move genome
	void _jumpToEnd();
	void _restartChromosomes(TWindow_base & window);
	void _moveChromosome(TWindow_base & window);
	bool _moveToNextWindowOnChr(TWindow_base & window);
	bool _moveToNextPredefinedWindow(TWindow_base & window);
	bool _moveWindow(TWindow_base & window);

	GenotypeLikelihoods::PMDType _getEnumPMDType(std::string pmdType);
	void _initializePostMortemDamage(TParameters & params);

	bool _readAlignment();
	bool _applyFilters();
	bool _fillAlignment(BAM::TAlignment & alignment);
	void _readAlignmentsIntoWindow(TWindow & window);
	void _applyWindowFilters(TWindow_base & window);
	void _recalibrate(BAM::TAlignment & alignment);
	void _adaptQualityWhenMerging(TBase & bestBase, TBase & worstBase, const bool & adaptQuality);

public:
	//alignment: goal is to make this private!
	bool doMasking, considerRegions;
	bool windowsPredefined;
	TBed* predefinedWindows;
	bool sitesProvided;
	TSiteSubset* subset;



	//std::string outname;

	//construction
	TAlignmentParser();
	TAlignmentParser(TParameters & params, BAM::TBamFile * BamFile, GenotypeLikelihoods::TGenotypeLikelihoodCalculator* GenotypeLikelihoodCalculator, TLog* Logfile);
	~TAlignmentParser();


	//getters
	unsigned int getWindowSize(){return windowSize;};
	int getMaxPhredInt(){return maxQualityAsPhredInt;};
	uint32_t getMaxDepth(){ return maxDepth; };

	//setters
	void setParsingToTrue(){_parse = true;};

	void addReference(BAM::TFastaBuffer* FastaBuffer);


	//blacklist
	void setUpdateBlacklistToTrue();
	void setWriteBlacklistToFileToTrue(const std::string filename);

	//functions to read and _parse



	//reading data requires windows
	bool readDataInNextWindow(TWindow & window);
	void downsampleWindow(TWindow_base & destination, TWindow & source, const double downsamplingProb, TRandomGenerator* randomGenerator);

	//reading data only requires alignments
	bool readNextAlignment(BAM::TAlignment & alignment); //to be used to go through bam file alignment by alignment

	//qualityTransformation
	//void initializeRecalibrationForQualityTransformation(TParameters & params);

	//merging
	//TODO: why are these functions not in the merger?
	void mergeAlignedBasesBamReads(BAM::TAlignment* fwdAlignment, BAM::TAlignment* revAlignment, bool adaptQuality);
	void mergeAlignedBasesBamReadsRandom(BAM::TAlignment* fwdAlignment, BAM::TAlignment* revAlignment, bool adaptQuality, TRandomGenerator* randomGenerator);
	void mergeAlignedBasesOneRead(BAM::TAlignment* fwdAlignment, BAM::TAlignment* revAlignment, bool adaptQuality, TRandomGenerator* randomGenerator);
};


#endif /* TALIGNMENTPARSER_H_ */
