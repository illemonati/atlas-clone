/*
 * TAlignmentParser.h
 *
 *  Created on: Nov 14, 2017
 *      Author: wegmannd
 */

#ifndef TALIGNMENTPARSER_H_
#define TALIGNMENTPARSER_H_

#include "TBamFile.h"
#include "bamtools/api/BamWriter.h"
#include "bamtools/utils/bamtools_fasta.h"
#include "TGenotypeMap.h"
#include "TLog.h"
#include "TPostMortemDamage.h"
#include "TAlignment.h"
#include "TWindow.h"
#include "TBed.h"
#include "TBedReader.h"
#include "TChromosomes.h"
#include "TFastaBuffer.h"
#include <vector>
#include "progressTools.h"


//-----------------------------------------------------
//TAlignmentParser
//-----------------------------------------------------
class TAlignmentParser{
private:
	TLog* logfile;
	bool _parse;

	TAlignment* oldAlignment;
	bool oldAlignmentInitialized;
	bool oldAlignmentMustBeConsidered;

	//counters
	bool hasWindowIndent;

	//read trimming
	bool trimReads;
	int trimmingLength3Prime;
	int trimmingLength5Prime;

	//window params
	unsigned int windowSize;
	unsigned int numWindowsOnChr;
	unsigned int windowNumber;
	double maxMissing;
	double maxRefN;

 	//reference
	bool hasReference;
	bool chrChangedWindow;
	TFastaBuffer* fastaBuffer;

	//masks
	TBedReader* mask;

	//filters
	bool applyQualityFilter;
	uint8_t minQualityAsPhredInt, maxQualityAsPhredInt;
	int curReadGroupID;
	bool applyContextFilter;
	std::map<BaseContext,int> ignoreTheseContexts;
	bool applyDepthFilter;
	uint32_t readUpToDepth, minDepth, maxDepth;

	//blacklist
	bool _updateBlacklist;
	TAlignmentBlacklist blacklist;

	/*
	bool _writeBlackList;
	std::map <std::string, int> blacklist;
	gz::ogzstream ignoredReads;
	*/

	//limit chr and windows
	long limitWindows;
	int skipWindows;
	bool doLimitReads;
	long limitReads;

	//BAM file
	BAM::TBamFile* bamFile;

	//recalibration
	GenotypeLikelihoods::TGenotypeLikelihoodCalculator* genotypeLikelihoodCalculator;

	//contructor functions
	void _openBamFile(std::string filename, bool indexNotRequired);
	void _setWindowParameters(TParameters & params);
	void _setFilters(TParameters & params);
	void _setWindowFilters(TParameters & params);
	void _setAlignmentFilters(TParameters & params);
	void _setSiteFilters(TParameters & params);
	void _setBaseFilters(TParameters & params);
	void _setQualityFilters(int MinPhredInt, int MaxPhredInt);
	void _setMappingQualityFilters(int MinMQ, int MaxMQ);
	void _setContextFilter(std::vector<std::string> contexts);
	void _setFragmentLengthFilter(TParameters & params);
	void _setAlignmentFiltersToKeepAll();
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
	bool _fillAlignment(TAlignment & alignment);
	void _readAlignmentsIntoWindow(TWindow & window);
	void _applyWindowFilters(TWindow_base & window);
	void _recalibrate(TAlignment & alignment);
	void _adaptQualityWhenMerging(TBase & bestBase, TBase & worstBase, const bool & adaptQuality);

public:
	//alignment: goal is to make this private!


	bool doMasking, considerRegions;
	bool windowsPredefined;
	TBed* predefinedWindows;
	bool sitesProvided;
	TSiteSubset* subset = NULL;
 	TReadGroups readGroups;

	//maps
	TGenotypeMap genoMap;
	TQualityMap qualMap;

	//std::string outname;

	//construction
	TAlignmentParser();
	TAlignmentParser(TParameters & params, TBamFile * BamFile, GenotypeLikelihoods::TGenotypeLikelihoodCalculator* GenotypeLikelihoodCalculator, TLog* Logfile);
	~TAlignmentParser();
	void init(TParameters & params, TBamFile * BamFile, GenotypeLikelihoods::TGenotypeLikelihoodCalculator* GenotypeLikelihoodCalculator, TLog* Logfile);

	//getters
	int numReadGroups(){ return readGroups.size(); };
	unsigned int getWindowSize(){return windowSize;};
	int getMaxPhredInt(){return maxQualityAsPhredInt;};
	uint32_t getMaxDepth(){ return maxDepth; };

	//setters
	void setOutName(std::string outputName);
	void setParsingToTrue(){_parse = true;};
	void fillReferenceSequence(TFastaBuffer* fastaBuffer, TAlignment & alignment);

	//blacklist
	void setUpdateBlacklistToTrue();
	void setWriteBlacklistToFileToTrue();

	//functions to read and _parse
	void addReference(TFastaBuffer* FastaBuffer);


	//reading data requires windows
	bool readDataInNextWindow(TWindow & window);
	void downsampleWindow(TWindow_base & destination, TWindow & source, const double downsamplingProb, TRandomGenerator* randomGenerator);

	//reading data only requires alignments
	bool readNextAlignment(TAlignment & alignment); //to be used to go through bam file alignment by alignment

	//qualityTransformation
	//void initializeRecalibrationForQualityTransformation(TParameters & params);
	void mergeAlignedBasesBamReads(TAlignment* fwdAlignment, TAlignment* revAlignment, bool adaptQuality);
	void mergeAlignedBasesBamReadsRandom(TAlignment* fwdAlignment, TAlignment* revAlignment, bool adaptQuality, TRandomGenerator* randomGenerator);
	void mergeAlignedBasesOneRead(TAlignment* fwdAlignment, TAlignment* revAlignment, bool adaptQuality, TRandomGenerator* randomGenerator);
};


#endif /* TALIGNMENTPARSER_H_ */
