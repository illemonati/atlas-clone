/*
 * TAlignmentParser.h
 *
 *  Created on: Nov 14, 2017
 *      Author: wegmannd
 */

#ifndef TALIGNMENTPARSER_H_
#define TALIGNMENTPARSER_H_

#include "bamtools/api/BamReader.h"
#include "bamtools/api/BamWriter.h"
#include "bamtools/api/SamSequenceDictionary.h"
#include "bamtools/utils/bamtools_fasta.h"
#include "TGenotypeMap.h"
#include "TReadGroups.h"
#include "TLog.h"
#include "TPostMortemDamage.h"
#include "TAlignment.h"
#include "TWindow.h"
#include "TBed.h"
#include "TBedReader.h"
#include "TChromosomes.h"
#include "TFastaBuffer.h"
#include <vector>


//-----------------------------------------------------
//TAlignmentParser
//-----------------------------------------------------
class TAlignmentParser{
private:
	TLog* logfile;
	bool _keepDuplicates;
	bool _filterSoftClips;
	bool _keepImproperPairs;
	bool _keepUnmappedReads;
	bool _keepFailedQC;
	bool _keepSecondary;
	bool _keepSupplementary;
	bool _parse;

	unsigned int previousAlignmentPos;
	int previousAlignmentChr; //negative at beginning to trigger a chromosome change
	TAlignment* oldAlignment;
	bool oldAlignmentInitialized;
	bool oldAlignmentMustBeConsidered;

	//counters
	int totalNumberAlignmentsRead;
	int64_t sizeOfBamFile;

	//read trimming
	bool trimReads;
	int trimmingLength3Prime;
	int trimmingLength5Prime;

	//window params
	unsigned int windowSize;
	unsigned int numWindowsOnChr;
	unsigned int windowNumber;

	unsigned int maxReadLength;
	double maxMissing;
	double maxRefN;

	//masks
	TBedReader* mask;

	//filters
	bool applyQualityFilter;
	bool applyContextFilter;
	std::map<BaseContext,int> ignoreTheseContexts;
	size_t readUpToDepth, minDepth, maxDepth;
	bool applyMQFilter;
	int minMQ, maxMQ;
	bool applyFragmentLengthFilter;
	int minFragmentLength, maxFragmentLength;
	bool applyFragmentLengthLongerThanInsertSizeFilter;
	bool useStrand[2];
	bool useMate[2];

	//blacklist
	bool _updateBlacklist, _writeBlackList;
	std::map <std::string, int> blacklist;
	gz::ogzstream ignoredReads;

	//limit chr and windows
	long limitWindows;
	int skipWindows;
	int indexOfLimitChr;
	bool doLimitReads;
	long limitReads;
//	bool* useChromosome;

	//contructor functions
	void openBamFile(std::string filename);
	void setWindowParameters(TParameters & params);
	void setFilters(TParameters & params);
	void setMasks(TParameters & params);
//	void initializeSiteSubset(TParameters & params);
	void initializeReadGroups(TParameters & params);
	void setParsingLimits(TParameters & params);
	void setChrPloidy(TParameters & params);

	//move genome
	void jumpToEnd();
	void restartChromosomes(TWindow_base & window);
	void moveChromosome(TWindow_base & window);
	bool moveToNextWindowOnChr(TWindow_base & window);
	bool moveToNextPredefinedWindow(TWindow_base & window);
	bool moveWindow(TWindow_base & window);

	PMDType getEnumPMDType(std::string pmdType);
	void initializePostMortemDamage(TParameters & params);
	void initializeRecalibration(TParameters & params);

	bool readAlignment();
	bool applyFilters();
	bool fillAlignment(TAlignment & alignment);
	void readAlignmentsIntoWindow(TWindow & window);
	void applyWindowFilters(TWindow_base & window);
	void adaptQualityWhenMerging(TBase & bestBase, TBase & worstBase, const bool & adaptQuality);

public:
	//alignment: goal is to make this private!
	int curReadGroupID;
	int minQualForPrinting, maxQualForPrinting;
	int minQualityAsPhredInt, maxQualityAsPhredInt;
	bool doMasking, considerRegions;
	bool applyDepthFilter;
	bool windowsPredefined;
	TBed* predefinedWindows;
	bool sitesProvided;
	TSiteSubset* subset = NULL;
 	TReadGroups readGroups;

	BamTools::BamAlignment bamAlignment;

	//maps
	TGenotypeMap genoMap;
	TQualityMap qualMap;

	//BAM file
	std::string filename;
	std::string outname;
	BamTools::BamReader bamReader;
	BamTools::BamRegion bamRegion;
 	BamTools::SamHeader bamHeader;
 	TChromosomes chromosomes;
// 	BamTools::SamSequenceIterator chrIterator;

 	//reference
	bool hasReference;
	bool chrChangedAlignment;
	bool chrChangedWindow;
	BamTools::Fasta* fastaReference;
	TFastaBuffer fastaBuffer;

 	//recalibration
	TRecalibration* recalObject;
	bool doRecalibration;
	bool recalObjectInitialized;

	//PMD
	bool hasPMD;
	TPMD* pmdObjects;

	//construction
	TAlignmentParser();
	TAlignmentParser(int MaxReadLength, TParameters & params, TLog* Logfile);
	~TAlignmentParser();
	void init(int MaxReadLength, TParameters & params, TLog* Logfile);

	//getters
	bool qualitiesScoresAreRecalibrated(){ return recalObject->recalibrationChangesQualities(); };
	int numReadGroups(){ return readGroups.size(); };
	std::string recalibrationType(){ return recalObject->type(); };
	unsigned int getWindowSize(){return windowSize;};
	int getMaxPhredInt(){return maxQualityAsPhredInt;};
	int getNumAlignmentsRead(){ return totalNumberAlignmentsRead; };
	double getPositionInFile(){ return (double) bamReader.tell() / (double) sizeOfBamFile; };
	uint32_t getMaxDepth(){ return maxDepth; };

	//setters
	void setOutName(std::string outputName);
	void setQualityFilters(int MinPhredInt, int MaxPhredInt);
	void setMappingQualityFilters(int MinMQ, int MaxMQ);
	void setQualityRangeForPrinting(int minQual, int maxQual);
	void setContextFilter(std::vector<std::string> contexts);
	void setReadTrimming(int trim3Prime, int trim5Prime);
	void setApplyFragmentLengthLongerThanInsertSizeFilter(bool filterYesNo);
	void setFragmentLengthFilter(int MinFragmentLength, int MaxFragmentLength);

	void keepDuplicates(){_keepDuplicates = true;};
	void filterSoftClips(){_filterSoftClips = true;};
	void keepImproperPairs(){_keepImproperPairs = true;};
	void keepUnmappedReads(){_keepUnmappedReads = true;};
	void keepFailedQC(){_keepFailedQC = true;};
	void keepSecondaryReads(){_keepSecondary = true;};
	void keepSupplementaryReads(){_keepSupplementary = true;};
	void setParsingToTrue(){_parse = true;};
	void fillReferenceSequence(TFastaBuffer* fastaBuffer, TAlignment & alignment);
	std::string chrNumberToName(uint16_t chrNumber);
	uint32_t calcReferenceLength();
	std::string getCurChrName();
	uint16_t getCurRefId();
	uint32_t getCurChrLength();
	uint8_t getCurChrPloidy();

	//blacklist
	void setUpdateBlacklistToTrue(){
		_updateBlacklist = true;
		logfile->list("Storing ignored reads in a blacklist");
	};

	void setWriteBlacklistToFileToTrue(){
		_writeBlackList = true;
		std::string ignoredReadsFile = outname + "_ignoredReads.txt.gz";
		logfile->list("Writing ignored reads to '" + ignoredReadsFile + "'");
		ignoredReads.open(ignoredReadsFile.c_str());
		if(!ignoredReads) throw "Failed to open output file '" + ignoredReadsFile + "'!";
	}

	void addToBlacklist(TAlignment & alignment, const std::string & errorMessage){
		//TODO: should check if read already exists in blackfile (could be case in paired-end data) -> remove
		blacklist.emplace(alignment.alignmentName, 1);
		if(_writeBlackList){
			if(alignment.isReverseStrand){
				ignoredReads << "Read " << alignment.alignmentName << ", rev : " << errorMessage << "\n";
			} else {
				ignoredReads << "Read " << alignment.alignmentName << ", fwd : " << errorMessage << "\n";
			}
		}
	};
	void addToBlacklist(BamTools::BamAlignment & alignment, const std::string & errorMessage){
		//TODO: should check if read already exists in blackfile (could be case in paired-end data) -> remove
		blacklist.emplace(alignment.Name, 1);
		if(_writeBlackList){
			if(alignment.IsReverseStrand()){
				ignoredReads << "Read " << alignment.Name << ", rev : " << errorMessage << "\n";
			} else {
				ignoredReads << "Read " << alignment.Name << ", fwd : " << errorMessage << "\n";
			}
		}

	};

	void addToBlacklist(std::string & alignmentName, const std::string & errorMessage){
		//TODO: should check if read already exists in blackfile (could be case in paired-end data) -> remove
		blacklist.emplace(alignmentName, 1);
		if(_writeBlackList){
			ignoredReads << "Read " << alignmentName << " : " << errorMessage << "\n";
		}
	};

	void removeFromBlacklist(TAlignment & alignment, const std::string & errorMessage){
		blacklist.erase(alignment.alignmentName);
		if(_writeBlackList){
			if(alignment.isReverseStrand){
				ignoredReads << "Read " << alignment.alignmentName << ", rev : " << errorMessage << "\n";
			} else {
				ignoredReads << "Read " << alignment.alignmentName << ", fwd : " << errorMessage << "\n";
			}
		}
	};

	bool isInBlacklist(const std::string & alignmentName){
		if(blacklist.count(alignmentName) > 0)
			return true;
		return false;
	};

	//functions to read and _parse
	void checkAndFillAlingment(BamTools::BamAlignment& bamAlignment, TAlignment & alignment);
	void addReference(BamTools::Fasta* reference);
	void recalibrate(TAlignment & alignment);

	//reading data requires windows
	bool readDataInNextWindow(TWindow & window);
	void downsampleWindow(TWindow_base & destination, TWindow & source, const double downsamplingProb, TRandomGenerator* randomGenerator);

	//reading data only requires alignments
	bool readNextAlignment(TAlignment & alignment); //to be used to go through bam file alignment by alignment

	//qualityTransformation
	//void initializeRecalibrationForQualityTransformation(TParameters & params);
	void addSitesToQualityTransformTable(TAlignment & alignment, TQualityTransformTables & QTtables);
	void addSitesToQualityTransformTable(TAlignment & alignment, TRecalibration* otherRecalObject, TQualityTransformTables & QTtables);
	void mergeAlignedBasesBamReads(TAlignment* fwdAlignment, TAlignment* revAlignment, bool adaptQuality);
	void mergeAlignedBasesBamReadsRandom(TAlignment* fwdAlignment, TAlignment* revAlignment, bool adaptQuality, TRandomGenerator* randomGenerator);
	void mergeAlignedBasesOneRead(TAlignment* fwdAlignment, TAlignment* revAlignment, bool adaptQuality, TRandomGenerator* randomGenerator);
};

//-----------------------------------------------------
// TBamProgressReporter
//-----------------------------------------------------
class TBamProgressReporter{
private:
	timeval start, end;
	TAlignmentParser* parser;
	TLog* logfile;
	int progressFrequency;
	int lastProgressPrinted;

	void _init(int Frequency, TAlignmentParser* Parser, TLog* Logfile);
	std::string _getRunTime();
	void _printProgress();

public:
	TBamProgressReporter(int Frequency, TAlignmentParser* Parser, TLog* Logfile);
	TBamProgressReporter(TAlignmentParser* Parser, TLog* Logfile);

	void printProgress();
	void printEnd();
};


#endif /* TALIGNMENTPARSER_H_ */
