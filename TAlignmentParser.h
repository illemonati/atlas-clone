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
#include "TRecalibration.h"
#include "TPostMortemDamage.h"
#include "TAlignment.h"
#include "TWindow.h"
#include "TBed.h"
#include "TBedReader.h"
#include <vector>

//-----------------------------------------------------
//TFastaBuffer
//-----------------------------------------------------
//a buffer class to speed up adding the reference sequence to each read
//This class makes use of the fact that bam files are sorted, hence the buffer can always start at the current position

class TFastaBuffer{
private:
	BamTools::Fasta* reference;
	int bufferSize;
	std::string referenceSequence;

	int curChr;
	long curStart, curEnd;

	void moveTo(const int & chr, const int32_t & pos);

public:
	TFastaBuffer(BamTools::Fasta* Reference);
	~TFastaBuffer(){};
	void fill(const int & chr, const int32_t & start, const int32_t end, std::string & ref);
	int getCurChr(){return curChr;};
	long getCurStart(){return curStart;};
	long getCurEnd(){return curEnd;};
};

//-----------------------------------------------------
//TAlignmentParser
//-----------------------------------------------------
class TAlignmentParser{
private:

	//TReadGroups* readGroupTable;

	TLog* logfile;
	bool _keepDuplicates;
	bool _parse;
	int previousAlignmentPos;
	int previousAlignmentChr;
	TAlignment* oldAlignment;
	bool oldAlignmentInitialized;
	bool oldAlignmentMustBeConsidered;

	BamTools::BamAlignment bamAlignment;

	//read trimming
	bool trimReads;
	int trimmingLength3Prime;
	int trimmingLength5Prime;

	//iterators
 	int chrNumber;
 	long chrLength;

	//window params
	int windowSize;
	int numWindowsOnChr;
	int windowNumber;

	unsigned int maxReadLength;
	double maxMissing;
	double maxRefN;

	//masks
	TBedReader* mask;

	//filters
	bool applyQualityFilter;
	size_t readUpToDepth, minDepth, maxDepth;
	int minPhredInt, maxPhredInt;
	bool applyFragmentLengthFilter;

	//blacklist
	bool _updateBlacklist;
	std::map <std::string, int> blacklist;
	gz::ogzstream ignoredReads;

	//limit chr and windows
	long limitWindows;
	int limitChr;
	bool* useChromosome;

	//move genome
	void jumpToEnd();
	void restartChromosomes(TWindow & window);
//	bool iterateChromosome(TWindow & window);
	void moveChromosome(TWindow & window);
	bool moveToNextWindowOnChr(TWindow & window);
	bool moveToNextPredefinedWindow(TWindow & window);
	bool moveWindow(TWindow & window);

	PMDType getEnumPMDType(std::string pmdType);
	void initializePostMortemDamage(TParameters & params);
	void initializeRecalibration(TParameters & params);

	bool readAlignment();
	bool applyFilters();
	void fillAlignment(TAlignment & alignment);
	void readAlignmentsIntoWindow(TWindow & window);
	void applyFilters(TWindow & window);

public:
	//alignment: goal is to make this private!
	int curReadGroupID;
	int minQualForPrinting, maxQualForPrinting;
	int minQual, maxQual;
	bool doMasking, considerRegions;
	bool doCpGMasking;
	bool applyDepthFilter;
	bool windowsPredefined;
	TBed* predefinedWindows;
	bool sitesProvided;
	TSiteSubset* subset = NULL;
 	TReadGroups readGroups;

	//maps
	TGenotypeMap genoMap;
	TQualityMap qualMap;

	//BAM file
	std::string filename;
	BamTools::BamReader bamReader;
	BamTools::BamRegion bamRegion;
 	BamTools::SamHeader bamHeader;
 	BamTools::SamSequenceIterator chrIterator;

 	//reference
	bool hasReference;
	bool chrChanged;
	BamTools::Fasta* fastaReference;
	TFastaBuffer* fastaBuffer;

 	//recalibration
	TRecalibration* recalObject;
	TRecalibration* recalObject2;
	bool doRecalibration;
	bool doRecalibration2;
	bool recalObjectInitialized;
	bool recalObjectInitialized2;

	//PMD
	bool hasPMD;
	TPMD* pmdObjects;

	//construction
	TAlignmentParser();
	TAlignmentParser(int MaxReadLength, TParameters & params, TLog* Logfile);
	~TAlignmentParser();
	void init(int MaxReadLength, TParameters & params, TLog* Logfile);

	//getters
	int getWindowSize(){return windowSize;}
	int getMaxPhredInt(){return maxPhredInt;}

	//setters
	void setQualityFilters(int minQual, int maxQual);
	void setQualityRangeForPrinting(int minQual, int maxQual);
	void setReadTrimming(int trim3Prime, int trim5Prime);
	void setApplyFragmentLengthFilter(bool filterYesNo);

	void keepDuplicates(){_keepDuplicates = true;};
	void setParsingToTrue(){_parse = true;};
	void fillReferenceSequence(TFastaBuffer* fastaBuffer, TAlignment & alignment);
	std::string chrNumberToName(int chrNumber);
	long calcReferenceLength();

	//blacklist
	void setUpdateBlacklistToTrue(){
		_updateBlacklist = true;
		std::string ignoredReadsFile = extractBeforeLast(filename, ".bam") + "_ignoredReads.txt.gz";
		logfile->list("Writing sequencing depth estimates to '" + ignoredReadsFile + "'");
		ignoredReads.open(ignoredReadsFile.c_str());
		if(!ignoredReads) throw "Failed to open output file '" + ignoredReadsFile + "'!";
	};
	void addToBlacklist(TAlignment & alignment, const std::string & errorMessage){
		//TODO: should check if read already exists in blackfile (could be case in paired-end data) -> remove
		blacklist.emplace(alignment.alignmentName, 1);
		ignoredReads << "Rea	d " << alignment.alignmentName << " isReverse=" << alignment.isReverseStrand << " : " << errorMessage << "\n";
	}
	void addToBlacklist(BamTools::BamAlignment & alignment, const std::string & errorMessage){
		//TODO: should check if read already exists in blackfile (could be case in paired-end data) -> remove
		blacklist.emplace(alignment.Name, 1);
		ignoredReads << "Read " << alignment.Name << " isReverse=" << alignment.IsReverseStrand() << " : " << errorMessage << "\n";
	}
	void addToBlacklist(std::string & alignmentName, const std::string & errorMessage){
		//TODO: should check if read already exists in blackfile (could be case in paired-end data) -> remove
		blacklist.emplace(alignmentName, 1);
		ignoredReads << "Read " << alignmentName << " : " << errorMessage << "\n";
	}
	void removeFromBlacklist(TAlignment & alignment, const std::string & errorMessage){
		blacklist.erase(alignment.alignmentName);
		ignoredReads << "Read " << alignment.alignmentName << " isReverse=" << alignment.isReverseStrand << " : " << errorMessage << "\n";
	}
	bool isInBlacklist(std::string & alignmentName){
		if(blacklist.count(alignmentName) > 0)
			return true;
		return false;
	}

	//functions to read and _parse
	void checkAndFillAlingment(BamTools::BamAlignment& bamAlignment, TAlignment & alignment);
	void addReference(BamTools::Fasta* reference);
	void recalibrate(TAlignment & alignment);

	//reading data requires windows
	bool readDataInNextWindow(TWindow & window);

	//reading data only requires alignments
	bool readNextAlignment(TAlignment & alignment); //to be used to go through bam file alignment by alignment
	bool readNextAlignmentWithBlacklist(TAlignment & alignment);

	//qualityTransformation
	void initializeRecalibrationForQualityTransformation(TParameters & params);
	void addSitesToQualityTransformTable(TAlignment & alignment, TRecalibration* recalObject, std::vector<TQualityTransformTable*> & QTtables, TLog* logfile);
	void addSitesToQualityTransformTable(TAlignment & alignment, TRecalibration* recalObject, TRecalibration* otherRecalObject, std::vector<TQualityTransformTable*> & QTtables, TLog* logfile);
	void mergeAlignedBasesBamReads(TAlignment* fwdAlignment, TAlignment* revAlignment, bool adaptQuality);



//	bool addReadToWindow(TAlignmentParser & alignemntParser, TPMD* pmdObjects);

};

#endif /* TALIGNMENTPARSER_H_ */
