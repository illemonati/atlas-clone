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

//-----------------------------------------------------
//TAlignmentParser
//-----------------------------------------------------
class TAlignmentParser{
private:
	//variables
	TGenotypeMap genoMap;
	TQualityMap qualityMap;
	//TReadGroups* readGroupTable;
 	TReadGroups readGroups;

	TLog* logfile;
	bool _keepDuplicates;
	bool parse;
	int previousAlignmentPos;
	int previousAlignmentChr;
	TAlignment* oldAlignment;
	bool oldAlignmentInitialized;
	bool oldAlignmentMustBeConsidered;

	//quality filter
	bool applyQualityFilter;

	//read trimming
	bool trimReads;
	int trimmingLength3Prime;
	int trimmingLength5Prime;

	//reference
	BamTools::Fasta* fastaReference;
	std::string referenceSequence;
	bool hasReference;

	//iterators
 	int chrNumber;
 	long chrLength;
 //	long curStart;
 //	long curEnd;

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
	size_t minDepth, maxDepth;
	int minPhredInt, maxPhredInt;

	//limit chr and windows
	long limitWindows;
	int limitChr;
	bool* useChromosome;

	//recal and pmd objects
	TPMD* pmdObjects;
	bool hasPMD;
	TRecalibration* recalObject2;
	bool doRecalibration;
	bool doRecalibration2;
	bool recalObjectInitialized;
	bool recalObjectInitialized2;

	//move genome
	void jumpToEnd();
	void restartChromosome(TWindow & window);
//	bool iterateChromosome(TWindow & window);
	void moveChromosome(TWindow & window);
	bool moveToNextWindowOnChr(TWindow & window);
	bool moveWindow(TWindow & window);

	void initializePostMortemDamage(TParameters & params);
	void initializeRecalibration(TParameters & params);

	bool readAlignment();
	void fillAlignment(TAlignment & alignment);
	void readAlignmentsIntoWindow(TWindow & window);
	void applyFilters(TWindow & window);

public:
	//alignment: goal is to make this private!
	BamTools::BamAlignment bamAlignment;
	int curReadGroupID;
	int minQualForPrinting, maxQualForPrinting;
	int minQual, maxQual;
	TFastaBuffer* fastaBuffer;
	bool doMasking, considerRegions;
	bool doCpGMasking;
	bool applyDepthFilter;
	bool windowsPredefined;
	TBed* predefinedWindows;

	TRecalibration* recalObject;


	//BAM file
	std::string filename;
	BamTools::BamReader bamReader;
	BamTools::BamRegion bamRegion;
 	BamTools::SamHeader bamHeader;
 	BamTools::SamSequenceIterator chrIterator;

	//construction
	TAlignmentParser();
	TAlignmentParser(TParameters & params, TLog* Logfile);
	~TAlignmentParser();
	void init(TParameters & params, TLog* Logfile);

	//setters
	void keepDuplicates(){_keepDuplicates = true;};
	void setParsingToTrue(){parse = true;};
	void fillReferenceSequence(TFastaBuffer* fastaBuffer, TAlignment & alignment);
	std::string chrNumberToName(int chrNumber);
	void setQualityFilters(int minQual, int maxQual);
	void setQualityRangeForPrinting(int minQual, int maxQual);
	void setReadTrimming(int trim3Prime, int trim5Prime);

	//functions to read and parse
	bool readNextAligment(TAlignment & alignment); //to be used to go through bam file alignment by alignment
	void checkAndFillAlingment(BamTools::BamAlignment& bamAlignment, TAlignment & alignment);
	void addReference(BamTools::Fasta* reference);

	//read data in windows
	bool readDataInNextWindow(TWindow & window);


//	bool addReadToWindow(TAlignmentParser & alignemntParser, TPMD* pmdObjects);

};

#endif /* TALIGNMENTPARSER_H_ */
