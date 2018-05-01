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

	void initializePostMortemDamage(TParameters & params);
	void initializeRecalibration(TParameters & params);

public:
	//alignment: goal is to make this private!
	BamTools::BamAlignment bamAlignment;
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
	~TAlignmentParser(){
		if(hasReference)
			delete fastaBuffer;
		if(doMasking)
			delete mask;
		if(windowsPredefined)
			delete predefinedWindows;
		if(useChromosome)
			delete[] useChromosome;
		if(recalObjectInitialized) delete recalObject;
		if(pmdObjects) delete[] pmdObjects;
	}
	void init(TParameters & params, TLog* Logfile);

	//setters
	void keepDuplicates(){_keepDuplicates = true;};
	void setParsingToTrue(){parse = true;};
	void fillReferenceSequence(TFastaBuffer* fastaBuffer, TAlignment & alignment);
	void setQualityFilters(int minQual, int maxQual);
	void setQualityRangeForPrinting(int minQual, int maxQual);
	void setReadTrimming(int trim3Prime, int trim5Prime);

	//functions to read and parse
	bool readAlignment(BamTools::BamReader & bamReader, TAlignment & alignment);
	void checkAndFillAlingment(BamTools::BamAlignment& bamAlignment, TAlignment & alignment);
	void addReference(BamTools::Fasta* reference);

	//read data in windows
	bool nextWindow(TWindow & window);
	bool readDataInWindow(TWindow & window);
	bool checkPosition(TAlignment & alignment, TWindow & window);
	bool readAlignmentsIntoWindow(TWindow & window, TReadGroups & readGroups);
	bool applyFilters(TWindow & window);
//	bool addReadToWindow(TAlignmentParser & alignemntParser, TPMD* pmdObjects);

};

#endif /* TALIGNMENTPARSER_H_ */
