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

//-----------------------------------------------------
//TAlignmentParser
//-----------------------------------------------------
class TAlignmentParser{
private:
	//variables
	TGenotypeMap genoMap;
	TQualityMap qualityMap;
	TReadGroups* readGroupTable;
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

	//tmp variables

	//reference
	BamTools::Fasta* fastaReference;
	std::string referenceSequence;
	bool hasReference;

	//window params
	bool windowsPredefined;
	TBed* predefinedWindows;
	int windowSize;
	int numWindowsOnChr;
	int windowNumber;
	unsigned int maxReadLength;
	double maxMissing;
	double maxRefN;

	//masks
	TBedReader* mask;
	bool doMasking, considerRegions;
	bool doCpGMasking;

	//filters
	bool applyDepthFilter;
	size_t minDepth, maxDepth;
	int minPhredInt, maxPhredInt;
	int minOutQual, maxOutQual;



public:
	//alignment: goal is to make this private!
	BamTools::BamAlignment bamAlignment;
	int minQualForPrinting, maxQualForPrinting;
	int minQual, maxQual;
	TFastaBuffer* fastaBuffer;

	//construction
	TAlignmentParser();
	TAlignmentParser(TReadGroups* readGroupTable, TParameters & params, TLog* Logfile);
	~TAlignmentParser(){
		if(hasReference)
			delete fastaBuffer;
		if(doMasking)
			delete mask;

	}
	void init(TReadGroups* readGroupTable, TParameters & params, TLog* Logfile);

	//setters
	void keepDuplicates(){_keepDuplicates = true;};
	void setParsingToTrue(){parse = true;};
	void fillReferenceSequence(TFastaBuffer* fastaBuffer, TAlignment & alignment);
	void setQualityFilters(int minQual, int maxQual);
	void setQualityRangeForPrinting(int minQual, int maxQual);
	void setReadTrimming(int trim3Prime, int trim5Prime);

	//functions to read and parse
	bool readAlignment(BamTools::BamReader & bamReader, TAlignment & alignment);
	void addReference(BamTools::Fasta* reference);

	//read data in windows
	bool addAlignementToWindow(TAlignment & alignment, TWindow & window);
	bool readData(BamTools::BamReader & bamReader, TWindow & window, TReadGroups & readGroups);
//	bool addReadToWindow(TAlignmentParser & alignemntParser, TPMD* pmdObjects);

};

#endif /* TALIGNMENTPARSER_H_ */
