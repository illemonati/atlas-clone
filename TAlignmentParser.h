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

//-----------------------------------------------------
//TAlignmentParser
//-----------------------------------------------------
class TAlignmentParser{
private:
	//variables
	unsigned int maxSize;
	TGenotypeMap genoMap;
	TQualityMap qualityMap;
	TReadGroups* readGroupTable;
	TLog* logfile;
	bool _keepDuplicates;
	bool parse;

	//quality filter
	bool applyQualityFilter;

	//read trimming
	bool trimReads;
	int trimmingLength3Prime;
	int trimmingLength5Prime;

	//tmp variables
	std::string referenceSequence;

	//reference
	bool hasReference;

public:
	//alignment: goal is to make this private!
	BamTools::BamAlignment bamAlignment;
	int minQualForPrinting, maxQualForPrinting;
	int minQual, maxQual;
	TFastaBuffer* fastaBuffer;

	//construction
	TAlignmentParser();
	TAlignmentParser(TReadGroups* readGroupTable, unsigned int MaxSize, TLog* Logfile);
	~TAlignmentParser(){
		if(hasReference)
			delete fastaBuffer;
	}
	void init(TReadGroups* readGroupTable, unsigned int MaxSize, TLog* Logfile);

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
};

#endif /* TALIGNMENTPARSER_H_ */
