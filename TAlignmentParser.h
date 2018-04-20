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
	void fill(const int & chr, const int32_t & start, const int32_t end, std::string & ref);
};


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

	//quality filter
	bool applyQualityFilter;
	int minQual, maxQual;
	int minQualForPrinting, maxQualForPrinting;

	//read trimming
	bool trimReads;
	int trimmingLength3Prime;
	int trimmingLength5Prime;

	//tmp variables
	std::vector<BamTools::CigarOp>::const_iterator cigarIter;
	std::vector<BamTools::CigarOp>::const_iterator cigarEnd;
	std::string tmpString;
	std::string tmpString2;
	std::string referenceSequence;
	std::string::iterator stringIt;

	//reference
	bool hasReference;
	TFastaBuffer* fastaBuffer;

public:
	//alignment: goal is to make this private!
	BamTools::BamAlignment bamAlignment;

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
	void fillReferenceSequence(TFastaBuffer* fastaBuffer, TAlignment & alignment);
	void setQualityFilters(int minQual, int maxQual);
	void setQualityRangeForPrinting(int minQual, int maxQual);
	void setReadTrimming(int trim3Prime, int trim5Prime);

	//functions to read and parse
	bool readAlignment(BamTools::BamReader & bamReader, TAlignment & alignment, BamTools::Fasta* reference, bool & parse, bool & addReference, bool & filterBaseQual, bool & trim);
	void addReference(BamTools::Fasta* reference);
};



#endif /* TALIGNMENTPARSER_H_ */
