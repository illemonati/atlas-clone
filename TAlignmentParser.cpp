/*
 * TAlignmentParser.cpp
 *
 *  Created on: Nov 14, 2017
 *      Author: wegmannd
 */

#include "TAlignmentParser.h"
#include "TWindow.h"

//-----------------------------------------------------
//TFastaBuffer
//-----------------------------------------------------
TFastaBuffer::TFastaBuffer(BamTools::Fasta* Reference){
	bufferSize = 10000000;
	reference = Reference;
	referenceSequence = "";
	curStart = -1;
	curChr = -1;
	curEnd = -1;
}


void TFastaBuffer::moveTo(const int & chr, const int32_t & pos){
	curChr = chr;
	curStart = pos;
	curEnd = pos + bufferSize;
	if(!reference->GetSequence(chr, curStart, curEnd, referenceSequence))
		throw "Problem reading " + toString(chr) + ":" + toString(curStart) + "-" + toString(curEnd) + " from fasta file!";
}

void TFastaBuffer::fill(const int & chr, const int32_t & start, const int32_t end, std::string & ref){
	//move buffer, if necessary
	if(chr != curChr || end > curEnd || start < curStart)
		moveTo(chr, start);

	//now copy to string
	ref.assign(referenceSequence, start - curStart, end - start + 1);
}

//-----------------------------------------------------
//TAlignmentParser
//-----------------------------------------------------
TAlignmentParser::TAlignmentParser(){
	maxSize = 0;
	readGroupTable = NULL;
	logfile = NULL;
	_keepDuplicates = false;
	parse = false;

	//filters
	applyQualityFilter = false;
	minQual = 33;
	maxQual = 126;
	minQualForPrinting = 33;
	maxQualForPrinting = 126;
	trimReads = false;
	trimmingLength3Prime = 0;
	trimmingLength5Prime = 0;

	//reference
	hasReference = false;
	fastaBuffer = NULL;
};

TAlignmentParser::TAlignmentParser(TReadGroups* ReadGroupTable, unsigned int MaxSize, TLog* Logfile){
	TAlignmentParser();
	init(ReadGroupTable, MaxSize, Logfile);
};

void TAlignmentParser::init(TReadGroups* ReadGroupTable, unsigned int MaxSize, TLog* Logfile){
	maxSize = MaxSize;
	readGroupTable = ReadGroupTable;
	logfile = Logfile;
};

void TAlignmentParser::setQualityFilters(int MinQual, int MaxQual){
	applyQualityFilter = true;
	minQual = MinQual;
	maxQual = MaxQual;
};

void TAlignmentParser::setQualityRangeForPrinting(int minQual, int maxQual){
	minQualForPrinting = minQual;
	maxQualForPrinting = maxQual;
};

void TAlignmentParser::setReadTrimming(int trim3Prime, int trim5Prime){
	trimmingLength3Prime = trim3Prime;
	trimmingLength5Prime = trim5Prime;
	trimReads = true;
};

void TAlignmentParser::addReference(BamTools::Fasta* reference){
	hasReference = true;
	fastaBuffer = new TFastaBuffer(reference);
};

void TAlignmentParser::fillReferenceSequence(TFastaBuffer* fastaBuffer, TAlignment & alignment){
	if(!hasReference) //is this check really necessary?
		throw "No reference provided!";

	fastaBuffer->fill(alignment.chrNumber, alignment.position, alignment.position + alignment.alignedPos[alignment.length-1], referenceSequence);
};


//------------------------------
//public functions
//------------------------------
bool TAlignmentParser::readAlignment(BamTools::BamReader & bamReader, TAlignment* alignment){
	//make sure container is empty
	alignment->clear();

	if(!bamReader.GetNextAlignment(bamAlignment))
		return false;

	//check read length
	if(bamAlignment.AlignedBases.size() > maxSize)
		throw "Alignment '" +  bamAlignment.Name + "' is longer than the max read length! Please change max read length to parse this data.";

	//fill alignment
	std::string readGroup;
	bamAlignment.GetTag("RG", readGroup);
	int readGroupId = readGroupTable->find(readGroup);
	alignment->fill(bamAlignment, readGroupId);

	//check if insert size is shorter than read, this means we are reading the adaptor sequence
	bool filtersPassed = true;
	if(bamAlignment.IsPaired() && abs(bamAlignment.InsertSize) <= bamAlignment.AlignedBases.length()){
		logfile->warning("The following alignment is longer than its insert size: " + bamAlignment.Name);
		filtersPassed = false;
	} else {
		//apply filters: read group in use and basic QC
		filtersPassed = readGroupTable->readGroupInUse(readGroupId)
						&& bamAlignment.IsMapped() && !bamAlignment.IsFailedQC()
						&& bamAlignment.IsPrimaryAlignment()
						&& (_keepDuplicates || !bamAlignment.IsDuplicate());
	}

	//set filter
	alignment.setFiltersPassed(filtersPassed);

	if(parse){
		alignment.parse(genoMap, qualityMap);
		if(hasReference)
			fillReferenceSequence(fastaBuffer, alignment);
		if(applyQualityFilter)
			alignment.filterForBaseQuality(minQual, maxQual);
	}

	return true;
}


//---------------------
//read data in windows
//---------------------


bool TAlignmentParser::readData(BamTools::BamReader & bamReader, TWindow & window){
	logfile->listFlush("Reading data ...");

	TAlignment* alignmentP;
	setParsingToTrue();

	//measure runtime
	struct timeval start, end;
	gettimeofday(&start, NULL);

	//decide which alignment to fill (from empty stack or create a new one)
	if(emptyAlignments.size() > 0)
		alignmentP = emptyAlignments.begin();
	else
		alignmentP = new TAlignment;

	//pass it to readAlignment
	readAlignment(bamReader, alignmentP);

	//if alignment.passedFilters & readGroupInUse:
		//addAlignementToWindow
	if(alignmentP->passedFilters && readGroups.readGroupInUse(alignmentP->readGroupId)){
		//and add to windows
		if(!addAlignementToWindows(alignment, window)){
			//read is beyond window and should be reconsidered
			oldAlignementMustBeConsidered = true;
			break;
		}
	}
	//TODO: stacks should be in TWindow, but parser checks if alignment is still in window

	//In addAlignementToWindow:
		//check if it still overlaps with current window
		//if no, return false, break readData loop
		//if yes, add all relevant positions$
		//if only beginning, add relevant positions and break readData loop; clear window
		//if only end,

	//In clear window:
		//move all alignments that have start+alignedLength < window end to emptyAlignments

	//apply all masks at the end of break readData loop




	gettimeofday(&end, NULL);
	logfile->write(" done (in " , end.tv_sec  - start.tv_sec, "s)!");

/*	if(windowPair.cur->numReadsInWindow > 0){
		//apply masks and filters
		if(doMasking){
			logfile->listFlush("Masking sites ...");
			windowPair.cur->applyMask(mask, considerRegions);
			logfile->done();
		} else if(considerRegions){
			logfile->listFlush("Masking sites outside regions ...");
			windowPair.cur->applyMask(mask, considerRegions);
			logfile->done();
		} else if(doCpGMasking){
			logfile->listFlush("Masking CpG sites ...");
			windowPair.cur->maskCpG(reference, chrNumber);
			logfile->done();
		} if(applyDepthFilter){
			windowPair.cur->applyDepthFilter(minDepth, maxDepth);
		} if(maxRefN < 1.0 && fastaReference == true){
			windowPair.cur->addReferenceBaseToSites(reference, chrNumber);
			windowPair.cur->calcFracN();
		}

		//calc sequencing depth
		windowPair.cur->calcDepth();

		//report
		logfile->conclude("read data from " + toString(windowPair.cur->numReadsInWindow) + " reads.");
		logfile->conclude("sequencing depth is " + toString(windowPair.cur->depth));
		logfile->conclude(toString(windowPair.cur->fractionsitesDepthAtLeastTwo * 100) + "% of all sites are covered at least twice");
		logfile->conclude(toString(windowPair.cur->fractionSitesNoData * 100) + "% of all sites have no data");
		if(windowPair.cur->fractionSitesNoData > maxMissing){
			logfile->conclude("Level of missing data > threshold of " + toString(maxMissing) + " -> skipping this window");
			return false;
		}
		if(maxRefN < 1.0 && fastaReference == true){
			logfile->conclude(toString(windowPair.cur->fractionRefIsN * 100) + "% of all reference bases are 'N'");
			if(windowPair.cur->fractionRefIsN > maxRefN){
				logfile->conclude("Fraction of 'N' in reference > threshold of " + toString(maxRefN) + " -> skipping this window");
				return false;
			}
		}
		return true;
	} else {
		logfile->conclude("No data in this window.");
		return false;
	}*/
};


bool TAlignmentParser::addAlignementToWindows(TAlignment & alignment, TWindow & window, int oldPos, int curStart, int curEnd){
	//check if bam file is sorted
	if(alignment.getPosition() < oldPos)
		throw "BAM file must be sorted by position!";
	oldPos = alignment.getPosition();

	//and add
	if(alignment.position >= curStart){
		//check if still within current window and add to window
		if(alignment.position >= curEnd){
			return false;
		}
		else {
			//check if there are emptyAlignments
			if(emptyAlignments.size() > 0) emptyAlignments.pop_back()
			else {

				window.addFromRead())
			}

		}
	}
	return true; //continue
}




/*bool TAlignmentParser::addFromRead(TAlignment & alignment){
	 Note:
	 * Function returns true if read also maps to next window and
	 * returns false if end of read is within this (or a previous) window


	//TODO: move to alignment parser instead?

	//check if alignment is inside window
	if(alignment.position >= end) return true;
	if(alignment.position + alignment.length < start) return false;

	//find which position to consider first
	++numReadsInWindow;
	int firstPos = alignment.position - start;

	//std::cout << "[" << start << "," << end <<  "]: firstPos = " << firstPos;

	int p = 0;

	if(firstPos < 0){
		while(p < alignment.length && (firstPos + alignment.alignedPos[p]) < 0)
			++p;
		if(p == alignment.length)
			return false;
	}
	int internalPos;

	 Note:
	 *  1) Reference is 5' -> 3'
	 *  2) distance is 0-based!
	 *  3) Ignoring indels in other mate when calculating distances
	 *  4) Function add needs first P(C->T), then P(G->A)


	for(; p < alignment.length; ++p){
		if(alignment.aligned[p] && alignment.bases[p].base != N){
			internalPos = firstPos + alignment.alignedPos[p];
			if(internalPos >= length)
				return true; //since part of the read maps to next window
			sites[internalPos].add(&alignment.bases[p]);
		}
	}

	return false;
}*/

