/*
 * TBamFile.cpp
 *
 *  Created on: May 18, 2020
 *      Author: wegmannd
 */

#include "TBamFile.h"

//-----------------------------------------------------
//TBamFile
//-----------------------------------------------------
TBamFile::TBamFile(){
	_maxReadLength = 200;
	_fileSize = 0;
	_numAlignmentRead = 0;
	_numAlignmentsPassedQC = 0;

	//current alignment
	_curReadGroupID = 0;
	_previousAlignmentPos = 0;
	_previousAlignmentChr = -1;
	_chrChanged = false;

	//writing
	_openForWriting = false;

	//set filters to default
	_QCFiltersPassed = false;
	_keepAll = false;

	//blacklist
	_updateBlacklist = false;
	_blacklist = nullptr;
};

void TBamFile::setAlignmentFilters(TParameters & params, TLog* logfile){
	if(params.parameterExists("keepAllReads")){
		_keepAll = true;
		logfile->list("Will keep all reads. (parameter 'keepDuplicates', overrules any other QC filter)");
	} else {
		//duplicates
		if(params.parameterExists("keepDuplicates")){
			_duplicateFilter.keep();
			logfile->list("Will keep duplicate reads. (parameter 'keepDuplicates')");
		} else {
			_duplicateFilter.filter("duplicate");
			logfile->list("Will filter out duplicate reads. (use 'keepDuplicates' to keep)");
		}

		//soft clips
		if(params.parameterExists("filterSoftClips")){
			_softClippedFilter.filter("soft clipped");
			logfile->list("Will filter out soft clipped reads. (parameter 'filterSoftClips')");
		} else {
			_softClippedFilter.keep();
			logfile->list("Will keep soft clipped reads. (use 'filterSoftClips' to filter out)");
		}

		//improper pairs
		if(params.parameterExists("keepImproperPairs")){
			_improperPairsFilter.keep();
			logfile->list("Will keep improper pairs. (parameter 'keepImproperPairs')");
		} else {
			_improperPairsFilter.filter("improper pair");
			logfile->list("Will filter out improper pairs. (use 'keepImproperPairs' to keep)");
		}

		//unmapped reads
		if(params.parameterExists("keepUnmappedReads")){
			_unmappedFilter.keep();
			logfile->list("Will keep unmapped reads. (parameter 'keepUnmappedReads')");
		} else {
			_unmappedFilter.filter("unmapped");
			logfile->list("Will filter out unmapped reads. (use 'keepUnmappedReads' to keep)");
		}

		//failed QC
		if(params.parameterExists("keepFailedQC")){
			_failedQCFilter.keep();
			logfile->list("Will keep reads that failed QC. (parameter 'keepFailedQC')");
		} else {
			_failedQCFilter.filter("failed QC");
			logfile->list("Will filter out reads that failed QC. (use 'keepFailedQC' to keep)");
		}

		//secondary reads
		if(params.parameterExists("keepSecondaryReads")){
			_secondaryFilter.keep();
			logfile->list("Will keep secondary reads. (parameter 'keepSecondaryReads')");
		} else {
			_secondaryFilter.filter("secondary alignment");
			logfile->list("Will filter out secondary reads. (use 'keepSecondaryReads' to keep)");
		}

		//supplementary reads
		if(params.parameterExists("keepSupplementaryReads")){
			_supplementaryFilter.keep();
			logfile->list("Will keep supplementary reads. (parameter 'keepSupplementaryReads')");
		} else {
			_supplementaryFilter.filter("supplementary alignment");
			logfile->list("Will filter out supplementary reads. (use 'keepSupplementaryReads' to keep)");
		}

		//fragment length
		if(params.parameterExists("keepReadsLongerThanFragment")){
			_longerThanFragmentFilter.keep();
			logfile->list("Will keep reads that are longer than the fragment size. (parameter 'keepReadsLongerThanFragment')");
		} else {
			_longerThanFragmentFilter.filter("longer than fragment");
			logfile->list("Will filter out reads that are longer than the fragment size. (use 'keepReadsLongerThanFragment' to keep)");
		}

		//strand
		if(params.parameterExists("keepOnlyFwd")){
			_fwdStrandFilter.keep();
			_revStrandFilter.filter("reverse strand");
			logfile->list("Will keep only forward mapping reads. (parameter 'keepOnlyFwd')");
		}
		else if(params.parameterExists("keepOnlyRev")){
			_fwdStrandFilter.filter("forward strand");
			_revStrandFilter.keep();
			logfile->list("Will keep only reverse mapping reads. (parameter 'keepOnlyRev')");
		} else {
			_fwdStrandFilter.keep();
			_revStrandFilter.keep();
			logfile->list("Will keep forward and reverse mapping reads. (use 'keepOnlyFwd' or 'keepOnlyRev' to limit)");
		}

		//mate
		if(params.parameterExists("keepOnlyFirst")){
			_firstMateFilter.keep();
			_secondMateFilter.filter("second mate");
			logfile->list("Will keep only the first mates. (parameter 'keepOnlyFirst')");
		}
		else if(params.parameterExists("keepOnlySecond")){
			_firstMateFilter.filter("first mate");
			_secondMateFilter.keep();
			logfile->list("Will keep only the second mates. (parameter 'keepOnlySecond')");
		} else {
			_firstMateFilter.keep();
			_secondMateFilter.keep();
			logfile->list("Will keep first and second mates. (use 'keepOnlyFirst' or 'keepOnlySecond' to limit)");
		}

		//Mapping quality filter
		if(params.parameterExists("minMQ") || params.parameterExists("maxMQ")){
			int MinMQ = params.getParameterInt("minMQ", 0);
			int MaxMQ = params.getParameterInt("maxMQ", 254);

			if(MinMQ < 0 || MinMQ > 254)
				throw "minMQ '" + toString(MinMQ) + "' is outside the accepted range [0,254]!";
			if(MaxMQ < 0 || MaxMQ > 254)
				throw "maxMQ '" + toString(MaxMQ) + "' is outside the accepted range [0,254]!";
			if(MaxMQ < MinMQ)
				throw "maxMQ must be <= minMQ";

			_mappingQualityFilter.filter(MinMQ, MaxMQ, "mapping quality outside [" + toString(MinMQ) + ", " + toString(MaxMQ) + "]");
			logfile->list("Will filter out reads with mapping quality outside the range [" + toString(MinMQ) + ", " + toString(MaxMQ) + "]. (parameters 'minMQ', 'maxMQ')");
		} else {
			_mappingQualityFilter.keep();
			logfile->list("Will keep reads regardless of their mapping quality. (use 'minMQ' and 'maxMQ' to limit)");
		}

		//Fragment length filter
		if(params.parameterExists("minFragmentLength") || params.parameterExists("maxFragmentLength")){
			int MinFragmentLength = params.getParameterIntWithDefault("minFragmentLength", 0);
			int MaxFragmentLength = params.getParameterIntWithDefault("maxFragmentLength", 1000);
			if(MinFragmentLength < 0)
				throw "minFragmentLength '" + toString(MinFragmentLength) + "' must be >0 0!";
			if(MinFragmentLength < 0)
				throw "minFragmentLength '" + toString(MinFragmentLength) + "' must be >0 0!";
			if(MinFragmentLength > MaxFragmentLength)
				throw "minFragmentLength must be <= maxFragmentLength!";



			logfile->list("Will filter out reads with fragment length outside the range [" + toString(MinFragmentLength) + ", " + toString(MaxFragmentLength) + "]. (parameters 'minFragmentLength', 'maxFragmentLength')");
		} else {
			_fragmentLengthfilter.keep();
			logfile->list("Will keep reads reagrless of their fragment length. (use 'minFragmentLength', 'maxFragmentLength' to limit)");
		}

		//read grup
		_readGroupFilter.filter("read group not in use");

	}
};

void TBamFile::limitReadGroups(std::string readGroupList, TLog* logfile){
	readGroups.filterReadGroups(readGroupList);
	logfile->startIndent("Will limit analysis to the following read groups:");
	readGroups.printReadgroupsInUse(logfile);
	_readGroupFilter.filter("read group not in use");
	logfile->endIndent();
};

//--------------------------------------------------------
// Functions for reading
//--------------------------------------------------------

void TBamFile::open(const std::string filename, const uint32_t maxReadLength, const bool indexNotRequired){
	//open BAM file
	if (!_bamReader.Open(filename))
		throw "Failed to open BAM file '" + filename + "'!";
	//load index file
	if(!_bamReader.LocateIndex() && !indexNotRequired)
		throw "No index file found for BAM file '" + filename + "'!";

	//initialize bam stuff
	_bamHeader = _bamReader.GetHeader();

	//set max read length, must be >=1 but smaller than 65536 (uint16_t)
	if(maxReadLength < 1)
		throw "Max read length must be at least 1 bp!";
	if(maxReadLength > 65536)
		throw "Atlas currently only supports read length up to 65536 bp. Contact us if you have longer reads / fragments";
	_maxReadLength = maxReadLength;

	//initialize chromosomes
	chromosomes.readChromosomes(&_bamHeader);

	//initialize read groups
	readGroups.fill(_bamHeader);

	//parse CIGAR
	_curCigar.clear();
	for(auto& it : _curBamAlignment.CigarData){
		_curCigar.add(it.Type, it.Length);
	}

	//get file size
	_bamReader.Jump(chromosomes.size() - 1, 0);
	BamTools::BamAlignment bamAlignment;
	_bamReader.GetNextAlignment(bamAlignment);
	_fileSize = _bamReader.tell();
	_bamReader.Rewind();
};

void TBamFile::_applyFilters(){
	//keep all?
	if(_keepAll){
		_QCFiltersPassed = true;
	} else {
		_QCFiltersPassed =  _duplicateFilter.pass(_curBamAlignment.IsDuplicate())
						 && _softClippedFilter.pass(_curCigar.lengthSoftClipped() > 0)
						 && _improperPairsFilter.pass(_curBamAlignment.IsPaired() && !_curBamAlignment.IsProperPair())
						 && _unmappedFilter.pass(!_curBamAlignment.IsMapped())
						 && _failedQCFilter.pass(_curBamAlignment.IsFailedQC())
						 && _secondaryFilter.pass(!_curBamAlignment.IsPrimaryAlignment())
						 && _supplementaryFilter.pass(_curBamAlignment.IsSupplementary())
						 && _readGroupFilter.pass(!readGroups.readGroupInUse(_curReadGroupID))
						 && _fwdStrandFilter.pass(!_curBamAlignment.IsReverseStrand())
						 && _revStrandFilter.pass(_curBamAlignment.IsReverseStrand())
						 && _firstMateFilter.pass(_curBamAlignment.IsFirstMate())
						 && _secondMateFilter.pass(_curBamAlignment.IsSecondMate())




						 && _longerThanFragmentFilter.pass()

	//standard QC
	_QCFiltersPassed = readGroups.readGroupInUse(_curReadGroupID)
					&& (_keepImproperPairs || (!_curBamAlignment.IsPaired() || _curBamAlignment.IsProperPair()))
					&& (_keepUnmappedReads || _curBamAlignment.IsMapped())
					&& (_keepFailedQC || !_curBamAlignment.IsFailedQC())
					&& (_keepSecondary || _curBamAlignment.IsPrimaryAlignment())
					&& (_keepSupplementary || !_curBamAlignment.IsSupplementary())
					&& chromosomes.inUse(_curBamAlignment.RefID)
					&& (_keepDuplicates || !_curBamAlignment.IsDuplicate())
					&& (_keepSoftClips || !_curBamAlignment.GetSoftClips(clipSizes, readPositions, genomePositions))
					&& _useStrand[_curBamAlignment.IsReverseStrand()]
					&& _useMate[_curBamAlignment.IsSecondMate()]
					&& (_curBamAlignment.MapQuality >= _minMQ && _curBamAlignment.MapQuality <= _maxMQ);

	if()


	//fragment length
	//calculate relevant fragment length
	if(_applyFragmentLengthFilter){
		uint16_t fragmentLength;
		if(_curBamAlignment.IsProperPair()){
			fragmentLength = abs(_curBamAlignment.InsertSize) + _curCigar.lengthInserted() - _curCigar.lengthDeleted();
		} else {
			fragmentLength = _curCigar.lengthSequenced();
		}
	}


	/////

	if(_curBamAlignment.IsPaired() && !_keepReadsLongerThanFragment && abs(_curBamAlignment.InsertSize) < _curBamAlignment.AlignedBases.length()){
		_QCFiltersPassed = false;
	}




	/* check if insert size is shorter than read-insertions+deletions=alignedBases.length() -> this means we are reading the adaptor sequence.
			Insert size is determined by mapping -> insertions are not in ref and should not count. If we don't add deletions, adapter at end could be sequenced but we still keep read
			(deletions in aligned bases are represented as dashes) */
			if(bamAlignment.IsPaired() && applyFragmentLengthLongerThanInsertSizeFilter && abs(bamAlignment.InsertSize) < bamAlignment.AlignedBases.length()){
				if(_updateBlacklist){
					addToBlacklist(bamAlignment, "longer than insert size (TLEN)");
				} else {
					logfile->warning("The following alignment is longer than its insert size: " + bamAlignment.Name);
				}
				filtersPassed = false;
			} else {
				//apply filters: read group in use and basic QC
				filtersPassed = applyFilters();
				if(_updateBlacklist && !filtersPassed){
					addToBlacklist(bamAlignment, "did not pass parser filters");
				}
			}

	////////
	if(_QCFiltersPassed){
		++_numAlignmentsPassedQC;
	}
};

bool TBamFile::readNextAlignment(){
	if(!_bamReader.GetNextAlignment(_curBamAlignment)){
		return false;
	}

	++_numAlignmentRead;

	//check if BAM file is sorted
	if(_curBamAlignment.Position < _previousAlignmentPos)
		throw "BAM file must be sorted by position! Alignment '" + _curBamAlignment.Name + "' is at position " + toString(_curBamAlignment.Position) + ", which is before the position of the previous alignment (" + toString(_previousAlignmentPos) + ")";
	_previousAlignmentPos = _curBamAlignment.Position;

	//check if chromosome changed
	if(_curBamAlignment.RefID != _previousAlignmentChr){
		_previousAlignmentPos = 0;
		_previousAlignmentChr = _curBamAlignment.RefID;
		_chrChanged = true;
	} else {
		_chrChanged = false;
	}

	//check read length
	if(_curBamAlignment.AlignedBases.size() > _maxReadLength)
		throw "Alignment '" +  _curBamAlignment.Name + "' is longer than the max read length " + toString(_maxReadLength) + "! Please change max read length to parse this data.";

	//store current read group ID
	std::string readGroup;
	_curBamAlignment.GetTag("RG", readGroup);
	_curReadGroupID = readGroups.find(readGroup);

	//apply filters
	_applyFilters();

	return true;
};

bool TBamFile::readNextAlignmentThatPassesFilters(){
	_QCFiltersPassed = false;
	while(!_QCFiltersPassed){
		if(!readNextAlignment()){
			return false;
		}
	}
	return true;
};

void TBamFile::fill(TAlignment & alignment){
	alignment.fill(_curBamAlignment.Name,
				   BAM::TSamFlags(_curBamAlignment.AlignmentFlag),
				   _curBamAlignment.RefID,
				   _curBamAlignment.Position,
				   _curBamAlignment.MapQuality,
				   _curCigar,
				   _curBamAlignment.MateRefID,
				   _curBamAlignment.MatePosition,
				   _curBamAlignment.InsertSize,
				   _curBamAlignment.QueryBases,
				   _curBamAlignment.Qualities,
				   _curReadGroupID);
};

bool TBamFile::readNextAlignment(TAlignment & alignment){
	if(!readNextAlignment()){
		return false;
	}

	fill(alignment);
	return true;
};

bool TBamFile::readNextAlignmentThatPassesFilters(TAlignment & alignment){
	if(!readNextAlignmentThatPassesFilters()){
		return false;
	}
	fill(alignment);
	return true;
};

bool TBamFile::jump(int chr, int position){
	_previousAlignmentPos = -1;
	return _bamReader.Jump(chr, position);
};

void TBamFile::rewind(){
	_bamReader.Rewind();
	_numAlignmentRead = 0;
	_numAlignmentsPassedQC = 0;
	_previousAlignmentPos = -1;
};

//--------------------------------------------------------
// Functions for writing
//--------------------------------------------------------
void TBamFile::openOutput(std::string filename){
	BamTools::RefVector references = _bamReader.GetReferenceData();
	if(!bamWriter.Open(filename, _bamHeader, references))
		throw "Failed to open BAM file '" + filename + "'!";
	_openForWriting = true;
};

void TBamFile::writeCurAlignment(){
	if(!_openForWriting){
		throw "BAM writer is not open!";
	}

	// write alignment
	if(!bamWriter.SaveAlignment(_curBamAlignment))
		throw "Read '" + _curBamAlignment.Name + "' could not be written!";
};

void TBamFile::writeAlignment(TAlignment & alignment, const TGenotypeMap & genoMap, const TQualityMap & qualityMap){
	if(!_openForWriting){
		throw "BAM writer is not open!"
	}

	//create bamAlignment and then write
	BamTools::BamAlignment _tmpBamAlignment;

	_tmpBamAlignment.Name = alignment.name;
	_tmpBamAlignment.AlignmentFlag = alignment.flags.asInt();
	_tmpBamAlignment.RefID = alignment.refID;
	_tmpBamAlignment.Position = alignment.position;

	_tmpBamAlignment.MapQuality = alignment.mappingQuality;
	_tmpBamAlignment.MateRefID = alignment.mateRefID;
	_tmpBamAlignment.MatePosition = alignment.matePosition;
	_tmpBamAlignment.InsertSize = alignment.insertSize_TLEN;

	//CIGAR
	for(auto& it : alignment.cigar){
		_tmpBamAlignment.CigarData.emplace_back(it.type, it.length);
	}

	//add sequences and qualities
	_tmpBamAlignment.QueryBases = alignment.getSequence(genoMap, qualityMap);
	_tmpBamAlignment.Qualities = alignment.getQualities(genoMap, qualityMap);

	//add read group information
	_tmpBamAlignment.AddTag("RG", "Z", readGroups.getName(alignment.readGroupID));

	//and now write
	if(!bamWriter.SaveAlignment(_tmpBamAlignment))
		throw "Read '" + _tmpBamAlignment.Name + "' could not be written!";

};


