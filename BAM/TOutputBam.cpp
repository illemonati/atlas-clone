/*
 * TOutputBam.cpp
 *
 *  Created on: Jun 9, 2020
 *      Author: phaentu
 */

#include "TOutputBam.h"

/*

namespace BAM{


//-----------------------------------------------------
//TOutputBamFile
//----------------------------------------------------
TOutputBamFile::TOutputBamFile(){
	_openForWriting = false;
	_readGroups = nullptr;
	_genoMap = nullptr;
	_qualityMap = nullptr;
};

TOutputBamFile::TOutputBamFile(const std::string Filename, const TBamFile & Original, TGenotypeMap* GenoMap, TQualityMap* QualityMap){
	_openForWriting = false;
	const TSamHeader & Header = Original.samHeader();
	open(Filename, Header, Original.chromosomes(), Original.readGroups(), GenoMap, QualityMap);
};

TOutputBamFile::~TOutputBamFile(){
	closeNoIndex();
};

void TOutputBamFile::open(const std::string Filename, const TSamHeader & Header, const TChromosomes & Chromosomes, const TReadGroups & ReadGroups, TGenotypeMap* GenoMap, TQualityMap* QualityMap){
	closeNoIndex();

	_outputFilename = Filename;
	_readGroups = &ReadGroups;
	_genoMap = GenoMap;
	_qualityMap = QualityMap;

	//construct new header /without chromosomes
	std::string header = Header.compileSamHeader(ReadGroups, Chromosomes);

	//fill bamtools chromosomes
	BamTools::RefVector ref;
	for(auto it = Chromosomes.cbegin(); it != Chromosomes.cend(); ++it){
		ref.emplace_back(it->name, it->length);
	}

	//open file for writing
	if(!_bamWriter.Open(_outputFilename, header, ref)){
		throw "Failed to open BAM file '" + _outputFilename + "'!";
	}

	_openForWriting = true;
};

void TOutputBamFile::open(const std::string Filename, const TBamFile & Original, TGenotypeMap* GenoMap, TQualityMap* QualityMap){
	open(Filename, Original.samHeader(), Original.chromosomes(), Original.readGroups(), GenoMap, QualityMap);
};

void TOutputBamFile::close(TLog* logfile){
	if(_openForWriting){
		_bamWriter.Close();

		logfile->listFlush("Creating index of BAM file '" + _outputFilename + "' ...");
		BamTools::BamReader reader;
		if(!reader.Open(_outputFilename))
			throw "Failed to open BAM file '" + _outputFilename + "' for indexing!";

		// create index for BAM file
		reader.CreateIndex(BamTools::BamIndex::STANDARD);

		//close BAM file
		reader.Close();
		logfile->done();

		_openForWriting = false;
	}
};

void TOutputBamFile::closeNoIndex(){
	if(_openForWriting){
		_bamWriter.Close();
		_openForWriting = false;
	}
};

void TOutputBamFile::_writeAlignment(BamTools::BamAlignment & alignment){
	if(!_openForWriting){
		throw "BAM writer is not open!";
	}

	//adjust qualities for printing
	_qualityMap->adjustQualitiesForWriting(alignment.Qualities);

	// write alignment
	if(!_bamWriter.SaveAlignment(alignment))
		throw "Read '" + alignment.Name + "' could not be written!";
};

void TOutputBamFile::_writeAlignment(const TAlignment & alignment){
	if(!_openForWriting){
		throw "BAM writer is not open!";
	}

	//create bamAlignment and then write
	BamTools::BamAlignment _tmpBamAlignment;

	_tmpBamAlignment.Name = alignment._name;
	_tmpBamAlignment.AlignmentFlag = alignment._flags.asInt();
	_tmpBamAlignment.RefID = alignment.refID();
	_tmpBamAlignment.Position = alignment.position();
	_tmpBamAlignment.InsertSize = alignment._insertSize_TLEN;
	_tmpBamAlignment.MapQuality = alignment._mappingQuality;

	if(alignment.isPaired()){
		_tmpBamAlignment.MateRefID = alignment.mateRefID();
		_tmpBamAlignment.MatePosition = alignment.matePosition();
	}

	//CIGAR
	for(auto& it : alignment._cigar){
		_tmpBamAlignment.CigarData.emplace_back(it.type, it.length);
	}

	//add sequences and qualities
	_tmpBamAlignment.QueryBases = alignment.sequence(*_genoMap, *_qualityMap);
	_tmpBamAlignment.Qualities = alignment.qualities(*_genoMap, *_qualityMap);

	//add read group information
	_tmpBamAlignment.AddTag("RG", "Z", _readGroups->getName(alignment._readGroupID));


	//and now write
	if(!_bamWriter.SaveAlignment(_tmpBamAlignment))
		throw "Read '" + _tmpBamAlignment.Name + "' could not be written!";
};

void TOutputBamFile::writeAlignment(const TAlignment & alignment){
	auto it = _futureAlignments.begin();
	if(it != _futureAlignments.end()){
		//write alignments BEFORE alignment to write next
		while(it != _futureAlignments.end() && *it < alignment){
			_writeAlignment(*it);
		}

		//remove written alignments
		_futureAlignments.erase(_futureAlignments.begin(), it);
	}

	//write next alignment
	_writeAlignment(alignment);
};

}; //end namespace

*/
