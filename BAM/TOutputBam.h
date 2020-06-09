/*
 * TOutputBam.h
 *
 *  Created on: Jun 9, 2020
 *      Author: phaentu
 */

#ifndef BAM_TOUTPUTBAM_H_
#define BAM_TOUTPUTBAM_H_

#include "TBamFile.h"

namespace BAM{

//----------------------------------------------------
//TOutputBamFile
//----------------------------------------------------
struct

class TFutureAlignments{
private:
	uint32_t _refID;
	uint32_t _position;
	mutable TAlignment _alignment;

public:
	TFutureAlignments(const uint32_t RefID, const uint32_t Pos, const TAlignment & Alignment):
		_refID(RefID),
		_position(_position),
		_alignment(Alignment){};

	uint32_t position() const{ return _position; };
	TAlignment& alignment() const{ return _alignment; };

	bool operator<(const TFutureAlignments & other) const{
		//are they on the same chromosome?
		if(this->_refID < other._refID){
			return true;
		} else if(this->_refID > other._refID){
			return false;
		} else {
			//on same chromosome: check posotion
			return this->_position < other._position;
		}
	};
};


//----------------------------------------------------
//TOutputBamFile
//----------------------------------------------------
class TOutputBamFile{
	friend TBamFile;

private:
 	std::string _outputFilename;
 	BamTools::BamWriter _bamWriter;
 	bool _openForWriting;
 	TReadGroups* const _readGroups;
 	TGenotypeMap* _genoMap;
 	TQualityMap* _qualityMap;

 	std::multiset<TAlignment, std::less<>> _futureAlignments;

 	void _writeAlignment(const TAlignment & alignment);
 	void _writeAlignment(BamTools::BamAlignment & alignment);

public:
 	TOutputBamFile();
 	TOutputBamFile(const std::string filename, TBamFile & original,  TGenotypeMap* GenoMap, TQualityMap* QualityMap);
 	~TOutputBamFile();

 	void open(const std::string Filename, const TSamHeader & Header, const TChromosomes & Chromosomes, const TReadGroups & ReadGroups,  TGenotypeMap* GenoMap, TQualityMap* QualityMap);
	void open(const std::string Filename, TBamFile & original,  TGenotypeMap* GenoMap, TQualityMap* QualityMap);
	bool isOpen() const{ return _openForWriting; };
	void close(TLog* logfile);
	void closeNoIndex();
	void writeAlignment(const TAlignment & alignment);
};

}; //end namespace

#endif /* BAM_TOUTPUTBAM_H_ */
