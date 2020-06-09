/*
 * TOutputBam.h
 *
 *  Created on: Jun 9, 2020
 *      Author: phaentu
 */

#ifndef BAM_TOUTPUTBAM_H_
#define BAM_TOUTPUTBAM_H_

#include "bamtools/api/BamReader.h"
#include "bamtools/api/BamWriter.h"

#include "TAlignment.h"
#include "TSamHeader.h"
#include <set>

namespace BAM{

class TBamFile;

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
 	TOutputBamFile(const std::string filename, const TBamFile & original, TGenotypeMap* GenoMap, TQualityMap* QualityMap);
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
